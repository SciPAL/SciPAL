#ifndef ASSEMBLE_SYSTEM_HH
#define ASSEMBLE_SYSTEM_HH

#include <step-28/timing.h>


#include <step-27/LaplaceForm.hh>

namespace step28 {

/*
// @sect5{Function: assemble_system}

// This is the function to parallelized in this project. The remainder of the
// @p DRSProblem class is in @p drs_simulation.hh .\n
// The following is the documentation given in the deal.ii docs:

// """ The following function assembles the
// linear system on the finesh level of the
// mesh. It is almost exactly the same as in
// step 6, with the exception that we don't
// eliminate hanging nodes and boundary
// values after assembling, but +while copying
// local contributions into the global
// matrix. This is not only simpler but also
// more efficient for large problems.
//
// This latter trick is something that only
// found its way into deal.II over time and
// wasn't used in the initial version of this
// tutorial program. There is, however, a
// discussion of this function in the
// introduction of step 27. """

// Most of the variables allocated at the start of the function and the
// Finite Element part of the problem is taken from the deal.ii version and
// is explained there.
template <int dim>
void DRSProblem<dim>::assemble_system (uint cycle)
{
    std::cout << __FUNCTION__ << " refactored" << std::endl;

    const QGauss<dim>  quadrature_formula(degree+1);

    LaplaceForm laplace(mapping, fe, quadrature_formula);



    // -------- BEM DATA BEGIN ---------------

    const QGauss<dim -1> regular_quad(numerics_prms.bem_quad_order);

    const unsigned char bc_id = *fem_bem.boundary_indicators.begin();

    fem_bem.solid_angles  = 2*numbers::PI;

    BemJNForm bem_form(mapping,
                       fe,
                       regular_quad,
                       fem_bem.support_points[bc_id],
                       fem_bem.dof_to_boundary_mapping);

    dealii::Vector<double> face_rhs(bem_form.dofs_per_cell);
    dealii::FullMatrix<double> face_mass_matrix(bem_form.dofs_per_cell, bem_form.dofs_per_cell);
    // -------- BEM DATA END ---------------
    //    const Coefficient<dim> coefficient;
    //    std::vector<double>    coefficient_values (n_q_points);

    typename MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();


    // Variables for CUDA computation. Note that they are allocated here
    // with the minimal sizes
    dealii::FullMatrix<double> DL_matrix_k(0, 0);
    dealii::FullMatrix<double> SL_matrix_k(0, 0);



    step28::FullMatrixAccessor local_DL_matrix_k(0,0, BemJNForm::is_col_major);
    step28::FullMatrixAccessor local_SL_matrix_k(0,0, BemJNForm::is_col_major);

    step28::CUDADriver<dim> run(bem_form.n_bem_points, bem_form.dofs_per_cell);


    std::vector<double> x_k;

    std::vector<double> q_k;
    std::vector<double> n_k;


    // If we actually use CUDA...
    if((this->arch == step28::cuda) || (this->arch == step28::both))
    {
        // ...then we need to resize the containers to be big enough to actually hold
        // out data
        local_DL_matrix_k.reinit(bem_form.n_bem_points, bem_form.dofs_per_cell);
        local_SL_matrix_k.reinit(bem_form.n_bem_points, bem_form.dofs_per_cell);


        x_k.resize(dim * bem_form.x.size());
        point_to_double_vector(x_k, bem_form.x);

        q_k.resize(dim * bem_form.q.size());
        n_k.resize(dim * bem_form.normals.size());

    }
    if(this->arch == step28::both){
        // and if we want to compare results, we need additional storage
        DL_matrix_k.reinit(bem_form.n_bem_points, bem_form.n_bem_points);
        SL_matrix_k.reinit(bem_form.n_bem_points, bem_form.n_bem_points);
    }



    // Construct a Timing object. Because the times get written to a
    // @p dealii::ConvergenceTable at the end of each call to
    // assemble_system(), the timer needs not outlive
    // this function
    step28::Timing timer(this->numerics_prms, bem_form.dofs_per_cell,
                         std::distance(cell, endc),
                         bem_form.n_bem_points,
                         bem_form.n_q_points,
                         cycle);

    dealii::Timer total_timer;
    dealii::Timer iteration_timer;

    // total_timer measures the time of the cell loop:
    total_timer.restart();

    // ------------------------- CELL LOOP BEGIN -----------------------
    double iteration_time = 0;
    uint iteration_counter = 0;

    for (; cell!=endc; ++cell)
    {
        laplace.assemble (cell);

        // Finally, distribute the local Laplace values into the global matrices.
        // Do this before doing the FEM-BEM boundary as the cell_rhs is needed for that.
        laplace.distribute_local_to_global(constraints, this->matrices(LaplacePhi), system_rhs);

        laplace.cell_rhs = 0.;
        // ------------- BEM stuff ------------
        if (cell->at_boundary() )
        {
            for (uint f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
            {
                // During the tests of the weak Dir BCs reinit has to be here
                 bem_form.reinit(cell, f, laplace.local_dof_indices);

                if(cell->face(f)->boundary_indicator() == bc_id )
                {
                  //  bem_form.reinit(cell, f, laplace.local_dof_indices);


                    // @p iteration_timer measures the time needed to construct the
                    // BEM matrices in each iteration.
                    iteration_timer.restart();

                    // If the @p arch flag is not @p cuda, then we either want to run the
                    // simple CPU implementation or  we want to both,
                    // i.e. CUDA and CPU for benchmarking purposes.
                    // At any rate, we have to run the unchanged CPU implementation.
                    if (this->arch != step28::cuda)
                    {
                        bem_form.assemble(cell, f,
                                          fem_bem.DL_matrix,
                                          fem_bem.SL_matrix);
                    }

                    // If the flag for running the assembly on the GPU using CUDA
                    // or @p both is selected:
                    if (this->arch != step28::cpu)
                    {
                        // Copy std::vector<Point<dim> > to std::vector<double>
                        point_to_double_vector(q_k, bem_form.q);
                        point_to_double_vector(n_k, bem_form.normals);

                        // Here we start the CUDA computation:\n
                        // Input:
                        //      - x_k, q_k    for G_k(i,a)
                        //      - also n_k,  for H_k(i,a)
                        //      - W(a,j)      for for SL, DL Matrix
                        //
                        // Output:
                        //      - DL_matrix_k
                        //      - SL_matrix_k
                        //
                        // Actually run the CUDA computation
                        run.assemble_bem_matrices(local_DL_matrix_k, local_SL_matrix_k,
                                                  x_k, q_k, n_k, bem_form.W);

                        // Write the calculated values into the glocal matrices
                        for (uint i=0; i < bem_form.n_bem_points; ++i)
                            for (uint j = 0; j <bem_form.dofs_per_cell; ++j)
                            {
                                uint J = bem_form.global_bc_dofs[j];
                                // loop over quad points only for valid boundary dofs
                                if (J != DoFHandler<dim>::invalid_dof_index)
                                {
                                    if (this->arch == step28::cuda)
                                    {
                                        fem_bem.DL_matrix(i,J) +=  local_DL_matrix_k(i,j);
                                        fem_bem.SL_matrix(i,J) +=  local_SL_matrix_k(i,j);
                                    }
                                // In case of @p both we have to store
                                // the CUDA-generated matrices somewhere else
                                    else {
                                        DL_matrix_k(i, J) +=  local_DL_matrix_k(i,j);
                                        SL_matrix_k(i, J) +=  local_SL_matrix_k(i,j);
                                    }
                                }
                            }
                     }

                     // Run both and compare the results, to check CUDA
                     if(this->arch == step28::both)
                     {
                         // Compare CPU and GPU results and output every differing value
                         for (uint i=0; i < bem_form.n_bem_points; ++i)
                             for (uint j = 0; j <bem_form.dofs_per_cell; ++j)
                             {
                                 uint J = bem_form.global_bc_dofs[j];

                                 // loop over quad points only for valid boundary dofs
                                 if (J != DoFHandler<dim>::invalid_dof_index)
                                 {
                                     if(std::fabs(fem_bem.DL_matrix(i,J) - DL_matrix_k(i, J) )
                                             > 1e-12)
                                         qDebug("DL(%d, %d)= %f \t\t  %f (CUDA)", i, J,
                                                fem_bem.DL_matrix(i,J), DL_matrix_k(i, J));


                                     if(std::fabs(fem_bem.SL_matrix(i,J) - SL_matrix_k(i, J) )
                                             > 1e-12)
                                         qDebug("SL(%d, %d)= %f \t\t %f (CUDA)", i, J,
                                                fem_bem.SL_matrix(i,J), SL_matrix_k(i, J));

                                 }
                             }
                    }

                    // Add up the iteration time and increment the iteration counter
                    // so that we actually know how many iterations we measured
                    iteration_time += iteration_timer.wall_time();
                    ++iteration_counter;
                }

                if(cell->face(f)->at_boundary())
                {
                bem_form.assemble_weak_Dir_BC(cell, f, *fem_bem.u_ref,
                                              face_mass_matrix,
                                              face_rhs);
                // Before we can copy the contributions of the weak form of the Dirichlet BC at the cavity surface
                                 // we have to copy it into the cell_rhs of the laplacian
                                 laplace.cell_rhs = face_rhs;
                                 laplace.cell_matrix.add(1., face_mass_matrix);

                }
            }
        }
        // For the ions we do not have
        // any other constraints than those due to hanging nodes.

        laplace.distribute_local_to_global(hanging_node_constraints,
                                           this->matrices(LaplaceIons),
                                           system_rhs_Dir_bc_in_weak_form);


    } // cell loop
    // cell loop ended, so we save the total time needed for that
    timer.averageIterationTime += iteration_time;
    timer.iteration_counter += iteration_counter;
    timer.total_time = total_timer.wall_time();

    // Last but not least, add solid angles to the diagonal of the double layer matrix
    fem_bem.alpha_plus_DL_matrix = fem_bem.DL_matrix;

    for (uint i = 0; i < fem_bem.n_bc_dofs; i++)
    {
        fem_bem.alpha_plus_DL_matrix(i,i) += fem_bem.solid_angles(i);
    }

    // compute boundary mass matrix for creating BEM contribution to FEM rhs

    {
        QGauss<dim-1> surface_quad(this->degree+1);


        fem_bem.mass_matrix.reinit(fem_bem.sparsity_pattern);
        Vector< double >  	dummy_rhs_vector(fem_bem.n_bc_dofs);
        ConstantFunction<dim>                dummy_fct(1., 1);
        fem_bem.interface[DRS::BCids::cavity()] = &dummy_fct;
        MatrixCreator::create_boundary_mass_matrix
                (mapping,
                 mg_dof_handler,
                 surface_quad,
                 fem_bem.mass_matrix,
                 fem_bem.interface,
                 dummy_rhs_vector,
                 fem_bem.dof_to_boundary_mapping
                 //     const Function< spacedim > *const  	weight = 0,
                 //     std::vector< uint > 	component_mapping = std::vector< uint >()
                 );
    }

    // Now copy dense matrices from BEM representation into sparse representation of the FEM problem.
    // The off-diagonal blocks are a bit more subtle, because the only exist as square matrices
    // but have to be padded with zeros to have the right size.
    //
    {
        // setup the matrix entries of the boundary term of the Poisson problem of the potential
       matrices.upper_off_block.reinit(sparsity_patterns.FEM_BEM_upper_off_block);

        // The sparsity pattern is (hopefully) already set. Thus we simply copy the value array.
        const unsigned int n_bulk_dofs = fem_bem.n_bulk_dofs;
        for (unsigned int r = n_bulk_dofs; r < fem_bem.n_fem_dofs; r++)
        {
            dealii::SparseMatrix<double>::const_iterator
                    col =  fem_bem.mass_matrix.begin(r - n_bulk_dofs),
                    endc = fem_bem.mass_matrix.end(r - n_bulk_dofs);

            for (; col != endc; ++col)
                matrices.upper_off_block.set(r, col->column(), col->value());
        }


        // The sparsity pattern is (hopefully) already set. Thus we simply copy the value array.
        matrices.lower_off_block.reinit(sparsity_patterns.FEM_BEM_lower_off_block);

        fem_bem.alpha_plus_DL_matrix.scatter_matrix_to(sparsity_patterns.lower_off_block_row_indices,
                                                       sparsity_patterns.lower_off_block_col_indices,
                                                       matrices.lower_off_block);
    }

    // And lastly, write out the measured times together with important
    // problem size parameters
    timer.add_values_to_table(timing_table);
} */

// @sect5{Function: assemble_system_generic}

// This is the function to be parallelized in this project. The remainder of the
// @p DRSProblem class is in @p drs_simulation.hh .\n
// The following is the documentation given in the deal.ii docs:

// """ The following function assembles the
// linear system on the finesh level of the
// mesh. It is almost exactly the same as in
// step 6, with the exception that we don't
// eliminate hanging nodes and boundary
// values after assembling, but while copying
// local contributions into the global
// matrix. This is not only simpler but also
// more efficient for large problems.
//
// This latter trick is something that only
// found its way into deal.II over time and
// wasn't used in the initial version of this
// tutorial program. There is, however, a
// discussion of this function in the
// introduction of step 27. """

// Most of the variables allocated at the start of the function and the
// Finite Element part of the problem is taken from the deal.ii version and
// is explained there.
template <int dim>
void DRSProblem<dim>::assemble_system_generic (uint cycle, step27::Architecture arch)
{
    // std::cout << __FUNCTION__ << " refactored" << std::endl;

    step27::LaplaceForm<dim> laplace(mapping, fe, *numerics_prms.fem_q_rule);



    // -------- BEM DATA BEGIN ---------------


    const unsigned char bc_id = *bem_form->boundary_indicators.begin();

    fem_bem.solid_angles  = 2*dealii::numbers::PI;



    // -------- BEM DATA END ---------------


    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();

    // Construct a Timing object. Because the times get written to a
    // @p dealii::ConvergenceTable at the end of each call to
    // assemble_system(), the timer needs not outlive
    // this function
    step28::Timing<dim> timer(this->numerics_prms,
                              bem_form->dofs_per_cell,
                         std::distance(cell, endc),
                         bem_form->n_bem_points,
                         cycle);

    dealii::Timer total_timer;
    dealii::Timer iteration_timer;

    // total_timer measures the time of the cell loop:
    total_timer.restart();

    // ------------------------- CELL LOOP BEGIN -----------------------
    double iteration_time = 0;
    uint iteration_counter = 0;

    for (; cell!=endc; ++cell)
    {
        laplace.assemble (cell );

        // Finally, distribute the local Laplace values into the global matrices.
        // Do this before doing the FEM-BEM boundary as the cell_rhs is needed for that.
        laplace.distribute_local_to_global(constraints, this->matrices(MtxInfo::LaplacePhi), system_rhs);

        laplace.cell_rhs = 0.;
        // ------------- BEM stuff ------------
        if (cell->at_boundary() )
        {

            for (uint f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
            {
                if(cell->face(f)->boundary_indicator() == bc_id )
                {
                    this->bem_form->reinit(cell, f, laplace.local_dof_indices);

                    // @p iteration_timer measures the time needed to construct the
                    // BEM matrices in each iteration.
                    iteration_timer.restart();


                    // At this point we finally have reached a face which contributes
                    // to the BEM part of the problem. Depending on the @p arch tag
                    // this function executes the assembly either on the CPU or
                    // using CUDA on the GPU.
                    this->bem_form->assemble_sl_dl(cell, f,
                                                   fem_bem.DL_matrix,
                                                   this->matrices.fm(MtxInfo::SLFullMatrix), // fem_bem.SL_matrix
                                                   arch);



                    iteration_time += iteration_timer.wall_time();
                    ++iteration_counter;
                }
            }
        }

        // laplace.get_function_values(this->old_solution.neutrals, old_sol_values);
        // laplace.assemble_tdep_rhs_contrib(timestep() , old_sol_values, tdep_cell_rhs);

        // For the ions we do not have
        // any other constraints than those due to hanging nodes.

        // laplace.cell_rhs = pure_vN_form.cell_rhs; // laplace.cell_rhs += tdep_cell_rhs;
        // In contrast to the potential we need a diffusion operator for the ions
        laplace.cell_matrix.add(1./this->numerics_prms.pseudo_time_step, laplace.cell_mass_matrix);
        laplace.distribute_local_to_global(hanging_node_constraints,
                                           this->matrices(MtxInfo::LaplaceIons) //,
                                        //  system_rhs_pure_vN_bc
                                           // system_rhs_Dir_bc_in_weak_form
                                           );

        // Additionally, we need the mass matrix itself.
        // Therefore, we have to copy it into the @p cell_matrix
        // in order to build the global mass matrix.
        laplace.cell_matrix = laplace.cell_mass_matrix;
        laplace.distribute_local_to_global(hanging_node_constraints,
                                                   this->matrices(MtxInfo::BulkMass) // ,
                                                 //  system_rhs_pure_vN_bc
                                                   // system_rhs_Dir_bc_in_weak_form
                                                   );


    } // cell loop


    // cell loop ended, so we save the total time needed for that
    timer.average_iteration_time += iteration_time;
    timer.iteration_counter += iteration_counter;
    timer.total_time = total_timer.wall_time();


    // Last but not least, add solid angles to the diagonal of the double layer matrix
    fem_bem.alpha_plus_DL_matrix = fem_bem.DL_matrix;

    for (uint i = 0; i < fem_bem.n_bc_dofs; i++)
        fem_bem.alpha_plus_DL_matrix(i,i) += fem_bem.solid_angles(i);


    assemble_boundary_mass_matrices();

    if (!this->drs_sol_comps_fem_part_all_vN.empty())
        reassemble_nonlinearity_cpu ();

    // And lastly, write out the measured times together with important
    // problem size parameters
    timer.add_values_to_table(timing_table);
}


template <int dim>
void DRSProblem<dim>::reassemble_nonlinearity_cpu ()
{
    // std::cout << __FUNCTION__ << " refactored" << std::endl;


    const dealii::QGauss<dim>  quadrature_formula(degree+1);

    DriftForm drift(mapping, fe, quadrature_formula,
                    this->old_solution(step27::Cations),
                    this->old_solution(step27::Anions),
                    this->old_solution(step27::Phi) );


    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();

    // reset global matrices
    this->matrices(MtxInfo::DriftCations) = 0.;
    this->matrices(MtxInfo::DriftAnions) = 0.;

    // ------------------------- CELL LOOP BEGIN -----------------------

    for (; cell!=endc; ++cell)
    {
        drift.assemble (cell );

        // Finally, distribute the local Laplace values into the global matrices.
        // Do this before doing the FEM-BEM boundary as the cell_rhs is needed for that.
        drift.distribute_local_to_global(hanging_node_constraints,
                                         this->matrices(MtxInfo::DriftCations),
                                         this->matrices(MtxInfo::DriftAnions),
                                         this->matrices(MtxInfo::DriftPot));
    } // cell loop

}


// TODO:
template <int dim>
void DRSProblem<dim>::assemble_boundary_mass_matrices()
{

    // compute boundary mass matrix for creating BEM contribution to FEM rhs
    dealii::QGauss<dim-1> surface_quad(this->degree+1);

    dealii::Vector< double >  	dummy_rhs_vector(fem_bem.n_bc_dofs);

    dealii::ConstantFunction<dim>                dummy_fct(1., 1);

    {
        this->matrices.reinit(MtxInfo::FEMBEMMass);

        fem_bem.interface[DRS::BCids::cavity()] = &dummy_fct;

        dealii::MatrixCreator::create_boundary_mass_matrix
                (mapping,
                 mg_dof_handler,
                 surface_quad,
                  this->matrices(MtxInfo::FEMBEMMass),
                 fem_bem.interface,
                 dummy_rhs_vector,
                 bem_form->dof_to_boundary_mapping
                 //     const Function< spacedim > *const  	weight = 0,
                 //     std::vector< uint > 	component_mapping = std::vector< uint >()
                 );
    }

    // Now copy dense matrices from BEM representation into sparse representation of the FEM problem.
    // The off-diagonal blocks are a bit more subtle, because the only exist as square matrices
    // but have to be padded with zeros to have the right size.
    //
    {
        // setup the matrix entries of the boundary term of the Poisson problem of the potential
        matrices(MtxInfo::FEMBEMUpperOffBlock);

        // The sparsity pattern is (hopefully) already set. Thus we simply copy the value array.
        const unsigned int n_bulk_dofs = fem_bem.n_bulk_dofs;
        for (unsigned int r = n_bulk_dofs; r < fem_bem.n_fem_dofs; r++)
        {
            dealii::SparseMatrix<double>::const_iterator
                    col =   this->matrices(MtxInfo::FEMBEMMass).begin(r - n_bulk_dofs),
                    endc =  this->matrices(MtxInfo::FEMBEMMass).end(r - n_bulk_dofs);

            for (; col != endc; ++col)
                matrices(MtxInfo::FEMBEMUpperOffBlock).set(r, col->column(), col->value());
        }


        // The sparsity pattern is (hopefully) already set. Thus we simply copy the value array.
        matrices(MtxInfo::BEMFEMLowerOffBlock);

        fem_bem.alpha_plus_DL_matrix.scatter_matrix_to(sparsity_patterns.lower_off_block_row_indices,
                                                       sparsity_patterns.lower_off_block_col_indices,
                                                       matrices(MtxInfo::BEMFEMLowerOffBlock));
    }


    // Boundary mass matrices for electrodes



    {

        dealii::SparseMatrix<double> tmp_cathode(sparsity_patterns(PatternInfo::BEMCathodeBoundaryMass)/*.cathode_mass_bc_dofs*/);

        dealii::types::boundary_id bc_id = DRS::BCids::cathode();

        std::vector<unsigned int> & cathode_dof_2_bc_map
                = sparsity_patterns.dof_2_bc_map[bc_id];

        typename dealii::FunctionMap<dim>::type cathode_indicator;
        cathode_indicator[bc_id] = &dummy_fct;

        dummy_rhs_vector.reinit(mg_dof_handler.n_boundary_dofs(cathode_indicator) );

        // Quote from the deal.II doc: "The size of the matrix is equal to the number of degrees of freedom that have support on the boundary,
        // i.e. it is not a matrix on all degrees of freedom, but only a subset."
        dealii::MatrixCreator::create_boundary_mass_matrix
                       (mapping,
                        mg_dof_handler,
                        surface_quad,
                        tmp_cathode,
                        cathode_indicator,
                        dummy_rhs_vector,
                        cathode_dof_2_bc_map
                        //     const Function< spacedim > *const  	weight = 0,
                        //     std::vector< uint > 	component_mapping = std::vector< uint >()
                        );

        std::map<uint, uint> & electrode_bc_local_2_global
                          = sparsity_patterns.electrode_bc_local_2_global_map[bc_id];

        std::map<uint, uint>::const_iterator bc_dof = electrode_bc_local_2_global.begin(),
                end_bc_dofs = electrode_bc_local_2_global.end();

        for (; bc_dof != end_bc_dofs; ++bc_dof)
        {
            unsigned int bc_r = bc_dof->first;
            unsigned int global_r = bc_dof->second;

            if (global_r != dealii::DoFHandler<dim>::invalid_dof_index)
            {
                dealii::SparseMatrix<double>::const_iterator
                        col =  tmp_cathode.begin(bc_r),
                        endc = tmp_cathode.end(bc_r);

                for (; col != endc; ++col) {
                    std::map<uint, uint>::const_iterator
                            c = electrode_bc_local_2_global.find(col->column());
                    Assert(c != end_bc_dofs,
                           dealii::ExcMessage("Illegal bc dof on electrode boundary"));
                    unsigned int global_col = c->second;
                    matrices(MtxInfo::CathodeBoundaryMass).set(global_r, global_col, col->value());
                }
            }
        }
    }



    {

        dealii::SparseMatrix<double> tmp_anode(sparsity_patterns(PatternInfo::BEMAnodeBoundaryMass)/*.anode_mass_bc_dofs*/);

        dealii::types::boundary_id bc_id = DRS::BCids::anode();

        std::vector<unsigned int> & anode_dof_2_bc_map
                = sparsity_patterns.dof_2_bc_map[bc_id];

        typename dealii::FunctionMap<dim>::type anode_indicator;
        anode_indicator[bc_id] = &dummy_fct;

        dummy_rhs_vector.reinit(mg_dof_handler.n_boundary_dofs(anode_indicator) );

        dealii::MatrixCreator::create_boundary_mass_matrix
                       (mapping,
                        mg_dof_handler,
                        surface_quad,
                        tmp_anode,
                        anode_indicator,
                        dummy_rhs_vector,
                        anode_dof_2_bc_map
                        //     const Function< spacedim > *const  	weight = 0,
                        //     std::vector< uint > 	component_mapping = std::vector< uint >()
                        );

        std::map<uint, uint> & electrode_bc_local_2_global
                          = sparsity_patterns.electrode_bc_local_2_global_map[bc_id];

        std::map<uint, uint>::const_iterator bc_dof = electrode_bc_local_2_global.begin(),
                end_bc_dofs = electrode_bc_local_2_global.end();

        for (; bc_dof != end_bc_dofs; ++bc_dof)
        {
            unsigned int bc_r = bc_dof->first;
            unsigned int global_r = bc_dof->second;

            if (global_r != dealii::DoFHandler<dim>::invalid_dof_index)
            {
                dealii::SparseMatrix<double>::const_iterator
                        col =  tmp_anode.begin(bc_r),
                        endc = tmp_anode.end(bc_r);

                for (; col != endc; ++col) {
                    std::map<uint, uint>::const_iterator
                            c = electrode_bc_local_2_global.find(col->column());
                    Assert(c != end_bc_dofs,
                           dealii::ExcMessage("Illegal bc dof on electrode boundary"));
                    unsigned int global_col = c->second;
                    matrices(MtxInfo::AnodeBoundaryMass).set(global_r, global_col, col->value());
                }
            }
        }
    }

}



// @sect4{DRSProblem::assemble_multigrid}

// The next function is the one that builds
// the linear operators (matrices) that
// define the multigrid method on each level
// of the mesh. The integration core is the
// same as above, but the loop below will go
// over all existing cells instead of just
// the active ones, and the results must be
// entered into the correct matrix. Note also
// that since we only do multi-level
// preconditioning, no right-hand side needs
// to be assembled here.
//
// Before we go there, however, we have to
// take care of a significant amount of book
// keeping:
template <int dim>
void DRSProblem<dim>::assemble_multigrid ()
{
    std::cout << __FUNCTION__ << " begin" << std::endl;

    // Objects for computing a cell's contribution to the FE part of the problem.
    dealii::QGauss<dim>  quadrature_formula(1+degree);

    step27::LaplaceForm<dim> laplace(mapping, fe, quadrature_formula);


    // -------- BEM DATA BEGIN ---------------

     // Objects for computing a cell's contribution to the BIE part of the problem.
    const dealii::QGauss<dim-1>  regular_quad(fe.degree+1);

    const unsigned char bc_id = *bem_form->boundary_indicators.begin();

//    step27::BemCollocationForm<dim, step27::FaceIntTraits<dim> > bem_form(mapping,
//                                fe,
//                                regular_quad,
//                                bem_form->support_points[bc_id],
//                                fem_bem.dof_to_boundary_mapping);

//    dealii::Vector<double> face_rhs(bem_form.dofs_per_cell);
//    dealii::FullMatrix<double> face_mass_matrix(bem_form.dofs_per_cell, bem_form.dofs_per_cell);

    // We have to backup the laplacian for the computation of the interface matrices.
    dealii::FullMatrix<double> cell_matrix_backup(laplace.dofs_per_cell,
                                                  laplace.dofs_per_cell),
            cell_laplace_only(laplace.dofs_per_cell,
                              laplace.dofs_per_cell);

    // -------- BEM DATA END ---------------

    // Next a few things that are specific to building the multigrid data structures
    // (since we only need them in the current function, rather than also elsewhere, we
    // build them here instead of the <code>setup_system</code>
    // function). Some of the following may be a bit obscure if you're not familiar
    // with the algorithm actually implemented in deal.II to support multilevel
    // algorithms on adaptive meshes; if some of the things below seem strange, take a
    // look at the @ref mg_paper.
    //
    // Our first job is to identify those degrees of freedom on each level that
    // are located on interfaces between adaptively refined levels, and those
    // that lie on the interface but also on the exterior boundary of the domain. As
    // in many other parts of the library, we do this by using boolean masks,
    // i.e. vectors of booleans each element of which indicates whether the
    // corresponding degree of freedom index is an interface DoF or not. The <code>MGConstraints</code>
    // already computed the information for us when we called initialize in <code>setup_system()</code>.
    std::vector<std::vector<bool> > Phi_interface_dofs
            = mg_Phi.constrained_dofs.get_refinement_edge_indices ();

    std::vector<std::vector<bool> > Phi_boundary_interface_dofs
            = mg_Phi.constrained_dofs.get_refinement_edge_boundary_indices ();



    std::vector<std::vector<bool> > Ions_interface_dofs
               = mg_Ions.constrained_dofs.get_refinement_edge_indices ();

       std::vector<std::vector<bool> > Ions_boundary_interface_dofs
               = mg_Ions.constrained_dofs.get_refinement_edge_boundary_indices ();


    // The indices just identified will later be used to decide where the assembled value
    // has to be added into on each level.
    // On the other hand, we also have to impose zero boundary conditions on the external boundary of
    // each level. But this the <code>MGConstraints</code> knows it. So we simply ask for them by calling
    // <code>get_boundary_indices ()</code>.
    // The third step is to construct constraints on all those degrees of
    // freedom: their value should be zero after each application of the level
    // operators. To this end, we construct ConstraintMatrix objects for each level,
    // and add to each of these constraints for each degree of freedom. Due to the way
    // the ConstraintMatrix stores its data, the function to add a constraint on a
    // single degree of freedom and force it to be zero is called
    // Constraintmatrix::add_line(); doing so for several degrees of freedom at once
    // can be done using Constraintmatrix::add_lines():
    std::vector<dealii::ConstraintMatrix> Phi_boundary_constraints (triangulation.n_levels());
    std::vector<dealii::ConstraintMatrix> Ions_boundary_constraints (triangulation.n_levels());

    std::vector<dealii::ConstraintMatrix> Phi_boundary_interface_constraints (triangulation.n_levels());
    std::vector<dealii::ConstraintMatrix> Ions_boundary_interface_constraints (triangulation.n_levels());

    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
        Phi_boundary_constraints[level].add_lines (Phi_interface_dofs[level]);

        Phi_boundary_constraints[level].add_lines (mg_Phi.constrained_dofs.get_boundary_indices()[level]);
        Phi_boundary_constraints[level].close ();

        // The ion species do not have to obey any Dirichlet boundary conditions
        // but rather inhomogeneous von Neumann conditions. Therefore, we can omit adding the
        // boundary dofs to the constraints object.
        Ions_boundary_constraints[level].add_lines (Ions_interface_dofs[level]);
        Ions_boundary_constraints[level].close ();


        Phi_boundary_interface_constraints[level].add_lines (Phi_boundary_interface_dofs[level]);
        Phi_boundary_interface_constraints[level].close ();

        // There are no boundary interface constraints for the ions
        Ions_boundary_interface_constraints[level].close();
    }

    // Now that we're done with most of our preliminaries, let's start the
    // integration loop. It looks mostly like the loop in
    // <code>assemble_system</code>, with two exceptions: (i) we don't need a right
    // hand side, and more significantly (ii) we don't just loop over all active cells,
    // but in fact all cells, active or not. Consequently, the correct iterator
    // to use is MGDoFHandler::cell_iterator rather than
    // MGDoFHandler::active_cell_iterator. Let's go about it:
    typename dealii::MGDoFHandler<dim>::cell_iterator
            cell = mg_dof_handler.begin(),
            endc = mg_dof_handler.end();



    for (; cell!=endc; ++cell)
    {
        unsigned int level = cell->level();

        laplace.assemble(cell );
        // The rest of the assembly is again slightly different. This starts with
        // a gotcha that is easily forgotten:
        // The indices of global degrees of freedom we want here are the ones
        // for current level, not for the global matrix. We therefore need the
        // function MGDoFAccessor::get_mg_dof_indices,
        // not MGDoFAccessor::get_dof_indices as used in the assembly of the
        // global system:
        cell->get_mg_dof_indices (laplace.local_dof_indices);

        // Next, we need to copy local contributions into the level
        // objects. We can do this in the same way as in the global assembly, using
        // a constraint object that takes care of constrained degrees (which here
        // are only boundary nodes, as the individual levels have no hanging
        // node constraints). Note that the <code>boundary_constraints</code>
        // object makes sure that the level matrices contain no contributions
        // from degrees of freedom at the interface between cells of different
        // refinement level.
        laplace.distribute_local_to_global(Phi_boundary_constraints[level],
                                           mg_Phi.matrices[level]);

        // After the bulk Laplacian we can take care of the boundary mass matrix
        // needed for the weak form of the Dirichlet BC.

        if (cell->at_boundary() )
         {
            // First backup what we just have computed
            cell_matrix_backup = laplace.cell_matrix;

            if (false)
            for (uint f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
            {
                if(cell->face(f)->boundary_indicator() == bc_id )
                {
                    // bem_form.reinit(cell, f, laplace.local_dof_indices);
                    //
                    // If the @p arch flag is not @p cuda, then we either want to run the
                    // simple CPU implementation or  we want to both,
                    // i.e. CUDA and CPU for benchmarking purposes.
                    // At any rate, we have to run the unchanged CPU implementation.
                    if (this->arch != step27::cuda)
                    {
#ifdef USE_BEM_MG_PC
                        bem_form.assemble_SL_mg_preconditioner(cell, f, fem_bem.SL_mg_matrices[level]);
#endif
                    }
                }

#ifdef USE_WEAK_DIR_BC
                if(cell->face(f)->at_boundary())
                {
                    bem_form.assemble_weak_Dir_BC(cell, f, *fem_bem.u_ref,
                                                  face_mass_matrix,
                                                  face_rhs);
                    // Before we can copy the contributions of the weak form of the Dirichlet BC at the cavity surface
                    // we have to add it to the cell_matrix of the laplacian
                    laplace.cell_matrix.add(1., face_mass_matrix);
                }
#endif
            }
        }


        // First backup what we just have computed
        cell_matrix_backup = laplace.cell_matrix;
        cell_laplace_only = laplace.cell_matrix;

        // DIFFUSION OP
        laplace.cell_matrix.add(1./this->numerics_prms.pseudo_time_step, laplace.cell_mass_matrix);

        laplace.distribute_local_to_global(Ions_boundary_constraints[level],
                                           mg_Ions.matrices[level]);



        // The next step is again slightly more obscure (but explained in the @ref
        // mg_paper): We need the remainder of the operator that we just copied
        // into the <code>mg_matrices</code> object, namely the part on the
        // interface between cells at the current level and cells one level
        // coarser. This matrix exists in two directions: for interior DoFs (index
        // $i$) of the current level to those sitting on the interface (index
        // $j$), and the other way around. Of course, since we have a symmetric
        // operator, one of these matrices is the transpose of the other.
        //
        // The way we assemble these matrices is as follows: since the are formed
        // from parts of the local contributions, we first delete all
        // those parts of the local contributions that we are not
        // interested in, namely all those elements of the local matrix for
        // which not $i$ is an interface DoF and $j$ is not. The result is one of
        // the two matrices that we are interested in, and we then copy it
        // into the <code>mg_interface_matrices</code>
        // object. The <code>boundary_interface_constraints</code>
        // object at the same time makes sure that we delete contributions from
        // all degrees of freedom that are not only on the interface but also on
        // the external boundary of the domain.
        //
        // The last part to remember is how to get the other matrix. Since it is
        // only the transpose, we will later (in the <code>solve()</code>
        // function) be able to just pass the transpose matrix where necessary.
        for (unsigned int i=0; i<laplace.dofs_per_cell; ++i)
            for (unsigned int j=0; j<laplace.dofs_per_cell; ++j)
                if( !(Phi_interface_dofs[level][laplace.local_dof_indices[i]]==true &&
                      Phi_interface_dofs[level][laplace.local_dof_indices[j]]==false))
                {
                                  laplace.cell_matrix(i,j) = 0;
                                  cell_laplace_only(i,j) = 0;
}
        // DIFFUSION OP
        laplace.distribute_local_to_global(Ions_boundary_interface_constraints[level],
                                          mg_Ions.interface_matrices[level]);

        // RESTORE
        laplace.cell_matrix = cell_laplace_only;


#ifdef USE_WEAK_DIR_BC
        // If necessary, restore original Laplacian
        if (cell->at_boundary() )
        {
            laplace.cell_matrix = cell_matrix_backup;
            for (unsigned int i=0; i<laplace.dofs_per_cell; ++i)
                for (unsigned int j=0; j<laplace.dofs_per_cell; ++j)
                    if( !(Phi_interface_dofs[level][laplace.local_dof_indices[i]]==true &&
                          Phi_interface_dofs[level][laplace.local_dof_indices[j]]==false))
                        laplace.cell_matrix(i,j) = 0;

        }
#endif

        laplace.distribute_local_to_global(Phi_boundary_interface_constraints[level],
                                           mg_Phi.interface_matrices[level]);

    }


     std::cout << __FUNCTION__ << " end" << std::endl;
}

// Test: $1/2 u(x_i) = +V t_P - K \Phi_P$
          //
          // $t_P = d_{n_P} \Phi_P = - d_n \Phi_P $
          //
          // For BIEs we need the normal in the opposite direction
          // $V (-t_P)$.
template <int dim>
void DRSProblem<dim>::assemble_bem_unit_test()
{
    const dealii::QGauss<dim-1>  regular_quad(fe.degree+1);

    const unsigned char bc_id = *fem_bem.boundary_indicators.begin();

    step27::BemCollocationForm<dim, step27::FaceIntTraits<dim> > bem_form(mapping,
                                fe,
                                regular_quad,
                                fem_bem.support_points[bc_id],
                                fem_bem.dof_to_boundary_mapping);

    // For testing the integrals
    std::vector<double>          u_values    (bem_form.n_q_points);
    std::vector<double>          u_h_values  (bem_form.n_q_points);
    std::vector<dealii::Tensor<1, dim> > u_gradients (bem_form.n_q_points);

    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();

    for (; cell!=endc; ++cell)
    {
        // ------------- BEM stuff ------------
        if (cell->at_boundary() )
            // Test: $1/2 u(x_i) = +V t_P - K \Phi_P$
            //
            // $t_P = d_{n_P} \Phi_P = - d_n \Phi_P $
            //
            // For BIEs we need the normal in the opposite direction
            // $V (-t_P)$.
            for (uint f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
                if(cell->face(f)->boundary_indicator() == bc_id )
                {
                    bem_form.reinit(cell, f);

                    fem_bem.u_ref->value_list (bem_form.q, u_values);

                    bem_form.get_function_values(fem_bem.solution_ref, u_h_values);

                    fem_bem.u_ref->gradient_list(bem_form.q, u_gradients);

                    for (uint i=0; i < bem_form.n_bem_points; ++i)

                        for (uint a=0; a < bem_form.n_q_points; ++a)
                        {
                            dealii::Point<dim> R = bem_form.x[i] - bem_form.q[a];

                            double G_ia = step27::LaplaceKernel<dim>::single_layer(R);

                            fem_bem.rhs_for_dirichlet_test(i) +=
                                    + G_ia  * (-bem_form.normals[a] * u_gradients[a]) * bem_form.JxW[a];
                        }
                }
    }
}


} // namespace Step16 END


#endif // ASSEMBLE_SYSTEM_HH
