#ifndef ASSEMBLE_SYSTEM_HH
#define ASSEMBLE_SYSTEM_HH

#include <step-27/timing.h>
#include <step-27/UnitTestProblemData.hh>


// @sect5{Function: assemble_system_generic}

// This is the function to be parallelized in this project. The remainder of the
// @p FEMBEMPoissonProblem class is in @p drs_simulation.hh .\n
// The following is the documentation given in the deal.ii docs:

// """ The following function assembles the
// linear system on the finest level of the
// mesh. It is almost exactly the same as in
// deal.step 6, with the exception that we don't
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
// introduction of deal.II's step 27. """

// Most of the variables allocated at the start of the function and the
// finite element part of the problem is taken from the deal.ii version and
// is explained there.
template <int dim>
void step27::FEMBEMPoissonProblem<dim>::assemble_system_generic (uint cycle, step27::Architecture arch)
{
    const unsigned char bc_id = *bem_form->boundary_indicators.begin();

    fem_bem.solid_angles  = 2*dealii::numbers::PI;


    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();

    // Construct a Timing object. Because the times get written to a
    // @p dealii::ConvergenceTable at the end of each call to
    // assemble_system(), the timer needs not outlive
    // this function
    step27::Timing<dim> timer(this->numerics_prms,
                              bem_form->dofs_per_cell,
                         std::distance(cell, endc),
                         bem_form->n_bem_points,
                         cycle);

    dealii::Timer total_timer;
    dealii::Timer iteration_timer;

    // totalTimer measures the time of the cell loop:
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
        // If a cell is at a boundary, we have check which one it is and if it is part of the dielectric interface
        // we have to assemble the contribution to the BEM matrices.
        if (cell->at_boundary() )
        {

            // The different faces of a cell can be parts of different boundaries.
            for (uint f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
            {
                if(cell->face(f)->boundary_indicator() == bc_id )
                {

                    this->bem_form->reinit(cell, f, laplace.local_dof_indices);

                    // @p iteration_timer measures the time needed to construct the
                    // BEM matrices in each iteration.
                    iteration_timer.restart();

                    // TO DO: pass pointers to the local assembers and put architecture dependence into the
                    // form objects.
                    // if(arch == step27::cpu)

                    // Run the unchanged CPU implementation.
                    bem_form->assemble_sl_dl(cell, f,
                                            fem_bem.DL_matrix,
                                             this->matrices.fm(MtxInfo::SLFullMatrix),
                                             arch
                                           // fem_bem.SL_matrix
                                             );


                    iteration_time += iteration_timer.wall_time();
                    ++iteration_counter;
                }
            }
        }

    } // cell loop


    // cell loop ended, so we save the total time needed for that
    timer.average_time_per_iteration += iteration_time;
    timer.n_iterations += iteration_counter;
    timer.total_time = total_timer.wall_time();


    // Last but not least, add solid angles to the diagonal of the double layer matrix
    fem_bem.alpha_plus_DL_matrix = fem_bem.DL_matrix;

    for (uint i = 0; i < fem_bem.n_bc_dofs; i++)
        fem_bem.alpha_plus_DL_matrix(i,i) += fem_bem.solid_angles(i);


    assemble_boundary_mass_matrices();



    // And lastly, write out the measured times together with important
    // problem size parameters
    timer.add_values_to_table(timing_table);
}




// TODO:
template <int dim>
void step27::FEMBEMPoissonProblem<dim>::assemble_boundary_mass_matrices()
{

    // compute boundary mass matrix for creating BEM contribution to FEM rhs
    dealii::QGauss<dim-1> surface_quad(this->degree+1);

    dealii::Vector< double >  	dummy_rhs_vector(fem_bem.n_bc_dofs);

    dealii::ConstantFunction<dim>                dummy_fct(1., 1);

    {
        this->matrices.reinit(MtxInfo::FEMBEMMass);
       // fem_bem.mass_matrix.reinit(fem_bem.sparsity_pattern);

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

    // Now copy dense matrices from the BEM representation into the sparse representation of the FEM problem.
    // The off-diagonal blocks are a bit more subtle, because they only exist as square matrices
    // but have to be padded with zeros to have the right size.
    //
    {
        // setup the matrix entries of the boundary term of the Poisson problem of the potential
        // TODO: why pass pointer to pattern_set to mtx_set?
        matrices.reinit(MtxInfo::FEMBEMUpperOffBlock);


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
        matrices.reinit(MtxInfo::BEMFEMLowerOffBlock); //.reinit(sparsity_patterns.FEM_BEM_lower_off_block);

        // For the lower off-diagonal block we can use a function already provided by deal.II.
        fem_bem.alpha_plus_DL_matrix.scatter_matrix_to(sparsity_patterns.lower_off_block_row_indices,
                                                       sparsity_patterns.lower_off_block_col_indices,
                                                       matrices(MtxInfo::BEMFEMLowerOffBlock));
    }


    // Boundary mass matrices for electrodes

#ifdef USE_EL_MASS_MAT

    // Assembly
    {

        dealii::SparseMatrix<double> tmp_cathode(sparsity_patterns.cathode_mass_bc_dofs);

        dealii::types::boundary_id bc_id = DRS::BCids::cathode();

        std::vector<unsigned int> & cathode_dof_2_bc_map
                = sparsity_patterns.dof_2_bc_map[bc_id];

        typename dealii::FunctionMap<dim>::type cathode_indicator;
        cathode_indicator[bc_id] = &dummy_fct;

        dummy_rhs_vector.reinit(mg_dof_handler.n_boundary_dofs(cathode_indicator) );

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

        dealii::SparseMatrix<double> tmp_anode(sparsity_patterns.anode_mass_bc_dofs);

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
#endif
}



// @sect4{FEMBEMPoissonProblem::assemble_multigrid}

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
void step27::FEMBEMPoissonProblem<dim>::assemble_multigrid ()
{
    std::cout << __FUNCTION__ << " begin" << std::endl;

    // -------- BEM DATA BEGIN ---------------
    const unsigned char bc_id = *bem_form->boundary_indicators.begin();
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
                    // bem_form->reinit(cell, f, laplace.local_dof_indices);
                    //
                    // If the @p arch flag is not @p cuda, then we either want to run the
                    // simple CPU implementation or  we want to both,
                    // i.e. CUDA and CPU for benchmarking purposes.
                    // At any rate, we have to run the unchanged CPU implementation.
#ifdef USE_BEM_MG_PC
                    {
                        bem_form->assemble_SL_mg_preconditioner(cell, f, fem_bem.SL_mg_matrices[level]);
                    }
#endif
                }
            }
        }


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
void step27::FEMBEMPoissonProblem<dim>::assemble_bem_unit_test()
{
    const dealii::QGauss<dim-1>  regular_quad(fe.degree+1);

    const unsigned char bc_id = *fem_bem.boundary_indicators.begin();

    BemCollocationForm<dim, FaceIntTraits<dim> > bem_form(mapping,
                                fe,
                                regular_quad,
                                fem_bem.support_points[bc_id],
                                fem_bem.dof_to_boundary_mapping);

    // For testing the integrals
    std::vector<double>          u_values    (bem_form->n_q_points);
    std::vector<double>          u_h_values  (bem_form->n_q_points);
    std::vector<dealii::Tensor<1, dim> > u_gradients (bem_form->n_q_points);

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
                    bem_form->reinit(cell, f);

                    fem_bem.u_ref->value_list (bem_form->q, u_values);

                    bem_form->get_function_values(fem_bem.solution_ref, u_h_values);

                    fem_bem.u_ref->gradient_list(bem_form->q, u_gradients);

                    for (uint i=0; i < bem_form->n_bem_points; ++i)

                        for (uint a=0; a < bem_form->n_q_points; ++a)
                        {
                            dealii::Point<dim> R = bem_form->x[i] - bem_form->q[a];

                            double G_ia = LaplaceKernel<dim>::single_layer(R);

                            fem_bem.rhs_for_dirichlet_test(i) +=
                                    + G_ia  * (-bem_form->normals[a] * u_gradients[a]) * bem_form->JxW[a];
                        }
                }
    }
}

#endif // ASSEMBLE_SYSTEM_HH
