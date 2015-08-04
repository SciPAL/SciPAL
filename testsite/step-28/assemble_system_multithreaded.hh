#ifndef ASSEMBLE_SYSTEM_MULTITHREADED_HH
#define ASSEMBLE_SYSTEM_MULTITHREADED_HH




namespace Step16 {

// @sect4{Class: Step16::CellWorker}
//
// This class is used to partition the cell loop onto multiple threads and
// calculate the contributions in parallel. To this end, this class inherits from
// QThread and overwrites the @p run() function. This way, if @p start() is called on
// a @p CellWorker object, a new thread is started and @p run() is executed in that thread.
//\n\n
// It is put into the Step16 namespace (the namespace of the drs simulation), because it is
// only used in Step16::DRSProblem functions.

template<int dim>
class CellWorker : public QThread {
public:

    CellWorker(typename dealii::MGDoFHandler<dim>::active_cell_iterator _start_cell,
               typename dealii::MGDoFHandler<dim>::active_cell_iterator _end_cell,
               const step27::NumericsParams<dim> &_numerics_prms,
               typename Step16::FemBemData<dim> &_fem_bem,
               dealii::ConstraintMatrix &_constraints, dealii::SparseMatrix<double> &_system_matrix,
               dealii::Vector<double> &_system_rhs, const step28::Architecture _arch,
               const typename step27::PhysicalParams &_phys_data,
               double &_iterationTime, int &_iteration_counter);


    void run();


private:
    typename dealii::MGDoFHandler<dim>::active_cell_iterator start_cell;
    typename dealii::MGDoFHandler<dim>::active_cell_iterator end_cell;



    const step27::NumericsParams<dim> &numerics_prms;

    typename Step16::FemBemData<dim> &fem_bem;
    dealii::ConstraintMatrix     &constraints;

    dealii::SparseMatrix<double> &system_matrix;
    dealii::Vector<double>       &system_rhs;

    step28::Architecture arch;

    const typename step27::PhysicalParams &phys_data;

    static QMutex mutex;

    double &total_iteration_time;
    int    &total_iteration_counter;

};



}// namespace Step16




namespace Step16 {


// @sect5{Function: assemble_system_multithreaded}
//
// Multi threaded assembly of the system.
// Since the calculations on each cell are independent, these are
// easily parallelizable on the CPU side as well, just by starting multiple
// threads and giving each thread some of the cells to work on.
//
template <int dim>
void DRSProblem<dim>::assemble_system_multithreaded (uint cycle)
{
    std::cout << __FUNCTION__ << std::endl;


    fem_bem.solid_angles  = 2*dealii::numbers::PI;


    typename dealii::MGDoFHandler<dim>::active_cell_iterator cell_begin = mg_dof_handler.begin_active();
    typename dealii::MGDoFHandler<dim>::active_cell_iterator cell_end   = mg_dof_handler.end();

    uint n_cells = std::distance(cell_begin, cell_end);

    // The bem_values object constructed here is only needed to get the required
    // Parameters to construct the step28::Timing object.
    const dealii::QGauss<dim -1> regular_quad(numerics_prms.bem_quad_order_1);
    dealii::FEFaceValues<dim> bem_values (mapping,
                                  fe,
                                  regular_quad, dealii::update_values);


    // Construct a Timing object. Because the times get written to a
    // deal.II convergence table owned be the DRSProblem
    // at the end of each call to assemble_system(), the timer needs not outlive
    // this function
    step28::Timing<dim> timer(this->numerics_prms, bem_values.dofs_per_cell,
                         n_cells,
                         fem_bem.support_points[*fem_bem.boundary_indicators.begin()].size(),
                         cycle);

    double total_iteration_time = 0;
    int total_iteration_counter = 0;



    // A timer to measure the complete cell loop time
    dealii::Timer total_timer;


    // Each thread get @p stride number of cells to work on
    uint stride = n_cells / numerics_prms.n_threads;

    // A vector to hold the CellWorker objects, each representing one
    // thread
    QVector< CellWorker<dim>* > worker_vector;


    // Here we start the CellWorker's.\n
    // total_timer measures the time of the cell loop, so the time needed for all
    // workers to finish
    total_timer.restart();
    for(uint cell_offset = 0; cell_offset < n_cells; cell_offset += stride){
        typename dealii::MGDoFHandler<dim>::active_cell_iterator thread_start_cell = cell_begin;
        // Since the @p active_cell_iterator's do not have an operator+, we need these
        // @p for() loops here to set the iterators to the right position for each thread
        for(uint i = 0; i < cell_offset; ++i)
            ++thread_start_cell;

        typename dealii::MGDoFHandler<dim>::active_cell_iterator thread_end_cell = thread_start_cell;
        for(uint i = 0; (i < stride) && (thread_end_cell != cell_end); ++i)
            ++thread_end_cell;

        worker_vector.append(new CellWorker<dim>(thread_start_cell, thread_end_cell,
                                                 numerics_prms,
                                                 fem_bem, constraints, matrices(LaplacePhi),
                                                 system_rhs, arch,
                                                 phys_data,
                                                 total_iteration_time, total_iteration_counter));

        // Use the Qt Signal-Slot mechanism to make sure the CellWorker objects
        // are deleted properly. This should work because, even though we have no
        // Qt event loop in the main thread, the signals and slots are both local to
        // the CellWorker object and should therefore be executed.
        QObject::connect( worker_vector.back(), SIGNAL(finished()),
                          worker_vector.back(), SLOT(deleteLater()) );
        // Start the threads
        worker_vector.back()->start();
    }

    // Wait until all threads are finished. The @p Q_FOREACH macro makes
    // this very elegant
    Q_FOREACH(CellWorker<dim> *worker, worker_vector){
        worker->wait();
    }
    timer.total_time = total_timer.wall_time();

    // The @p total_iteration_counter und @p total_iteration_time
    // are added up in the threads, so now they really hold the total values
    timer.iteration_counter = total_iteration_counter;
    timer.average_iteration_time = total_iteration_time;




    // Last but not least, add solid angles to the diagonal of the double layer matrix
    fem_bem.alpha_plus_DL_matrix = fem_bem.DL_matrix;
    fem_bem.alpha_plus_DL_matrix *= -1;
    for (uint i = 0; i < fem_bem.n_bc_dofs; i++){
        fem_bem.alpha_plus_DL_matrix(i,i) += fem_bem.solid_angles(i);
    }

    // compute boundary mass matrix for creating BEM contribution to FEM rhs

    {
        dealii::QGauss<dim-1> surface_quad(this->degree+1);

        fem_bem.sparsity_pattern.reinit(fem_bem.n_bc_dofs,
                                        fem_bem.n_bc_dofs,
                                        mg_dof_handler.max_couplings_between_dofs());


        dealii::DoFTools::make_boundary_sparsity_pattern(
                    static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                 fem_bem.interface,
                                                 fem_bem.dof_to_boundary_mapping,
                                                 fem_bem.sparsity_pattern
                                                 );

        fem_bem.sparsity_pattern.compress();

        fem_bem.mass_matrix.reinit(fem_bem.sparsity_pattern);
        dealii::Vector< double >  	dummy_rhs_vector(fem_bem.n_bc_dofs);
        dealii::ConstantFunction<dim>                dummy_fct(1., 1);
        fem_bem.interface[DRS::BCids::cavity()] = &dummy_fct;
        dealii::MatrixCreator::create_boundary_mass_matrix
                (mapping,
                 mg_dof_handler,
                 surface_quad,
                 fem_bem.mass_matrix,
                 fem_bem.interface,
                 dummy_rhs_vector,
                 fem_bem.dof_to_boundary_mapping
                 );
    }

    // And lastly, write out the measured times together with imporatant
    // problem size parameters
    timer.add_values_to_table(timing_table);


}



} // namespace Step16 END



// @sect4{Class: CellWorker}
// @sect5{Constructor: CellWorker}
//
// The constructor is responsible for setting up the various parameters
// needed to execute part of the cell loop.
template<int dim>
Step16::CellWorker<dim>::CellWorker(typename dealii::MGDoFHandler<dim>::active_cell_iterator _start_cell,
                                     typename dealii::MGDoFHandler<dim>::active_cell_iterator _end_cell,
                                     const step27::NumericsParams<dim> &_numerics_prms,
                                     typename Step16::FemBemData<dim> &_fem_bem,
                                     dealii::ConstraintMatrix &_constraints, dealii::SparseMatrix<double> &_system_matrix,
                                     dealii::Vector<double> &_system_rhs, const step28::Architecture _arch,
                                     const typename step27::PhysicalParams &_phys_data,
                                     double &_iterationTime, int &_iteration_counter)
    :
      start_cell(_start_cell),
      end_cell(_end_cell),
      numerics_prms(_numerics_prms),
      fem_bem(_fem_bem),
      constraints(_constraints),
      system_matrix(_system_matrix),
      system_rhs(_system_rhs),
      arch(_arch),
      phys_data(_phys_data),
      total_iteration_time(_iterationTime),
      total_iteration_counter(_iteration_counter)
{}


// @sect5{Function: run}
//
// overrides QThreads @p run() function. This is the function that gets executed
// in a thread when @p start() is called on a CellWorker object.
template<int dim>
void Step16::CellWorker<dim>::run()
{

    dealii::FE_Q<dim> fe(numerics_prms.fe_degree);
    dealii::MappingQ<dim> mapping(numerics_prms.fe_mapping_degree);

    const dealii::QGauss<dim>  quadrature_formula(fe.degree+1);


    typename DRSProblem<dim>::LaplaceForm laplace(mapping, fe,
                                                      quadrature_formula);


    // -------- BEM DATA BEGIN ---------------

    const dealii::QGauss<dim -1> regular_quad(numerics_prms.bem_quad_order_1);

    dealii::FEFaceValues<dim> bem_values (mapping,
                                  fe,
                                  regular_quad,
                                  dealii::update_values|dealii::update_normal_vectors|dealii::update_quadrature_points|dealii::update_JxW_values);

    const unsigned char bc_id       = *fem_bem.boundary_indicators.begin();
    const uint  dofs_per_be =  bem_values.dofs_per_cell; // THIS WOULD BE MORE INTUITIVE: fe.dofs_per_face;



    // Return the normal vectors at the quadrature points.
    // For a face, these are the outward normal vectors to the cell.
    // For a cell of codimension one, the orientation is
    // given by the numbering of vertices.
    const std::vector<dealii::Point<dim> > & normals = bem_values.get_normal_vectors ();

    const std::vector<dealii::Point<dim> > & q  = bem_values.get_quadrature_points ();

    const std::vector<dealii::Point<dim> > & x  = fem_bem.support_points[bc_id];
    const uint n_bem_quad_points= regular_quad.size(); //== q.size()
    const uint n_bem_points     = x.size();

    std::vector<uint> global_bc_dofs(dofs_per_be);

    const std::vector<double> & JxW     = bem_values.get_JxW_values();


    static const bool is_col_major = true;

    FullMatrixAccessor<double> W(n_bem_quad_points, dofs_per_be, is_col_major);


    dealii::Vector<double>       bem_local_rhs    (dofs_per_be);


    std::vector<double> norm_gradients(n_bem_quad_points);


    // -------- BEM DATA END ---------------


    step28::FullMatrixAccessor local_DL_matrix_k(n_bem_points, dofs_per_be, is_col_major);
    step28::FullMatrixAccessor local_SL_matrix_k(n_bem_points, dofs_per_be, is_col_major);

    step28::CUDADriver<dim> run(n_bem_points, dofs_per_be);


    std::vector<double> x_k;

    std::vector<double> q_k;
    std::vector<double> n_k;

    std::vector<double> rhs_for_dirichlet_test_k;

    // If we actually use CUDA...
    if(this->arch == step28::cuda)
    {
        // ...then we need to resize the containers to be big enough to actually hold
        // out data

        x_k.resize(dim * x.size());
        Step16::DRSProblem<dim>::point_to_double_vector(x_k, x);

        q_k.resize(dim * q.size());
        n_k.resize(dim * normals.size());

        rhs_for_dirichlet_test_k.resize(n_bem_points); //(==x.size())
    }


    dealii::Timer iteration_timer;
    // ------------------------- CELL LOOP BEGIN -----------------------
    double iteration_time = 0;
    uint iteration_counter = 0;



    typename dealii::MGDoFHandler<dim>::active_cell_iterator cell = start_cell;
    for (; cell!=end_cell; ++cell)
    {
        laplace.cell_matrix = 0;
        laplace.cell_rhs = 0;

        laplace.assemble (cell );


        // ------------- BEM stuff ------------
        if (cell->at_boundary() )
        {
            // Test: $1/2 u(x_i) = +V t_P - K \Phi_P$
            //
            // $t_P = d_{n_P} \Phi_P = - d_n \Phi_P $
            //
            // For BIEs we need the normal in the opposite direction
            // $V (-t_P)$.

            for (uint f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++){
                if(cell->face(f)->boundary_indicator() == bc_id ){

                    bem_local_rhs      = 0.;

                    bem_values.reinit (cell, f);


                    // Evaluate cell stuff only once
                    for (uint j = 0; j <dofs_per_be; ++j){
                        global_bc_dofs[j] = fem_bem.dof_to_boundary_mapping[laplace.local_dof_indices[j]];

                        for (uint a=0; a < n_bem_quad_points; ++a){
                            W(a,j) = bem_values.shape_value(j, a) * JxW[a];

                        }
                    }


                    // iterationTimer measures the time needed to construct the
                    // BEM matrices in each iteration.
                    iteration_timer.restart();
                    // switch-case for choosing the architecture at runtime
                    switch(this->arch)
                    {
                    // Run on GPU using CUDA:
                    case step28::cuda:
                    {
                        // Copy std::vector<Point<dim> > to std::vector<double>

                        Step16::DRSProblem<dim>::point_to_double_vector(q_k, q);
                        Step16::DRSProblem<dim>::point_to_double_vector(n_k, normals);

                        // Run the CUDA implementation
                        run.assemble_bem_matrices(local_DL_matrix_k, local_SL_matrix_k,
                                                  x_k, q_k, n_k, W);
                        // and write the values to global matrices
                        mutex.lock();
                        for (uint i=0; i < n_bem_points; ++i){

                            for (uint j = 0; j <dofs_per_be; ++j){
                                uint J = global_bc_dofs[j];
                                if (J != dealii::DoFHandler<dim>::invalid_dof_index){
                                    fem_bem.DL_matrix(i,J) +=  local_DL_matrix_k(i,j);
                                    fem_bem.SL_matrix(i,J) +=  local_SL_matrix_k(i,j);
                                }
                            }
                        }
                        mutex.unlock();
                        break;

                    }
                        // Unchanged CPU implementation
                    case step28::cpu:
                    {
                        // Reinitialize to zero the matrices and the vector
                        local_DL_matrix_k.reinit(n_bem_points, dofs_per_be);
                        local_SL_matrix_k.reinit(n_bem_points, dofs_per_be);

                        for (uint i=0; i < n_bem_points; ++i)
                        {

                            for (uint a=0; a < n_bem_quad_points; ++a)
                            {
                                dealii::Point<dim> R = x[i] - q[a];


                                double G_ia = LaplaceKernel<dim>::single_layer(R);

                                double H_ia = normals[a] * LaplaceKernel<dim>::double_layer(R);


                                for (uint j = 0; j < dofs_per_be; ++j)
                                {
                                    local_DL_matrix_k(i, j) +=  H_ia * W(a, j);
                                    local_SL_matrix_k(i, j) += -G_ia * W(a, j);
                                }
                            }
                        }

                        mutex.lock();
                        for (uint i=0; i < n_bem_points; ++i){

                            for (uint j = 0; j <dofs_per_be; ++j){
                                uint J = global_bc_dofs[j];
                                if (J != dealii::DoFHandler<dim>::invalid_dof_index){
                                    fem_bem.DL_matrix(i,J) +=  local_DL_matrix_k(i,j);
                                    fem_bem.SL_matrix(i,J) +=  local_SL_matrix_k(i,j);
                                }
                            }
                        }
                        mutex.unlock();
                        break;
                    }
                        // Run both and compare the results, to check CUDA
                    case step28::both:
                    {
                        // Not implemented here
#ifdef QT_NO_DEBUG
        AssertThrow(false,
        dealii::StandardExceptions::ExcMessage("Architecture \"both\" not implemented in multithreaded mode. "
        "Either reduce the number of threads to 1 or choose just one Architecture.") );
#else
        Assert(false,
        dealii::StandardExceptions::ExcMessage("Architecture \"both\" not implemented in multithreaded mode. "
        "Either reduce the number of threads to 1 or choose just one Architecture.") );
#endif
                        break;
                    }
                    }
                    // Add up the iteration time and increment the iteration counter
                    // so that we actually know how many iterations we measured

                    iteration_time += iteration_timer.wall_time();
                    ++iteration_counter;

                }
            }


        }


        mutex.lock();
        constraints.distribute_local_to_global (laplace.cell_matrix,
                                                laplace.cell_rhs,
                                                laplace.local_dof_indices,
                                                /*Laplace_Phi*/system_matrix, system_rhs);

        total_iteration_counter += iteration_counter;
        total_iteration_time += iteration_time;
        mutex.unlock();

    }

}



#endif // ASSEMBLE_SYSTEM_MULTITHREADED_HH
