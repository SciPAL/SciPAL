/*This file is part of SciPAL.

    SciPAL is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SciPAL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
*/


// STL header
#include <iostream>
#include <vector>

// QT
#include <QThread>
#include <QTime>

// deal.II-components
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

#include <QDir>


// Drivers for the GPU part.
// Include all other header files needed.
#include <step-1/SimParams.h>

#include <step-1/cuda_driver_step-1.h>
#include <step-1/cuda_driver_step-1.hh>

// Headers from SciPal.
#include <lac/FullMatrixAccessor.h>
#include <lac/MatrixCreator.h>


namespace step1 {

// @sect3{Class: CholeskyTest}
//
// This class is responsible for executing the Cholesky factorization in a separate thread.
// To this end we have to inherit from QThread and overwrite the run()-Function.
// @param Number : type of matrix entries.
template<typename Number>
class CholeskyTest : public QThread
{
    // We only want to print small matrices to the console.
    static const unsigned int max_display_size = 20;

public:

    CholeskyTest(int n_r,
                 dealii::ConvergenceTable& s_table,
                 const SimParams &_params);

protected:

    void setup_and_assemble_test_matrix();

    void factorize();

    void check_results();

    void run();



    int n_rows;

    dealii::ConvergenceTable & speedup_table;

    dealii::FullMatrix<Number> A, L, L_T;

private:

    const SimParams * params;
};


// @sect3{Class: Cholesky}
//
// For performance comparisons we also need a CPU-based Cholesky and LU factorization.
class Cholesky {

public:
    template<typename T> static void cpu(std::vector<std::vector<T> > & A);

    template<typename T> static void cpu_tiled(T* A, int tile_size);

    template<typename T> static void LLtMult(T * A, const T * L, int n_rows);
};


// @sect4{Constructor: CholeskyTest}
//
template<typename Number>
CholeskyTest<Number>::CholeskyTest(int n_r,
                                   dealii::ConvergenceTable& s_table,
                                   const SimParams &_params)
    :
      n_rows(n_r),
      speedup_table(s_table),
      params(&_params)
{}



// @sect4{Function: setup_and_assemble_test_matrix}
//
template<typename Number>
void CholeskyTest<Number>::setup_and_assemble_test_matrix()
{
    this->A.reinit(n_rows, n_rows);
    this->L.reinit(n_rows, n_rows);
    this->L_T.reinit(n_rows, n_rows);

    QTime t;

    this->speedup_table.add_value("n rows", n_rows);

    std::cout << "Initial Cholesky factor deal.II :" << std::endl;
    std::cout << "----------------------------------------------------"
              << std::endl;

    // For debugging purposes it is useful to
    // have the possibility of factorizing a symmetric and orthogonal matrix like
    // Hadamard matrices.
#ifdef USE_HADAMARD
    MatrixCreator::extended_hadamard(n_rows, A);
#else
    t.start();

    // Because of the design of the Tmmult-function of the class FullMatrix provided by deal.II
    // we have to compute the transpose of the reference
    // Cholesky factor $L^T_{rc} = (r+2)(c+2),
    // \quad r = 0,\ldots , n_{rows}-1\,, c = r,\ldots , n_{rows}-1$
    for (unsigned int r = 0; r < n_rows; ++r)
        for (unsigned int c = r; c <n_rows; ++c)
            L(r,c) = 1e-0*(r+2)*(c+2);


    if(false)
    qDebug("Time for CPU-based setup of Cholesky factor : %f s",
           t.elapsed()/1000.);

    if ( L.n_rows() < max_display_size)
        L.print(std::cout, 10, 5);
    if ( L.n_rows() < max_display_size)
        std::cout << std::endl;std::cout << std::endl;

    {
        // At this place it is instrcutive to measure the time
        // needed for copying matrices. To do this, we reinitialize the timer for measuring the copy time.
        t.restart();
        L_T.copy_from(L);
        if(false)
        qDebug("Time for CPU-based copying of Cholesky factor : %f s",
               t.elapsed()/1000.);

        // Then we construct the matrix which is to be factorized by computing $A = L\cdot L^T$.
        t.restart();
        L_T.Tmmult(A, L, false);
        if(false)
        qDebug("Time for CPU-based multiplication of Cholesky factor"
               " : %f s",
               t.elapsed()/1000.);
        // The Matrix @p A is backed up, as its recomputation is pretty expensive.
        L_T.copy_from(A);
    }

    if ( A.n_rows() < max_display_size)
    {
        std::cout << "Matrix to factorize :" << std::endl;
        std::cout << "----------------------------------------------------"
                  << std::endl;

        A.print(std::cout, 10, 5);

        std::cout << std::endl;std::cout << std::endl;
    }
#endif

}

// @sect4{Function: run}
//
// Initialization and factorization.
template<typename Number>
void CholeskyTest<Number>::run()
{
    // Prepare the test data.
    this->setup_and_assemble_test_matrix();

    // Timer
    QTime t;
    double cpu_time, gpu_time;

    // Execute CPU-based Cholesky factorization.
    FullMatrixAccessor<Number> A_h_cpu(A, true);
    {
        t.restart();

        // We need the temporary object @p A_h to get access to the value
        // array of @p A.

        Cholesky::cpu_tiled<Number>(A_h_cpu.val(), A_h_cpu.n_rows() );


        cpu_time =  t.elapsed()/1000.;
        if (false)
        qDebug("Time for CPU-based Cholesky factorization : %f s",
               cpu_time);

        // Save timing results in a deal.II convergence table.
        this->speedup_table.add_value("CPU factorization", cpu_time);
        this->speedup_table.set_precision("CPU factorization", 10);


        std::cout << "CPU-factorized Matrix (lower triangle contains "
                  << "transposed Cholesky factor) :" << std::endl;
        std::cout << "----------------------------------------------------"
                  << std::endl;

        if ( A.n_rows() < max_display_size)
            A_h_cpu.print(); //std::cout, 10, 5);

    }

    std::cout << std::endl;std::cout << std::endl;


    Assert( A.n_rows() == A.n_cols(),
            dealii::ExcMessage("Matrix not square! Cholesky impossible"));


    double kernel_time = 0;

    // Compute the Cholesky factorization on the GPU.
    CUDADriver<Number> run(*params);
    t.restart();
    {
        FullMatrixAccessor<Number> A_h(A, true);
        kernel_time = run.factorize(A_h);

        gpu_time =  t.elapsed()/1000.;

        // For the impatient user dump the results to screen.
        if (false)
        qDebug("Time for GPU-based Cholesky factorization %f"
               " including data transfer : %f s\n"
               "speed up factor factorization : %f netto : %f n_rows : %d\n",
               kernel_time,
               gpu_time,
               cpu_time/kernel_time,
               cpu_time/gpu_time,
               n_rows);


        std::cout << "GPU-factorized Matrix (lower triangle contains "
                  << "transposed Cholesky factor) :" << std::endl;
        std::cout << "----------------------------------------------------"
                  << std::endl;

        if ( A_h.n_rows() < max_display_size)
            A_h.print(); //std::cout, 10, 5);

        // Both factorizations should lead to the same result:
        // the original matrix still remains in the upper and the Cholesky
        // factor has appeared in the lower triangle. Thus, taking the norm
        // of their difference should yield a numerical zero.
        A_h_cpu -= A_h;

        std::cout << "difference of factorized matrices "
                        << " :" << std::endl;
              std::cout << "----------------------------------------------------"
                        << std::endl;

              if ( A_h_cpu.n_rows() < max_display_size)
                  A_h_cpu.print(); //std::cout, 10, 5);

              double F_norm =  A_h_cpu.frobenius_norm();

        std::cout << "||A_cpu - A_d ||_F = " << F_norm << "\n"
                  << "||A_cpu - A_d ||_F/n_el = " << F_norm/A_h_cpu.n_elements() << "\n"
                  << "||A_cpu - A_d ||_F/||A_d||_F = " << F_norm/A_h.frobenius_norm() << std::endl;

    }

    this->speedup_table.add_value("pure GPU fac", kernel_time);
    this->speedup_table.set_precision("pure GPU fac", 10);

    this->speedup_table.add_value("GPU fac incl data transfer", gpu_time);
    this->speedup_table.set_precision("GPU fac incl data transfer", 10);

    // Timing of individual components of the factorization in order
    // to compare manual version against CUBLAS-based variant
    FullMatrixAccessor<Number> A_h(A, true);
    FullMatrixAccessor<Number> A_original = A_h;

    {
        typename CUDADriver<Number>::TimerName2Value times;

        run.chol_fac(A_h, times);


        typename CUDADriver<Number>::TimerName2Value::const_iterator
                e=times.begin(),
                end_t=times.end();

        for( ; e != end_t ; ++e)
        {
            this->speedup_table.add_value(e->first, e->second);
            this->speedup_table.set_precision(e->first,10);
        }
    }
    return;
}


// @sect4{Function: cpu_tiled}
//
// @p T must be float or double.
//
// @param A : Matrix to factorize. factorization is in-place.
// Entries must be in the upper triangle including the diagonal.
//
// The Cholesky factor will be stored in the lower triangle and overwrites the diagonal.
// This function is identical to the Chol::__singe_thread() kernel.
template<typename T>
void Cholesky::cpu_tiled(T* A,
                         int tile_size)
{

    for (int r = 0; r < tile_size; ++r)
    {
        // Compute diagonal entry.
        T sum = 0.;
        int idx;
        int idx_c;

        for (int u = 0; u < r; ++u)
        {
            idx = r*tile_size + u;
            sum += A[idx] * A[idx];
        }
        idx = r*tile_size + r;
        A[idx] = sqrt(A[idx] - sum);

        for (int c = r+1; c < tile_size; ++c)
        {
            T tmp = 0.;

            for (int u = 0; u < r; ++u)
            {
                idx_c = c*tile_size + u;
                idx   = r*tile_size + u;
                tmp += A[idx_c]*A[idx];
            }

            idx_c = c*tile_size + r;
            idx   = r*tile_size + c;
            A[idx_c]  = A[idx] - tmp;
            A[idx_c] /= A[r*tile_size + r];
        }
    }
}


// @sect4{Function: LLtMult}
//
// Computes the original matrix.
// @param A : Pointer to value array of matrix
// @param L : Pointer to value array of Cholesky factor
// @param n_rows : Number of rows
template<typename T>
void Cholesky::LLtMult(T * A, const T * L, int n_rows)
{
    for (unsigned int r = 0; r < n_rows; ++r)
        for (unsigned int c = 0; c <=r; ++c)
        {
            unsigned int idx = c + (r*(r+1))/2;
            unsigned int k_max = std::min(r,c);

            A[idx] = 0.;

            for (unsigned int k = 0; k < k_max; ++k)
            {
                unsigned int idx_k   = k + (r*(r+1))/2;
                unsigned int idx_k_T = k + (c*(c+1))/2;

                A[idx] += L[idx_k]*L[idx_k_T];
            }
        }
}


// @sect3{Class: MyFancySimulation}
//
// The final class which drives the simulation.
// This class is primarily intended
// to manage the user's input.
template<typename Number>
class MyFancySimulation {

public:

    MyFancySimulation(SimParams &p);

    void run();

    static std::string precision_id();

private:
    const SimParams * params;

};



// @sect4{Constructor}
//
// The constructor of the simulation class sets the pointer to the runtime parameters
// and which GPU ("device") to use.
template <typename Number>
step1::MyFancySimulation<Number>::MyFancySimulation(SimParams &p)
    :
      params(&p)
{
    cudaSetDevice(params->device); 
}



// @sect4{Function: precision_id}
//
// Returns a string for identifying the precision.
template<>
std::string MyFancySimulation<float>::precision_id()
{
    return "float";
}

template<>
std::string MyFancySimulation<double>::precision_id()
{
    return "double";
}

}


// @sect4{Function: run}
//
// Compute the factorization for different matrix sizes.
template<typename Number>
void step1::MyFancySimulation<Number>::run()
{   

    // Before we start the test we setup the naem of the file
    // which at the end contains the runtimes.
    std::ostringstream filename;

    filename << "chol_fac_times_" << params->matrix_low << "_" << precision_id().c_str() << ".dat";


    // The results of the factorization are stored in a convergence table.
    dealii::ConvergenceTable factorization_times;

    // Loop over all matrix sizes in the specified range.
    for (int n = params->matrix_low; n < params->matrix_high; n+=params->step_size)
    {
        CholeskyTest<Number> driver(n, factorization_times, *params);

        driver.start();
        // For debugging purposes it is sometimes useful to disable
        // the inheritance from QThread and to call the run()-function directly.
        /* driver.run();*/
        driver.wait();

        // To avoid data loss we save the results after each factorizatrion.
        std::ofstream out(filename.str().c_str());
        factorization_times.write_text(out);
    }

    std::cout << "Done." << std::endl;
}




// @sect4{Funktion: main}
//
// Instantiate and execute the simulation.
int main(int argc, char *argv[])
{
    using namespace step1;

    SimParams params;

    // First we declare the parameters to expect ...
    dealii::ParameterHandler prm_handler;

    // Get the current working directory
    QDir cwd = QDir::current();

    // and backup the location where the program has been started.
    // Here, we assume the we use QTCreators shadow-biuld mechanism
    // whoch puts the build directory at the same level as the directory <i>step-1</i>
    // containing the source code.
    const QDir launch_dir = cwd;
    cwd.setPath("../step-1");

    // By default, the parameter file has the same name as the binary
    // and is supposed to be in a subdirectory prm of the directory,
    // where the program has been started.
    std::string prm_filename;
    if (argc == 1)
    {
        std::string tmp = argv[0];
        int found=tmp.find_last_of('/');
        prm_filename = tmp.substr(found+1);
        prm_filename += "-Decomp.prm";

        cwd.setPath("./prm");
    }
    else
    {
        QFileInfo tmp(argv[1]);

        // Subdivide the given filename into its path and filename
        // so that the corresponding subdirectories can be created.
        QString prm_path = tmp.absolutePath();
        cwd.setPath(prm_path);
        cwd.makeAbsolute();
        prm_filename = tmp.fileName().toStdString();

        std::cout << "chosen prm file : " << tmp.absoluteFilePath().toStdString().c_str() << std::endl;
    }

    // Before the parameter file can be read, we have to make sure that
    // its directory exists
    if (!cwd.exists() )
        launch_dir.mkpath( cwd.absolutePath() );

    QDir::setCurrent(cwd.absolutePath());
    SimParams::declare(prm_handler);
    prm_handler.read_input (prm_filename);

    QDir::setCurrent(launch_dir.absolutePath());

    params.get(prm_handler);

    // Create the toplevel run directory.
    cwd.setPath(params.run_dir.absolutePath());
    // The following lets a directory make its own path.
    if (!cwd.exists())
        cwd.mkpath( "." );

    // Now, change to the run dir
    QDir::setCurrent(cwd.absolutePath());

    cwd.setPath("./log");
    cwd.makeAbsolute();
    if (!cwd.exists())
        cwd.mkpath(".");

    // Create the log directory and write what has been actually read
    // into log file. Basically, this is just another parameter file
    // and can thus be used again as input to another run.
    QDir::setCurrent(cwd.absolutePath());

    prm_filename += ".log";
    std::ofstream log_out_text(prm_filename.c_str());
    prm_handler.print_parameters (log_out_text,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    QDir::setCurrent(params.run_dir.absolutePath());


    // Now, run the comparison of GPU vs. CPU for the selected precision.
    if (!params.use_double) {

        MyFancySimulation<float> machma_float(params);

        machma_float.run();
    }
    else {

        MyFancySimulation<double> machma_double(params);

        machma_double.run();
    }
}
