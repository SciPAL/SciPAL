#ifndef STEP_4_HH
#define STEP_4_HH


// STL header
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>

// deal.II
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

// QT
#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>

// SciPal
#include <base/CUDATimer.h>
#include <lac/blas++.h>
#include <lac/MatrixCreator.h>

// Driver for GPU part of the program
#include <cuda_driver_step-4.h>
#include <cuda_driver_step-4.hh>


// We start the commenting walk through the code of this program
// by some structures needed
// for retrieving information about the GPU computing capabilities of the system
// and defining a simple file-based user interface (UI).
//
// TODO: replace by GPUInfo.hh
// @sect3{Struct: GlobalData}
//
// This structure collects general information about the available graphics cards
// and upon request displays the most important ones. For more details of the
// @p cudaDeviceProp structure have a look a the CUDA reference manual.
struct GlobalData {

    int n_CUDA_devices;

    int current_device_id;

    cudaDeviceProp prop;


    void cublanc_gpu_info()
    {
        const int KB = 1024;
        const int MB = KB*KB;

        // Ask for the properties of a GPU with id @p current_device_id.
        // This attribute must be set before calling this function.
        cudaGetDeviceProperties(&prop, this->current_device_id);

        printf("Currently used GPU: %s \n",prop.name);
        printf("Compute Capability: %d.%d \n",prop.major,prop.minor);
        printf("ClockRate: %uMHz \n",prop.clockRate/1000);
        printf("Warpsize: %d \n",prop.warpSize);
        printf("Number of Multiprocessors: %d \n",prop.multiProcessorCount);

        printf("Shared Memory: %lKB\n",prop.sharedMemPerBlock/KB);
        printf("Constant Memory: %lKB \n",prop.totalConstMem/KB);
        printf("Global Memory: %l"
               "MB \n",prop.totalGlobalMem/MB);
        printf("the device %s can concurrently copy memory "
               "between host and device while executing a kernel\n",
               (prop.deviceOverlap? "can": "cannot"));
    }

};

// For simplicity the GPU data is instantiated as global object.
GlobalData global_data;

// @sect3{Struct: PrecisionName}
//
// This is an auxiliary structure for generating precision-dependent
// filenames for the output of results.
template<typename T> struct PrecisionName { static std::string name(); };

template<>
std::string PrecisionName<float>::name()  { return "float"; }

template<>
std::string PrecisionName<double>::name() { return "double"; }

template<typename T>
std::string PrecisionName<T>::name()      { return "unknown_T"; }


// The whole project has its own namespace so that it becomes possible to
// reuse each step in another one even if classes have the same names.
// Reuse of individual steps into more complex applications is, for instance.
// demonstrated in step-6 and step-7.
namespace step4 {

// @sect3{Struct: TestUIParamsBase}
//
// This structure collects general information about the range of matrix sizes
// for which a factorization method should be tested. This is not only
// needed for QR but also for SVD or other matrix factorizations.
// To reuse this class simply inherit from it.
struct TestUIParamsBase {

    TestUIParamsBase() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    // The following attributes indicate the bounds of the number of colums
    // of the test matrix.
    int  min_n_cols;

    int  max_n_cols;

    // To sample the problem sizes we will let the number of rows grow by this
    // factor when looping over the problem sizes in the factorization test.
    double n_rows_growth_rate;

    // The following attribute indicates how many times a test should be
    // repeated in order to compute a reliable average runtime.
    int n_repetitions;

    template <typename T, typename blas>
    void save_results(const std::string test_name,
                      const dealii::ConvergenceTable & results_table);

private:
    // To keep the compiler from generating unwanted copy constructors
    // or assignment operators we provide dummy implementations and
    // declare them as private. This will result in a compile-time error
    // if one of them accidentally is invoked.
    TestUIParamsBase(const TestUIParamsBase & ) {}

    TestUIParamsBase & operator = (const TestUIParamsBase & /*other*/)
    {
        return *this;
    }
};

} // namespace step4 END

// @sect4{Function: declare}
//
// Parameters will be read by an object of type dealii::ParameterHandler.
// This class requires to declare which parameters should be read which is done
// by this function. For details have a look at the manual of deal.II.
// This function must be static so that parameters can be declared before the object
// is instantiated which will hold them.
// @param prm : Object of type dealii::ParameterHandler
void step4::TestUIParamsBase::declare(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Dimensions of test problems.");

    prm.declare_entry ("log2(min n cols)",
                       "2",
                       dealii::Patterns::Double(2,12),
                       "Binary logarithm of minimal number of columns of "
                       "upper triangular matrix R. "
                       "Allowed range : [2,12]");

    prm.declare_entry ("log2(max n cols)",
                       "3",
                       dealii::Patterns::Double(2,13),
                       "Binary logarithm of minimal number of columns "
                       "of upper triangular matrix R. "
                       "Allowed range : [3,13]");

    prm.declare_entry ("n rows growth rate",
                       "1.2",
                       dealii::Patterns::Double(1.1,10),
                       "In each test instance the number of rows of R is "
                       "increased by this factor. "
                       "Allowed range : [1.1,10]");

    prm.declare_entry("n repetitions",
                      "1",
                      dealii::Patterns::Integer(1),
                      "Repeat the test for a given matrix size this many times.");

    prm.leave_subsection();
}


// @sect4{Function: get}
//
// @param prm : Object of type dealii::ParameterHandler
void step4::TestUIParamsBase::get(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Dimensions of test problems.");

    min_n_cols = pow(2., prm.get_double ("log2(min n cols)") );

    max_n_cols = pow(2., prm.get_double ("log2(max n cols)") );

    n_rows_growth_rate = prm.get_double("n rows growth rate");

    n_repetitions = prm.get_integer("n repetitions");

    prm.leave_subsection();
}

// @sect4{Function: save_results}
//
// Factorization tests basically create data indicating the relationship
// between runtime, possibly differentiating w.r.t. different parts of an algorithm,
// and problem size which can be collected in a deal.II ConvergenceTable.
// At the end the results will have to be saved to some file. To distinguish
// different runs the filename should contain some tags indicating what has been tested.
// For instance, precision or type of GPU used. The analysis is simplified if the names
// are created according to a fixed format as this makes the composition of
// the name predictable. In the body of the fucntion you can also see, that
// QT's string classes are as easy to use as the ones provided by python.
//
// @param test_name : name of the test, e.g. "QR_test".
// @param results_table : table containing the runtime data.
template <typename T, typename blas>
void
step4::TestUIParamsBase::save_results(const std::string test_name,
                               const dealii::ConvergenceTable & results_table)
{
    std::string filename(test_name);

    filename += ("_results_");

    filename += blas::name();

    filename += "_";

    if (blas::name() == "cublas")
    {
        printf("Currently used GPU: %s \n", global_data.prop.name);
        filename += QString(global_data.prop.name).replace(" ", "_").toStdString().c_str();

        filename += "_";
    }

    filename += PrecisionName<T>::name();

    filename += ".out";

    std::ofstream out(filename.c_str());

        // At the beginning we dump the parameters indicating the range of problem sizes.
        // This makes tests repeatable.
    out << "# max n_cols  : " << this->max_n_cols << std::endl
        << "# min n_cols  : " << this->min_n_cols << std::endl
        << "# growth rate : " << this->n_rows_growth_rate << std::endl
        << std::endl;

    results_table.write_text(out);

}


namespace step4 {

// @sect3{Struct: QRTestUIParams}
//
// The parameters for testing the QR factorization are completed by
// the attributes of this structure.
// A bonus of a file-based user interface is, that it can be easily integrated
// into a continuous integration testsuite which runs overnight and provides
// a per-day check whether changes of external libaries have introduced
// any errors in a program.
struct QRTestUIParams : public TestUIParamsBase {

    QRTestUIParams() : TestUIParamsBase() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);


    // The following flags indicate which matrices should be dumped
    // to screen for visual inspection.
    bool  print_Q_orig;

    bool  print_R_orig;

    bool  print_A_orig;

    bool  print_Q ;

    bool  print_R;

    bool  print_QR;

    bool  print_QtQ;

    bool  print_Q_col_errors;


private:
    QRTestUIParams (const QRTestUIParams & /*other*/) : TestUIParamsBase() {}

    QRTestUIParams & operator = (const QRTestUIParams & /*other*/) { return *this; }

};

} // namespace step4 END

// @sect4{Function: declare}
//
// @param prm : Object of type dealii::ParameterHandler
void
step4::QRTestUIParams::declare(dealii::ParameterHandler & prm)
{
    // Start with declaring the basic parameters, that is the range of
    // problem sizes.
    TestUIParamsBase::declare(prm);

    // In a separate subsection we declare the run-time parameters
    // specific to this test.
    prm.enter_subsection("Visualization flags.");

    prm.declare_entry ("print original Q",
                       "false",
                       dealii::Patterns::Bool(),
                       "");

    prm.declare_entry ("print original R",
                       "false",
                       dealii::Patterns::Bool(),
                       "");

    prm.declare_entry ("print original A",
                       "false",
                       dealii::Patterns::Bool(),
                       "The original A is given multiplying the original "
                       "Q with the original R.");


    prm.declare_entry ("print Q from factorization",
                       "false",
                       dealii::Patterns::Bool(),
                       "");

    prm.declare_entry ("print R from factorization",
                       "false",
                       dealii::Patterns::Bool(),
                       "");

    prm.declare_entry ("print QR",
                       "false",
                       dealii::Patterns::Bool(),
                       "Multiplying Q and R as obtained from the factorization "
                       "should give something that closely resembles the original A");

    prm.declare_entry ("print QtQ from factorization",
                       "false",
                       dealii::Patterns::Bool(),
                       "Multiplying the Q obtained from the factorization with "
                       "its own transpose should give a unit matrix.");

    prm.declare_entry ("print Q column errors",
                       "false",
                       dealii::Patterns::Bool(),
                       "Computing ||Q_orig - Q||_2 must be done column-wise "
                       "as the individual columns of Q are reproduced only "
                       "up to a sign (the algorithm cannot distinguish "
                       "the direction of a column vector).");

    prm.leave_subsection();

}


// @sect4{Function: get}
//
// @param prm : Object of type dealii::ParameterHandler
void
step4::QRTestUIParams::get(dealii::ParameterHandler & prm)
{

    this->TestUIParamsBase::get(prm);


    prm.enter_subsection("Visualization flags.");

    print_Q_orig       = prm.get_bool ("print original Q");

    print_R_orig       = prm.get_bool ("print original R");

    print_A_orig       = prm.get_bool ("print original A");


    print_Q            = prm.get_bool ("print Q from factorization");

    print_R            = prm.get_bool ("print R from factorization");

    print_QR           = prm.get_bool ("print QR");

    print_QtQ          = prm.get_bool ("print QtQ from factorization");

    print_Q_col_errors = prm.get_bool ("print Q column errors");

    prm.leave_subsection();
}


namespace step4 {

// @sect3{Class: QRTest}
//
// This class manages the performance test. It is responsible for generating
    // test matrices $A\in \mathbb{R}^{m \times n}$ from a given orthogonal
    // matrix $Q \in \mathbb{R}^{m \times m}$ and an
    // upper triangular matrix $R \in \mathbb{R}^{m \times n}$.
    // The test is passed if $A$ can be recovered
    // from the results $Q_{num}$ and $R_{num}$ of the factorization up
    // to numerical accuracy, i.e.
    // \f{eqnarray*}
    // \| A - Q_{num} R_{num} \| \le 10^{-16}C(mn)
    // \f}
    // where $C(mn)$ is a problem-size dependent constant.
    // Division by the problem size then gives an estimate for the
    // element-wise quality of the factorization.
    // The all-govering parameter for the test the number of rows $m$,
    // as the set of numbers of columns ${n}$ is chosen in dependence on that.
    // Then, for each pair $(m,n)$ a test matrix is generated and the
    // factorization error and runtime is measured.
    // The second template argument makes this class independent of the
    // BLAS implementation.
template <typename T, typename blas>
class QRTest {

public:

    typedef typename blas_pp<T, blas>::FullMatrixAccessor FullMatrixAccessor;

    typedef typename blas_pp<T, blas>::Matrix       Matrix;

//    typedef          SciPAL::transpose<Matrix>              tr;

    typedef typename blas_pp<T, blas>::SubMatrix    SubMatrix;

    typedef typename blas_pp<T, blas>::MatrixSubCol MatrixSubCol;

    typedef typename blas_pp<T, blas>::Vector       Vector;



    QRTest(dealii::ParameterHandler & prm);

    void run();

protected:

    void check_results(const dealii::FullMatrix<T> & A,
                       const step4::CudaQRDecomposition<T, blas>& QRf,
                       T elapsed_time);

private:
    void setup_and_assemble(unsigned int nr, unsigned int nc);

    void factorize();

    void save_results();

    dealii::FullMatrix<T> Q, R, H;

    QRTestUIParams params;

    dealii::ConvergenceTable results_table;
};

} // namespace step4 END

// @sect4{Constructor: QRTest}
//
// @param prm : ParameterHandler containing the parameter values read
// from the prm file.
template <typename T, typename blas>
step4::QRTest<T, blas>::QRTest(dealii::ParameterHandler & prm)
{
    params.get(prm);
}


// @sect4{Function: run}
//
// Executes all the individual tests.
// It is made sure, that the number of rows grows at least by 1.
// Otherwise we would end up in an infinite loop.
// The number of columns simply gets doubled until the allowed maximum is
// reached. The number of rows starts at the current number of columns and
// grows by the prescribed growth factor or at least by 1 until the
// maximum number of columns is reached.
template <typename T, typename blas>
void
step4::QRTest<T, blas>::run()
{
    for (int nc = this->params.min_n_cols;
         nc < this->params.max_n_cols; nc*=2)
        for (int nr = nc; nr <= this->params.max_n_cols;
            nr=std::max(nr+1, int(nr*this->params.n_rows_growth_rate) ) )
        {
            setup_and_assemble(nr, nc);
            factorize();
        }

    setup_and_assemble(this->params.max_n_cols,
                       this->params.max_n_cols);

    factorize();


    save_results();
}



// @sect4{Function: setup_and_assemble}
//
// The matrix @p Q is initialized by a Hadamard-Matrix $H_{nr}$
// and @p R with an upper triangular matrix of width @p nc.
// @param nr : dimension of Q, i.e. the number of constraints to fulfill in a least-squares problem
// @param nc : columns of @p R, i.e. the number of unknowns.
template <typename T, typename blas>
void
step4::QRTest<T, blas>::setup_and_assemble(unsigned int nr, unsigned int nc)
{

    Assert(nr >= nc, dealii::ExcMessage("n_cols must not exceed n_rows."));

    int log2_nc = std::log(nc)/std::log(2);

    // As orthogonal matrix we take a normalized
    // <a href="http://en.wikipedia.org/wiki/Hadamard_matrix">Hadamard
    // matrix</a>.
    // Redundant columns are filled with unit matrix.
    MatrixCreator::hadamard(log2_nc, this->H);

    this->H /= std::sqrt(this->H.n_rows() );

    this->Q.reinit(nr, nr);
    // The upper left diagonal block is filled with the Hadamard matrix ...
    for (unsigned int r = 0; r < nc; ++r)
        for (unsigned int c = 0; c < nc; ++c)
            Q(r,c) = H(r,c);
    // and the lower right diagonal block is filled by a unit matrix.
    // In case of a square matrix this block does not exist.
    for (unsigned int r = nc; r < nr; ++r)
        Q(r,r) = 1;

    // The upper triangular matrix which is to be reconstructed
    // is filled row-wise with $1+sin(n)$, $n\in \mathbb{N}$
    // where n is continuously increased.
    // Thus, the matrix entries are in the interval $(0,2)$.
    this->R.reinit(nr, nc);
    int tmp = 1;
    for (unsigned int r = 0; r < nr; ++r)
        for (unsigned int c = r; c < nc; ++c)
            R(r,c) = 1+std::sin(tmp++);


    // Upon request the original matrices are dumped to screen prior to the
    // factorization.
    if (this->params.print_Q_orig) {
        std::cout
            << "Initial Q (Hadamard-Matrix, siehe wikipedia): \n";
        std::cout
            << "----------------------------------------------" << std::endl;
        Q.print(std::cout);
    }

    if (this->params.print_R_orig) {
        std::cout
            << "\nInitial R : \n";
        std::cout
            << "----------------------------------------------" << std::endl;
        R.print(std::cout, 12, 5);
    }

    if (this->params.print_Q_orig || this->params.print_R_orig)
        std::cout << std::endl;
}


// @sect4{Function: factorize}
//
// Create the test matrix @p A, apply the Householder-based
// QR-factorization to it and check the results.
template <typename T, typename blas>
void
step4::QRTest<T, blas>::factorize()
{

    dealii::FullMatrix<T> A(R.n_rows(), R.n_cols() );

    dealii::FullMatrix<T> Q_orig(Q);

    Q_orig.mmult(A, R);

    step4::CudaQRDecomposition<T, blas> QR;

    T elapsed_time = QR.householder(A);

    check_results(A, QR, elapsed_time);
}



// @sect4{Function: check_results}
//
//
template <typename T, typename blas>
void
step4::QRTest<T, blas>::check_results(const dealii::FullMatrix<T> & A,
                                    const step4::CudaQRDecomposition<T, blas>& QRf,
                                    T elapsed_time)
{
    int n_rows = QRf.R().n_rows();
    int n_cols = QRf.R().n_cols();

    std::cout
            << "n_rows : " << n_rows
            << ", n_cols : " << n_cols
            << std::endl;

    // Check whether $A = QR$, $Q^TQ = I$ is fulfilled.
    // As BLAS usually addresses matrices in column-major order
    // we transpose the input matrix while copying it because in deal.II
    // matrices are stored row-major internally.
    // The FullmatrixAccessor class is needed, as deal.II's FullMatrix class
    // does not grant access to the internal memory so that we need some proxy
    // in order to be able to copy it onto the GPU when using CUBLAS.
    FullMatrixAccessor A_t_ugly(A,
                                true /*column-major copy*/);

    Matrix A_backup;


    A_backup = A_t_ugly;

    // We use the verification that the original matrix can
    // be recovered from ${\bf Q}\cdot {\bf R}$ as an opportunity to test
    // the multiplication of a matrix view with a matrix.
    const SubMatrix Q_num(const_cast<Matrix&>(QRf.Q()), 0, 0);
    const SubMatrix R_num(const_cast<Matrix&>(QRf.R()), 0, 0);

    Matrix QR(QRf.Q().n_rows(), QRf.R().n_cols());

    SubMatrix QR_num(QR, 0, 0);

    QR_num = Q_num * R_num;

    // At the end, @p QtQ should be a diagonal matrix provided @p QRf.Q()
    // is orthogonal.

    const Matrix & Q_tmp = QRf.Q();


    Matrix QtQ;
    QtQ = transpose(QRf.Q()) * Q_tmp;

    // Compute deviations due to factorization
    Matrix  A_m_QR;
    A_m_QR = A_backup;
    A_m_QR -= QR;

    T l2_A_QR_error = A_m_QR.l2_norm();

    std::cout << "||A -QR||_2 : " << l2_A_QR_error << std::endl;

    // Up to numerical accuracy $Q^TQ$ should recover an identity matrix.
    Matrix I;
    dealii::IdentityMatrix I_h(n_rows);
    I = I_h;
    QtQ -= I;

    T l2_QtQ_error = QtQ.l2_norm();

    std::cout << "||Q^T Q - I||_2 : " << l2_QtQ_error << std::endl;

    // To determine $\|Q-Q_{orig}\|_2$ we have to take care of the signs of
    // the entries.
    // We first compute $\|Q-Q_{orig}\|_2$ column-wise
    // and then we compute the l2-Norm of the vector of the per-column errors.
    FullMatrixAccessor Q_t_ugly(Q,
                                   true /*column-major copy*/);
    Matrix Q_orig; Q_orig = Q_t_ugly;

    T l2_Q_error = 0;
    for (int c = 0; c < n_cols; ++c)
    {
         MatrixSubCol Q_col(const_cast<Matrix&>(QRf.Q()),0,c);
         MatrixSubCol Q_orig_col(Q_orig,0,c);

         Q_orig_col -= Q_col;

         T c_m_err = Q_orig_col.l2_norm();

         // Undo and change sign.
         Q_orig_col += Q_col;
         Q_orig_col += Q_col;
         T c_p_err = Q_orig_col.l2_norm();

         l2_Q_error += std::min(c_m_err, c_p_err);
         if (this->params.print_Q_col_errors)
         std::cout
                 << "Error in col " << c
                 << " - : " << c_m_err << ", + : " << c_p_err
                 << std::endl;
    }

    l2_Q_error = std::sqrt(l2_Q_error);

    std::cout << "||Q_orig - Q||_2 : " << l2_Q_error << std::endl;

    // Store errors for later plotting.
    this->results_table.add_value("n_rows", n_rows);
    this->results_table.add_value("n_cols", n_cols);

    this->results_table.add_value("||A -QR||_2", double(l2_A_QR_error) );
    this->results_table.set_precision	("||A -QR||_2", 16);

    this->results_table.add_value("||Q^T Q - I||_2", double(l2_QtQ_error) );
    this->results_table.set_precision	("||Q^T Q - I||_2", 16);

    this->results_table.add_value("||Q - Q_orig||_2", double(l2_Q_error) );
    this->results_table.set_precision	("||Q - Q_orig||_2", 16);

    this->results_table.add_value("elapsed time", double(elapsed_time) );
    this->results_table.set_precision	("elapsed time", 16);

    // Print the matrices to compare it in a visual way
    // Ausgabe der Matrizen fuer den visuellen Vergleich.
    if (this->params.print_Q_orig || this->params.print_R_orig)
        std::cout << "\nAt this point R and Q are done :" << std::endl;

    if (this->params.print_Q_orig) {
        Q_orig = Q_t_ugly;
        std::cout << "\nQ_orig :" << std::endl; Q_orig.print();
    }

    if (this->params.print_Q) {
        std::cout << "\nQ_num :" << std::endl; QRf.Q().print();
    }

    if (this->params.print_R) {
        SubMatrix R_upper_square(const_cast<Matrix&>(QRf.R()),
                                 0, n_cols, 0, n_cols);
        std::cout << "\nR :" << std::endl; R_upper_square.print();
    }


    if (this->params.print_QtQ) {
        std::cout << "\nQ^TQ :" << std::endl; QtQ.print();
    }


    if (this->params.print_A_orig) {
        SubMatrix A_upper_square(A_backup, 0, n_rows, 0, n_cols);
        std::cout << "\nA_orig :" << std::endl; A_upper_square.print();
    }

    if (this->params.print_QR) {
        SubMatrix QR_upper_square(QR, 0, n_rows, 0, n_cols);
        std::cout << "\nQ_num * R_num :" << std::endl; QR_upper_square.print();
    }
}

// @sect4{Function: save_results}
//
// Save the content of @p results_table in a file in  ascii format.
template <typename T, typename blas>
void
step4::QRTest<T, blas>::save_results()
{
    this->params.template save_results<T, blas>(std::string("QRTest"),
                                                results_table);
}

namespace step4 {

// @sect3{Struct: SimParams}
//
    // This structure contains all parameters necessary for controling
    // the global test properties, i.e. precision and what BLAS to use.
struct SimUIParams {

    SimUIParams() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);


    bool run_cublas_float, run_cublas_double;
    bool run_cpu_blas_float, run_cpu_blas_double;

};

} // namespace step4 END

// @sect4{Function: declare}
//
// This function informs the ParameterHandler what parameters it
// should read from the parameter file.
void
step4::SimUIParams::declare(dealii::ParameterHandler & prm)
{

    prm.enter_subsection("Global parameters");


    prm.declare_entry("Run cublas float", "false",
                      dealii::Patterns::Bool(),
                      "Single precision is offered by any CUDA-capable GPU.");

    prm.declare_entry("Run cublas double", "false",
                      dealii::Patterns::Bool(),
                      "Only available on graphics cards of compute capability "
                      "1.3 and higher");

    prm.declare_entry("Run CPU-BLAS float", "false",
                      dealii::Patterns::Bool(),
                      "Depending on the platform this is either the reference "
                      "BLAS from www.netlib.org or a tuned vendor-specific "
                      "version. For instance Intel's MKL, the Accelerate "
                      "framework by Apple or AMDs ACML.");

    prm.declare_entry("Run CPU-BLAS double", "false",
                      dealii::Patterns::Bool(),
                      "Cf. test for floats");


    prm.leave_subsection();

}


// @sect4{Function: get}
//
// This function must be called in order to transfer the parameters read
// by the ParameterHandler into the attrbutes of this object.
void
step4::SimUIParams::get(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Global parameters");


    run_cublas_float    = prm.get_bool("Run cublas float");

    run_cublas_double   = prm.get_bool("Run cublas double");

    run_cpu_blas_float  = prm.get_bool("Run CPU-BLAS float");

    run_cpu_blas_double = prm.get_bool("Run CPU-BLAS double");


    prm.leave_subsection();
}


namespace step4 {

        // @sect4{Function: run_qr_tests}
        //
        // Process user input and run the tests.
void run_qr_tests(int /*argc*/, char *argv[])
{
    // Declare all parameter and read it from a file.
    // The name of the file is the same like this program
    // extended with the ending '.prm'.
    // The read parameter will be written into a log-file.

    // First we declare the parameters to expect ...
    dealii::ParameterHandler prm_handler;

    SimUIParams::declare(prm_handler);

    QRTestUIParams::declare(prm_handler);

    // ... then read them from a file residing in the same directory
    // as the binary and with a name starting with the program's name.
    std::string prm_filename(argv[0]);
    prm_filename += "-QR.prm";
    prm_handler.read_input (prm_filename);

    prm_filename += ".log";
    std::ofstream log_out_text(prm_filename.c_str());
    prm_handler.print_parameters (log_out_text,
                                  dealii::ParameterHandler::Text);

    SimUIParams params;

    params.get(prm_handler);


    for(int DevNo = 0; DevNo < global_data.n_CUDA_devices; DevNo++)
    {
        cudaSetDevice(DevNo);

        global_data.current_device_id = DevNo;

        global_data.cublanc_gpu_info();


        if (params.run_cublas_float) {
            std::cout << "Householder-QR<float> using cublas :\n"
                    << "------------------------------------------------"
                    << std::endl;

            // Test cublas
            QRTest<float, cublas> cublas_qr_test_float(prm_handler);

            cublas_qr_test_float.run();
            std::cout << "\nHouseholder-QR<float> using cublas DONE \n"
                    << "------------------------------------------------\n"
                    << std::endl;
        }


    if (params.run_cublas_double) {
        std::cout
            << "\nHouseholder-QR<double> using cublas :\n"
            << "------------------------------------------------"
            << std::endl;
        QRTest<double, cublas> cublas_qr_test_double(prm_handler);

        cublas_qr_test_double.run();

        std::cout
            << "\nHouseholder-QR<double> using cublas DONE \n"
            << "------------------------------------------------\n"
            << std::endl;
    }

}

// An almost verbatim copy of the preceeding lines gives the tests
    // for the CPU version.
    if (params.run_cpu_blas_float) {
        std::cout
                << "\nHouseholder-QR<float> using ATLAS/CBLAS :\n"
                << "------------------------------------------------"
                << std::endl;

        // Test  blas
        QRTest<float, blas> blas_qr_test_float(prm_handler);

        blas_qr_test_float.run();
        std::cout
                << "\nHouseholder-QR<float> using ATLAS/CBLAS DONE \n"
                << "------------------------------------------------\n"
                << std::endl;
    }



    if (params.run_cpu_blas_double) {
        std::cout
                << "\nHouseholder-QR<double> using ATLAS/CBLAS :\n"
                << "------------------------------------------------"
                << std::endl;
        QRTest<double, blas> blas_qr_test_double(prm_handler);

        blas_qr_test_double.run();
        std::cout
                << "\nHouseholder-QR<double> using ATLAS/CBLAS DONE \n"
                << "------------------------------------------------\n"
                << std::endl;
    }
}

} // namespace step4 END
#endif // STEP_4_HH
