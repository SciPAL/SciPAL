//@sect3{File: cusolverDn_wrapper.h}
#ifndef CUSOLVERDN_WRAPPER_STEP43_HH
#define CUSOLVERDN_WRAPPER_STEP43_HH

//STL includes
#include <string>
#include <iostream>

//CUDA includes
#include <cusolver_common.h>
#include <cusolverDn.h>

//BASE includes
#include <base/ParallelArch.h>
#include <base/ArchTraits.h>
#include <base/VTraits.h>

//LAC includes
#include <lac/cublas_wrapper.hh>

namespace SciPAL
{
static cusolverDnHandle_t cusolver_handle;

// @sect4{struct cusolverDn: cusolverDn S,D,C,Z wrapper functions}
//
// We remove the type dependence of the function of cusolverDn by using C++'s polymorphism and templates.
// For the numerical purpose of the functions have a look at the cusolverDn documentation.
//

struct cusolverDn
{

    // @sect5{Funktion: name (cusolverDn)}
    //
    // Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "cusolverDn"; }

    // Compile-time variable for identifying the side of PCIe bus this solver works on.
    static const ParallelArch arch = gpu_cuda;

    typedef typename archTraits<arch>::BlasType BW;

    // @sect5{Funktion: check_status}
    //
    // Evaluate the @p cusolverStatus argument and throws
    // an exception in case an error occurs.
    //
    // \param cusolverStatus : status flag returned by a cuSolver function
    //
#ifdef DEBUG
    static void  check_status(const cusolverStatus_t & status)
    {
        std::string cusolverDn_errors(" ");

        if (status != CUSOLVER_STATUS_SUCCESS)
        {
            if(status == CUSOLVER_STATUS_ALLOC_FAILED)
                cusolverDn_errors += "CUSOLVER_STATUS_ALLOC_FAILED";
            if(status == CUSOLVER_STATUS_ARCH_MISMATCH)
                cusolverDn_errors += "CUSOLVER_STATUS_ARCH_MISMATCH";
            if(status == CUSOLVER_STATUS_EXECUTION_FAILED)
                cusolverDn_errors += "CUSOLVER_STATUS_EXECUTION_FAILED";
            if(status == CUSOLVER_STATUS_INTERNAL_ERROR)
                cusolverDn_errors += "CUSOLVER_STATUS_INTERNAL_ERROR";
            if(status == CUSOLVER_STATUS_INVALID_LICENSE)
                cusolverDn_errors += "CUSOLVER_STATUS_INVALID_LICENSE";
            if(status == CUSOLVER_STATUS_INVALID_VALUE)
                cusolverDn_errors += "CUSOLVER_STATUS_INVALID_VALUE";
            if(status == CUSOLVER_STATUS_MAPPING_ERROR)
                cusolverDn_errors += "CUSOLVER_STATUS_MAPPING_ERROR";
            if(status == CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED)
                cusolverDn_errors += "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
            if(status == CUSOLVER_STATUS_NOT_INITIALIZED)
                cusolverDn_errors += "CUSOLVER_STATUS_NOT_INITIALIZED";
            if(status == CUSOLVER_STATUS_NOT_SUPPORTED)
                cusolverDn_errors += "CUSOLVER_STATUS_NOT_SUPPORTED";
            if(status == CUSOLVER_STATUS_SUCCESS)
                cusolverDn_errors += "CUSOLVER_STATUS_SUCCESS";
            if(status == CUSOLVER_STATUS_ZERO_PIVOT)
                cusolverDn_errors += "CUSOLVER_STATUS_ZERO_PIVOT";

            if (cusolverDn_errors == " ")
                cusolverDn_errors = "unknown cusolverDn error state";
#ifdef QT_NO_DEBUG
            AssertThrow(false, dealii::ExcMessage(cusolverDn_errors.c_str() ) );
#else
            Assert(false, dealii::ExcMessage(cusolverDn_errors.c_str() ) );
#endif
        }
    }
#endif
    // @sect5{Funktion: Init}
    //
    // Initialize cusolverDn.
    // If it fails an exception is thrown. To get a more readable output
    // we use deal.II's Assert and AssertThrow macros.
    static void  Init()
    {
        cusolverStatus_t s = cusolverDnCreate(&cusolver_handle);

#ifdef DEBUG
        if (s == CUSOLVER_STATUS_SUCCESS)
            std::cout << "cusolverDn init succeeded" << std::endl;
#endif

#ifndef DEBUG
        AssertThrow(s == CUSOLVER_STATUS_SUCCESS,
                    dealii::ExcMessage("cusolverDn init failed"));
#else
        Assert(s == CUSOLVER_STATUS_SUCCESS,
               dealii::ExcMessage("cusolverDn init failed"));
#endif
    }


    // @sect5{Funktion: Shutdown}
    //
    // Shut down cusolverDn. If it fails an exception is thrown.
    static void Shutdown()
    {

        cusolverStatus_t s = cusolverDnDestroy(cusolver_handle);

#ifdef DEBUG
        if (s == CUSOLVER_STATUS_SUCCESS)
            std::cout << "cusolverDn shutdown succeeded" << std::endl;
#endif

#ifdef QT_NO_DEBUG
        AssertThrow(s == CUSOLVER_STATUS_SUCCESS,
                    dealii::ExcMessage("cusolverDn shutdown failed"));
#else
        Assert(s == CUSOLVER_STATUS_SUCCESS,
                    dealii::ExcMessage("cusolverDn shutdown failed"));
#endif
    }

    //---------------------------------------------------------
    // @sect5{SVD functions (gpu)}
    // @sect6{Function (gpu): SVD}
    // Calculate singular-value-decomposition of $m \times n$ matrix $A$ of type s, d, c, z with
    // $A=U \cdot S \cdot V^T$
    // for cublas matrices via cuSolver. Correspondingly the vector of $S$ is of type s or d.
    // \param A : input matrix, contents are overwritten during lapack function call
    // \param U : output matrix of dimension $m \times \min(m,n)$, contains left singular vectors
    // \param S : output vector of length $\min(m,n)$, contains singular values
    // \param Vt : output matrix of dimension $\min(m,n) \times n$, contains right singular vectors
    template <typename T>
    static cusolverStatus_t SVD(SciPAL::Matrix<T, BW>&                                                  A,
                                SciPAL::Matrix<T, BW>&                                                  U,
                                SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&  S,
                                SciPAL::Matrix<T, BW>&                                                  Vt)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Allocate space for U,S,Vt:
        // dimension of S can be chosen to be $\min(m,n)$ (see cuSolver documentation: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-gesvd).
        //
        // if $m<n$:
        // - U : $m \times m$ matrix (quadratic)
        // - S : $m$-dim vector
        // - Vt: $m \times n$ matrix
        //
        // if $m \geq n$:
        // - U : $m \times n$ matrix
        // - S : $n$-dim vector
        // - Vt : $n \times n$ matrix (quadratic)

        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        // Workaround for the not-so-finished implementation of NVidias SVD:
        //
        // Since the CUDA SVD function only accepts matrices with $m \geq n$ by now
        // (see remarks of the cuSolver documentation http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-gesvd)
        // we have to catch the case $m<n$ and transpose the matrix A. This leads to transposed and
        // swapped matrices U and Vt. To get the desired results these matrices have to be
        // transposed and swapped explicitly in the end.

        // Catch the case which does not work by now:
        if(m < n)
        {
            // Define temporary workspaces:
            SciPAL::Matrix<T, BW> A_tmp(n, m);
            SciPAL::Matrix<T, BW> U_tmp;
            SciPAL::Matrix<T, BW> Vt_tmp;

            // Create Adjoint(A):
            A_tmp = SciPAL::adjoint(A);

            // Calculate the SVD for Adjoint(A):
            cusolverStatus_t stat = SVD(A_tmp, U_tmp, S, Vt_tmp);

            // Determine the desired U and Vt matrices:
            U = SciPAL::adjoint<SciPAL::Matrix<T, BW> >(Vt_tmp);
            Vt = SciPAL::adjoint<SciPAL::Matrix<T, BW> >(U_tmp);

            return stat;
        }

        // Allocate device side SVD workspace:
        int Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo),sizeof(int));

        // Get work size:
        cusolverStatus_t stat = cusolverDngesvd_bufferSize(m, n, &Lwork, T(0));
#ifdef DEBUG
        check_status(stat);
#endif
        // Allocate work space:
        SciPAL::Vector<T, BW> work(Lwork);
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> Rwork(5 * min_mn);

        // Create full (but temporary) U matrix:
        SciPAL::Matrix<T, BW> U_tmp(m, m);

        // Do the SVD:
        stat = cusolverDngesvd('A', 'A', m, n, A.data(), A.leading_dim, S.data(), U_tmp.data(), U_tmp.leading_dim, Vt.data(), Vt.leading_dim,
                               work.data(), Lwork, Rwork.data(), devInfo);
        cudaDeviceSynchronize();

        // Check error information (see cuSolver documentation.):
#ifdef DEBUG
        SVD_check_devInfo(stat, devInfo);
        check_status(stat);
#endif
        // Determine the desired U:
        // Later a view can used here.
        BW::copy(m * n, U_tmp.data(), U_tmp.stride, U.data(), U.stride);

        cudaFree(devInfo);

        return stat;
    }

    // @sect6{Function (gpu): SVD (const A)}
    // Same as above, but keeps A unchanged by creating a working copy A_tmp = A.
    template <typename T>
    static cusolverStatus_t SVD(const SciPAL::Matrix<T, BW>&                                            A,
                                SciPAL::Matrix<T, BW>&                                                  U,
                                SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&  S,
                                SciPAL::Matrix<T, BW>&                                                  Vt)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return SVD(A_tmp, U, S, Vt);
    }

    //---------------------------------------------------------
    // @sect5{LUD functions (gpu)}

    // @sect6{Function (gpu): LUD}
    // Calculate LU-decomposition of $m \times n$ Matrix $A$ of type s, d, c, z with $A=P^{-1}\cdot L\cdot U$.
    //
    // $P$ is a permutation matrix.
    // $L$ is lower triangular (lower trapezoidal if $m > n$) with unit diagonal elements.
    // $U$ is upper triangular (upper trapezoidal if $m < n$).
    // \param A : input/output matrix, contains $L$ (except for its unit diagonal) and $U$ after function call.
    // \param P : output vector, describes permutation matrix as pivot indices; for $1 \leq i \leq min(m,n)$, row $i$ of the matrix was interchanged with row $P(i)$.
    template <typename T>
    static cusolverStatus_t LUD(SciPAL::Matrix<T, BW>& A,
                                SciPAL::Vector<int, BW>& P)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Allocate space for P
        P.reinit(min_mn);

        // Calculate LU factorization for A:
        int Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo),sizeof(int));

        // Get work size:
        cusolverStatus_t stat = cusolverDngetrf_bufferSize(m, n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        // Allocate work space:
        SciPAL::Vector<T, BW> work(Lwork);

        // Do the LU decomposition:
        stat = cusolverDngetrf(m, n, A.data(), A.leading_dim, work.data(), P.data(), devInfo);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif
        cudaFree(devInfo);

        return stat;
    }

    // @sect6{Function (gpu): LUD (L,U seperate)}
    // Same as above but explicitly copies the content of $A$ into matrices $L$ and $U$.
    // \param A : input/output matrix, contains $L$ (except for its unit diagonal) and $U$ after function call.
    // \param P : output vector, describes permutation matrix as pivot indices; for $1 \leq i \leq min(m,n)$, row $i$ of the matrix was interchanged with row $P(i)$.
    // \param L : output matrix; lower triangular (lower trapezoidal if $m > n$) of dimension $m \times min(m,n)$ with unit diagonal elements.
    // \param U : output matrix; upper triangular (upper trapezoidal if $m < n$) of dimension $min(m,n) \times n$.
    template <typename T>
    static cusolverStatus_t LUD(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Vector<int, BW>&        P,
                                SciPAL::Matrix<T, BW>&          L,
                                SciPAL::Matrix<T, BW>&          U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Calculate LU factorization for A:
        cusolverStatus_t stat = LUD(A, P);

        //Resize output:
        L.reinit(m, min_mn);
        U.reinit(min_mn, n);

        // Set the diagonal entries of L to 1.0,
        // copy strict lower triangular part of A into L, and
        // copy upper triangular part of A into U:
        lower_upper_triangle(A, L, U);

        return stat;
    }

    // @sect6{Function (gpu): LUD (const A, L/U seperate)}
    // Same as above, but keeps $A$ unchanged by creating a working copy A_tmp = A.
    template <typename T>
    static cusolverStatus_t LUD(const SciPAL::Matrix<T, BW>&    A,
                                SciPAL::Vector<int, BW>&        P,
                                SciPAL::Matrix<T, BW>&          L,
                                SciPAL::Matrix<T, BW>&          U)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return LUD(A_tmp, P, L, U);
    }

    //---------------------------------------------------------
    // @sect5{QRF functions (gpu)}

    // @sect6{Function (gpu): QRF}
    // Calculate QR-factorization of $m \times n$ Matrix $A$ of type s, d, c, z with $A=Q\cdot R$.

    // $Q$ orthogonal matrix of dimension $m \times m$.
    // $R$ is upper triangular matrix of dimension $m \times n$ (upper trapezoidal if $m < n$).

    // Lapack doc: 'The matrix $Q$ is represented as a product of elementary reflectors $Q = H(1) H(2) ... H(k)$ , where $k = min(m,n)$.
    // Each $H(i)$ has the form $H(i) = I - TAU \cdot v \cdot v^T$ where $TAU$ is a real scalar, and $v$ is a real vector with $v(1:i-1) = 0$ and $v(i) = 1; v(i+1:m)$ is stored on exit in $A(i+1:m,i)$, and $TAU$ in $TAU(i)$.'

    // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
    // \param TAU : output vector of length $min(m,n)$ contains the scalar factors of the elementary reflector.
    //
    // work and Lwork are referenced in order to be reused in xxmqr.
    template <typename T>
    static cusolverStatus_t QRF(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Vector<T, BW>&          tau,
                                SciPAL::Vector<T, BW>&          work = SciPAL::Vector<T, BW>(),
                                int&                            Lwork = 0)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Allocate space for tau:
        tau.reinit(min_mn);

        Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // Get work size:
        cusolverStatus_t stat = cusolverDngeqrf_bufferSize(m, n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        // Allocate work space:
        work.reinit(Lwork);

        // Do the QR decomposition:
        stat = cusolverDngeqrf(m, n, A.data(), A.leading_dim, tau.data(), work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        cudaFree(devInfo);
        return stat;
    }

    // @sect6{Function (gpu): QRF (Q,R seperate)}
    // Same as above but explicitly processec the output content of $A$ into matrices $Q$ and $R$.
    // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
    // \param Q : output vector, orthogonal matrix of dimension $m \times min(m,n)$.
    // \param R : output matrix; upper triangular (upper trapezoidal if $m < n$) of dimension $min(m,n) \times n$.
    template <typename T>
    static cusolverStatus_t QRF(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Matrix<T, BW>&          Q,
                                SciPAL::Matrix<T, BW>&          R)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Allocate space for Q, R:
        Q.reinit(m, min_mn);
        R.reinit(min_mn, n);

        // Allocate work space:
        SciPAL::Vector<T, BW> tau(1);
        SciPAL::Vector<T, BW> work(1);
        int Lwork;

        // Do the QR decomposition:
        cusolverStatus_t stat = QRF(A, tau, work, Lwork);

        // Create an unit matrix:
        //
        // TO DO: This way is quite unefficient. Thus, find a better way. Deal.II IdentityMatrix did not work.
        for (unsigned int i = 0; i < min_mn; i++)
            Q(i, i, T(1));

        // Copy upper triangle part of A to R:
        for(unsigned int i = 0; i < n; i++)
            BW::copy(((i+1<min_mn)?i+1:min_mn), A.data() + i * A.leading_dim, A.stride, R.data() + i * R.leading_dim, R.stride);

        cudaDeviceSynchronize();

        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // Compute Q from the Householder vectors stored in the (strict) lower triangle part of A and tau:
        // It seems, as if the parameters m and n are not the rows/columns of A but of C (or in this case Q). This is not what is written in the CUDA-documentation: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-ormqr
        stat = cusolverDnormqr(CUBLAS_SIDE_LEFT, CUBLAS_OP_N, Q.n_rows(), Q.n_cols(), min_mn, A.data(), A.leading_dim, tau.data(), Q.data(), Q.leading_dim, work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        cudaFree(devInfo);
        return stat;
    }

    // @sect6{Function (gpu): QRF (const A, Q/R seperate)}
    // Same as above, but keeps $A$ unchanged by using a working copy A_tmp = A
    template <typename T>
    static cusolverStatus_t QRF(const SciPAL::Matrix<T, BW>&    A,
                                SciPAL::Matrix<T, BW>&          Q,
                                SciPAL::Matrix<T, BW>&          R)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return QR(A_tmp, Q, R);
    }

    //---------------------------------------------------------
    // @sect5{LDL functions (gpu)}
    // @sect6{Function (gpu): LDL}
    // Calculates the LDL decomposition for a $n \times n$ matrix A with $P \cdot A \cdot P^T = L \cdot D \cdot L^T$ or
    // $P \cdot A \cdot P^T = L^T \cdot D \cdot L = U \cdot D \cdot U^T$ for cublas matrices via cuSolver.
    // It seems that the cuSolver functions cusolverDn<X>sytrf do not work properly. Multiplying
    // the output leads to a results different to the input. Possilby there is a mistake in the
    // the wrapper function.
    // TO DO: repair this wrapper function or wait until cusolverDn<X>sytrf works properly.
    template <typename T>
    static cusolverStatus_t LDL(SciPAL::Matrix<T, BW>&      A,
                                SciPAL::Vector<int, BW>&    P,
                                SciPAL::Matrix<T, BW>&      L,
                                SciPAL::Matrix<T, BW>&      D,
                                cublasFillMode_t            uplo = CUBLAS_FILL_MODE_LOWER)
    {
        // dimension of A
        unsigned int n = A.n_cols();
#ifdef DEBUG
        unsigned int m = A.n_rows();
#ifdef QT_NO_DEBUG
            AssertThrow(m==n, dealii::ExcMessage("LDL decomposition: matrix not quadratic"));
#else
            Assert(m==n, dealii::ExcMessage("LDL decomposition: matrix not quadratic"));
#endif
#endif

        // Allocate space for P,L,D:
        P.reinit(n);
        L.reinit(n, n);
        D.reinit(n, n);

        int Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // Get work size:
        cusolverStatus_t stat = cusolverDnsytrf_bufferSize(n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        // Allocate work space:
        SciPAL::Vector<T, BW> work(Lwork);

        // Do the LDL decomposition:
        stat = cusolverDnsytrf(uplo, n, A.data(), A.leading_dim, P.data(), work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
#ifdef DEBUG
        check_status(stat);
#endif

        cudaFree(devInfo);

        BW::copy(n, A.data(), n + 1, D.data(), n + 1);
        if(uplo == CUBLAS_FILL_MODE_LOWER)
        {
            for(unsigned int i = 0; i < n; i++)
            {
                L(i, i, 1.0);
                if(i < n - 1)
                    BW::copy(n - i - 1, A.data() + i * n + i + 1, 1, L.data() + i * n + i + 1, 1);
            }
        }
        else
        {
#ifdef DEBUG
#ifdef QT_NO_DEBUG
            AssertThrow(false, dealii::ExcNotImplemented());
#else
            Assert(false, dealii::ExcNotImplemented());
#endif
#endif
        }
        cudaDeviceSynchronize();
        return stat;
    }

private:
    //---------------------------------------------------------
    // @sect5{Copy functions}
    // @sect6{Function: lower_upper_triangle}
    // Helper function. (TO DO: Move it to the Matrix class?) Copies lower and upper triangular
    // (trapezoidal) parts of $A$ into $L$ and $U$ using cublas copy.
    // TO DO: Does not yet supports leading dimension nor stride.
    // \param A : input matrix of dimension $m \times n$.
    // \param L : output matrix of dimension $m \times n$, contains strict lower triangular (trapezoidal) values of $A$ and unit diagonal.
    // \param U : output matrix of dimension $m \times n$, contains upper triangular (trapezoidal) values of $A$.
    template <typename T>
    static void lower_upper_triangle(const SciPAL::Matrix<T, BW>& A, SciPAL::Matrix<T, BW>& L, SciPAL::Matrix<T, BW>& U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // Set diagonal entries of L to 1.0 (TO DO: This way is quite unefficient.),
        // copy strict lower triangular part of A into L, and
        // copy upper triangular part of A into U.
        for(unsigned int i = 0, ip1 = 1, im = 0; i < min_mn - 1; i++, ip1++, im += m)
        {
            unsigned int impip1 = im + ip1;
            L(i, i, 1.0);                   // unefficient
            BW::copy(m - ip1, A.data() + impip1, 1, L.data() + impip1    , 1);
            BW::copy(    ip1, A.data() + im    , 1, U.data() + i * min_mn, 1);
        }
        L(min_mn - 1, min_mn - 1, 1.0);
        if(n >= m)
        {
            BW::copy(m * (n - m + 1), A.data() + m * (m - 1), 1, U.data() + m * (m - 1), 1);
        }
        else
        {
            BW::copy(m - n, A.data() +m * (n - 1) + n, 1, L.data() + m * (n - 1) + n, 1);
            BW::copy(n, A.data() + m * (n - 1), 1, U.data() + n * (n - 1), 1);
        }
        cudaDeviceSynchronize();
    }

    //---------------------------------------------------------
    // @sect5{cusolverDn function wrappers}
    // @sect6{Wrapper: SVD cusolverDn-wrapper}
    // Wraps the general matrix svd cusolverDn function for the four types s, d, c and z.
    // cuSolver doc: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-gesvd
    inline static cusolverStatus_t cusolverDngesvd(char jobu, char jobvt, int m, int n, float *A, int lda, float *S, float *U, int ldu, float *VT, int ldvt, float *Work, int Lwork, float *rwork, int *devInfo)
    {
        return cusolverDnSgesvd(cusolver_handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, Work, Lwork, rwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngesvd(char jobu, char jobvt, int m, int n, double *A, int lda, double *S, double *U, int ldu, double *VT, int ldvt, double *Work, int Lwork, double *rwork, int *devInfo)
    {
        return cusolverDnDgesvd(cusolver_handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, Work, Lwork, rwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngesvd(char jobu, char jobvt, int m, int n, SciPAL::CudaComplex<float> *A, int lda, float *S, SciPAL::CudaComplex<float> *U, int ldu, SciPAL::CudaComplex<float> *VT, int ldvt, SciPAL::CudaComplex<float> *Work, int Lwork, float *rwork, int *devInfo)
    {
        return cusolverDnCgesvd(cusolver_handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, Work, Lwork, rwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngesvd(char jobu, char jobvt, int m, int n, SciPAL::CudaComplex<double> *A, int lda, double *S, SciPAL::CudaComplex<double> *U, int ldu, SciPAL::CudaComplex<double> *VT, int ldvt, SciPAL::CudaComplex<double> *Work, int Lwork, double *rwork, int *devInfo)
    {
        return cusolverDnZgesvd(cusolver_handle, jobu, jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, Work, Lwork, rwork, devInfo);
    }

    inline static cusolverStatus_t cusolverDngesvd_bufferSize(int m, int n, int *Lwork, float /*Zero*/)
    {
        return cusolverDnSgesvd_bufferSize(cusolver_handle, m, n, Lwork);
    }
    inline static cusolverStatus_t cusolverDngesvd_bufferSize(int m, int n, int *Lwork, double /*Zero*/)
    {
        return cusolverDnDgesvd_bufferSize(cusolver_handle, m, n, Lwork);
    }
    inline static cusolverStatus_t cusolverDngesvd_bufferSize(int m, int n, int *Lwork, SciPAL::CudaComplex<float> /*Zero*/)
    {
        return cusolverDnCgesvd_bufferSize(cusolver_handle, m, n, Lwork);
    }
    inline static cusolverStatus_t cusolverDngesvd_bufferSize(int m, int n, int *Lwork, SciPAL::CudaComplex<double> /*Zero*/)
    {
        return cusolverDnZgesvd_bufferSize(cusolver_handle, m, n, Lwork);
    }
#ifdef DEBUG
    inline static void SVD_check_devInfo(cusolverStatus_t stat, int *devInfo)
    {
        // check error information (see cuSolver documentation)
        if(stat != CUSOLVER_STATUS_SUCCESS){
            std::cout << "cusolver crashed" << std::endl;
            std::cout << "device information:" << std::endl;
            int Info;
            cudaMemcpy(&Info, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
            if(Info == 0)
                std::cout << "\tThe operation was successful" << std::endl;
            if(Info > 0)
                std::cout << "\t" << Info << " superdiagonal(s) did not converge to zero" << std::endl;
            if(Info < 0)
                std::cout << "\tthe " << -Info << "-th parameter is wrong" << std::endl;
        }
    }
#endif
    //---------------------------------------------------------
    // @sect6{Wrapper: LUD cusolverDn-wrapper}
    // Wraps the general matrix LU decomposition cusolverDn function for the four types s, d, c and z.
    // cuSolver doc: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-getrf
    inline static cusolverStatus_t cusolverDngetrf(int m, int n, float *A, int lda, float *Workspace, int *devIpiv, int *devInfo)
    {
        return cusolverDnSgetrf(cusolver_handle, m, n, A, lda, Workspace, devIpiv, devInfo);
    }
    inline static cusolverStatus_t cusolverDngetrf(int m, int n, double *A, int lda, double *Workspace, int *devIpiv, int *devInfo)
    {
        return cusolverDnDgetrf(cusolver_handle, m, n, A, lda, Workspace, devIpiv, devInfo);
    }
    inline static cusolverStatus_t cusolverDngetrf(int m, int n, SciPAL::CudaComplex<float> *A, int lda, SciPAL::CudaComplex<float> *Workspace, int *devIpiv, int *devInfo)
    {
        return cusolverDnCgetrf(cusolver_handle, m, n, A, lda, Workspace, devIpiv, devInfo);
    }
    inline static cusolverStatus_t cusolverDngetrf(int m, int n, SciPAL::CudaComplex<double> *A, int lda, SciPAL::CudaComplex<double> *Workspace, int *devIpiv, int *devInfo)
    {
        return cusolverDnZgetrf(cusolver_handle, m, n, A, lda, Workspace, devIpiv, devInfo);
    }

    inline static cusolverStatus_t cusolverDngetrf_bufferSize(int m, int n, float *A, int lda, int *Lwork)
    {
        return cusolverDnSgetrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngetrf_bufferSize(int m, int n, double *A, int lda, int *Lwork)
    {
        return cusolverDnDgetrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngetrf_bufferSize(int m, int n, SciPAL::CudaComplex<float> *A, int lda, int *Lwork)
    {
        return cusolverDnCgetrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngetrf_bufferSize(int m, int n, SciPAL::CudaComplex<double> *A, int lda, int *Lwork)
    {
        return cusolverDnZgetrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }

    //---------------------------------------------------------
    // @sect6{Wrapper: QRF cusolverDn-wrapper}
    // Wraps the general matrix QR factorization cusolverDn function for the four types s, d, c and z.
    // cuSolver doc: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-geqrf
    inline static cusolverStatus_t cusolverDngeqrf(int m, int n, float *A, int lda, float *TAU, float *Workspace, int Lwork, int *devInfo)
    {
        return cusolverDnSgeqrf(cusolver_handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngeqrf(int m, int n, double *A, int lda, double *TAU, double *Workspace, int Lwork, int *devInfo)
    {
        return cusolverDnDgeqrf(cusolver_handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngeqrf(int m, int n, SciPAL::CudaComplex<float> *A, int lda, SciPAL::CudaComplex<float> *TAU, SciPAL::CudaComplex<float> *Workspace, int Lwork, int *devInfo)
    {
        return cusolverDnCgeqrf(cusolver_handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDngeqrf(int m, int n, SciPAL::CudaComplex<double> *A, int lda, SciPAL::CudaComplex<double> *TAU, SciPAL::CudaComplex<double> *Workspace, int Lwork, int *devInfo)
    {
        return cusolverDnZgeqrf(cusolver_handle, m, n, A, lda, TAU, Workspace, Lwork, devInfo);
    }

    inline static cusolverStatus_t cusolverDngeqrf_bufferSize(int m, int n, float *A, int lda, int *Lwork)
    {
        return cusolverDnSgeqrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngeqrf_bufferSize(int m, int n, double *A, int lda, int *Lwork)
    {
        return cusolverDnDgeqrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngeqrf_bufferSize(int m, int n, SciPAL::CudaComplex<float> *A, int lda, int *Lwork)
    {
        return cusolverDnCgeqrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDngeqrf_bufferSize(int m, int n, SciPAL::CudaComplex<double> *A, int lda, int *Lwork)
    {
        return cusolverDnZgeqrf_bufferSize(cusolver_handle, m, n, A, lda, Lwork);
    }

    inline static cusolverStatus_t cusolverDnormqr(cublasSideMode_t side, cublasOperation_t trans, int m, int n, int k, const float *A, int lda, const float *tau, float *C, int ldc, float *work, int lwork, int *devInfo)
    {
        return cusolverDnSormqr(cusolver_handle, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnormqr(cublasSideMode_t side, cublasOperation_t trans, int m, int n, int k, const double *A, int lda, const double *tau, double *C, int ldc, double *work, int lwork, int *devInfo)
    {
        return cusolverDnDormqr(cusolver_handle, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnormqr(cublasSideMode_t side, cublasOperation_t trans, int m, int n, int k, const SciPAL::CudaComplex<float> *A, int lda, const SciPAL::CudaComplex<float> *tau, SciPAL::CudaComplex<float> *C, int ldc, SciPAL::CudaComplex<float> *work, int lwork, int *devInfo)
    {
        return cusolverDnCunmqr(cusolver_handle, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnormqr(cublasSideMode_t side, cublasOperation_t trans, int m, int n, int k, const SciPAL::CudaComplex<double> *A, int lda, const SciPAL::CudaComplex<double> *tau, SciPAL::CudaComplex<double> *C, int ldc, SciPAL::CudaComplex<double> *work, int lwork, int *devInfo)
    {
        return cusolverDnZunmqr(cusolver_handle, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, devInfo);
    }

    //---------------------------------------------------------
    // @sect6{Wrapper: LDL cusolverDn-wrapper}
    // Wraps the general matrix LDL decomposition cusolverDn function for the four types s, d, c and z.
    // cuSolver doc: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-sytrf
    inline static cusolverStatus_t cusolverDnsytrf(cublasFillMode_t uplo, int n, float *A, int lda, int *ipiv, float *work, int lwork, int *devInfo)
    {
        return cusolverDnSsytrf(cusolver_handle, uplo, n, A, lda, ipiv, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnsytrf(cublasFillMode_t uplo, int n, double *A, int lda, int *ipiv, double *work, int lwork, int *devInfo)
    {
        return cusolverDnDsytrf(cusolver_handle, uplo, n, A, lda, ipiv, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnsytrf(cublasFillMode_t uplo, int n, SciPAL::CudaComplex<float> *A, int lda, int *ipiv, SciPAL::CudaComplex<float> *work, int lwork, int *devInfo)
    {
        return cusolverDnCsytrf(cusolver_handle, uplo, n, A, lda, ipiv, work, lwork, devInfo);
    }
    inline static cusolverStatus_t cusolverDnsytrf(cublasFillMode_t uplo, int n, SciPAL::CudaComplex<double> *A, int lda, int *ipiv, SciPAL::CudaComplex<double> *work, int lwork, int *devInfo)
    {
        return cusolverDnZsytrf(cusolver_handle, uplo, n, A, lda, ipiv, work, lwork, devInfo);
    }

    inline static cusolverStatus_t cusolverDnsytrf_bufferSize(int n, float *A, int lda, int *Lwork)
    {
        return cusolverDnSsytrf_bufferSize(cusolver_handle, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDnsytrf_bufferSize(int n, double *A, int lda, int *Lwork)
    {
        return cusolverDnDsytrf_bufferSize(cusolver_handle, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDnsytrf_bufferSize(int n, SciPAL::CudaComplex<float> *A, int lda, int *Lwork)
    {
        return cusolverDnCsytrf_bufferSize(cusolver_handle, n, A, lda, Lwork);
    }
    inline static cusolverStatus_t cusolverDnsytrf_bufferSize(int n, SciPAL::CudaComplex<double> *A, int lda, int *Lwork)
    {
        return cusolverDnZsytrf_bufferSize(cusolver_handle, n, A, lda, Lwork);
    }
};  //END struct cusolverDn
} //END namespace SciPAL
#endif // CUSOLVERDN_WRAPPER_STEP43_HH
