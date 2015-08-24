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

// @sect3{struct cusolverDn: cusolverDn S,D,C,Z wrapper functions}
//!
//! We remove the type dependence of the function of cusolverDn by using C++'s polymorphism and templates.
//! (...)
//! For the numerical purpose of the functions have a look at the cusolverDn doc.
//!

struct cusolverDn
{

    // @sect4{Funktion: name}
    //!
    //! Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "cusolverDn"; }

    //! Compile-time variable for identifying the side of PCIe bus this solver works on.
    static const ParallelArch arch = gpu_cuda;

    typedef typename archTraits<arch>::BlasType BW;

    // @sect4{Funktion: check_status}
    //!
    //! Evaluate the @p cusolverStatus argument and throws
    //! an exception in case an error occurs.
    static void  check_status(const cusolverStatus_t & status)
    {
        //! The following is for tracking cusolverDn usage
        //! std::cout << "Checking cusolver status" << std::endl;
#ifdef DEBUG
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
#endif
    }

    // @sect4{Funktion: Init}
    //!
    //! Initialize cusolverDn.
    //! If it fails an exception is thrown. To get a more readable output
    //! we use deal.II's Assert and AssertThrow macros.
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


    // @sect4{Funktion: Shutdown}
    //!
    //! Shut down cusolverDn. If it fails an exception is thrown.
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
    //SVD functions:
    // Calculates the SVD for a Matrix A (m x n) with A=U*S*Vt for cublas matrices via cuSolver.
    // The content of A is destroyed during the process.
    // A is of matrix type s,d,c or z.  Correspondingly the vector of S is of type s or d.
    // (...)
    template <typename T>
    static cusolverStatus_t SVD(SciPAL::Matrix<T, BW>&                                                  A,
                                SciPAL::Matrix<T, BW>&                                                  U,
                                SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&  S,
                                SciPAL::Matrix<T, BW>&                                                  Vt)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // allocate space for U,S,Vt:
        // dimension of S can be chosen to be min(m,n) (see cuSolver doc.)
        // if m<n:  U   m x m matrix (quadratic)
        //          S   m-dim vector
        //          Vt  m x n matrix
        // if m>=n: U   m x n matrix (quadratic)
        //          S   n-dim vector
        //          Vt  n x n matrix
        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        // Workaround for the not-so-finished implementation of Nvidias SVD:
        // CUDA svd function only accepts matrices with m>=n by now. (see cuSolver doc.)
        // Hence, we first have to transpose a matrix A with m<n which leads to transposed and
        // swapped matrices U and Vt. To get the desired results these have to be transposed
        // and swapped explicitly in the end.

        // catch the case which does not work by now
        if(m < n)
        {
            // define temporary workspaces:
            SciPAL::Matrix<T, BW> A_tmp(n, m);
            SciPAL::Matrix<T, BW> U_tmp;
            SciPAL::Matrix<T, BW> Vt_tmp;

            // create Adjoint(A):
            A_tmp = SciPAL::adjoint(A);

            // calculate SVD for Adjoint(A):
            cusolverStatus_t stat = SVD(A_tmp, U_tmp, S, Vt_tmp);

            // determine the desired U and Vt matrices:
            U = SciPAL::adjoint<SciPAL::Matrix<T, BW> >(Vt_tmp);

            BW::copy(m * n, U_tmp.data(), U_tmp.stride, A_tmp.data(), A_tmp.stride);                  // recycling A_tmp

            Vt = SciPAL::adjoint<SciPAL::Matrix<T, BW> >(A_tmp);

            return stat;
        }

        // --- device side SVD workspace and matrices
        int Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo),sizeof(int));

        // get work size
        cusolverStatus_t stat = cusolverDngesvd_bufferSize(m, n, &Lwork, T(0));
        check_status(stat);

        // allocate work space
        SciPAL::Vector<T, BW> work(Lwork);

        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> Rwork(5 * min_mn);

        // create full (but temporary) U matrix:
        SciPAL::Matrix<T, BW> U_tmp(m, m);

        // do the SVD
        stat = cusolverDngesvd('A', 'A', m, n, A.data(), A.leading_dim, S.data(), U_tmp.data(), U_tmp.leading_dim, Vt.data(), Vt.leading_dim,
                               work.data(), Lwork, Rwork.data(), devInfo);
        cudaDeviceSynchronize();

        // check error information (see cusolver doc.)
        SVD_check_devInfo(stat, devInfo);
        check_status(stat);

        // determine the desired U
        // Later a view can used here.
        BW::copy(m * n, U_tmp.data(), U_tmp.stride, U.data(), U.stride);

        cudaFree(devInfo);

        return stat;
    }
    // Calculates the SVD for a Matrix A (m x n) with A=U*S*Vt for cublas matrices via cuSolver.
    // A is of matrix type s,d,c or z.  Correspondingly the vector of S is of type s or d.
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
    //LUD functions:
    // Calculates the LU decomposition for a Matrix A (m x n) with P*A=L*U for cublas matrices via cuSolver.
    // The result is stored in A.
    // A is of matrix type s,d,c or z.
    template <typename T>
    static cusolverStatus_t LUD(SciPAL::Matrix<T, BW>& A,
                                SciPAL::Vector<int, BW>& P)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // allocate space for P
        P.reinit(min_mn);

        // calculate LU factorization for A:
        int Lwork = 0;
        int *devInfo;

        cudaMalloc((void**)(&devInfo),sizeof(int));

        // get work size
        cusolverStatus_t stat = cusolverDngetrf_bufferSize(m, n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
        check_status(stat);

        // allocate work space
        SciPAL::Vector<T, BW> work(Lwork);

        // do the LU decomp
        stat = cusolverDngetrf(m, n, A.data(), A.leading_dim, work.data(), P.data(), devInfo);
        cudaDeviceSynchronize();
        check_status(stat);
        cudaFree(devInfo);

        return stat;
    }
    // Calculates the LU decomposition for a Matrix A (m x n) with P*A=L*U for cublas matrices via cuSolver.
    // The content of A is destroyed during the process.
    // A is of matrix type s,d,c or z.
    template <typename T>
    static cusolverStatus_t LUD(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Vector<int, BW>&        P,
                                SciPAL::Matrix<T, BW>&          L,
                                SciPAL::Matrix<T, BW>&          U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // calculate LU factorization for A:
        cusolverStatus_t stat = LUD(A, P);

        // allocate space for L,U
        L.reinit(m, min_mn);
        U.reinit(min_mn, n);

        // set the diagonal entries of L to 1.0
        // copy strict lower triangular part of A into L
        // copy upper triangular part of A into U
        lower_upper_triangle(A, L, U);

        return stat;
    }
    // Calculates the LU decomposition for a Matrix A (m x n) with P*A=L*U for cublas matrices via cuSolver.
    // A is of matrix type s,d,c or z.
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
    //QRF functions:
    // Calculates the QR decomposition for a Matrix A (m x n) with A=Q*R for cublas matrices via cuSolver.
    // The result is stored in A.
    // A is of matrix type s,d,c or z.
    // work and Lwork are referenced in order to be reused in xxmqr
    template <typename T>
    static cusolverStatus_t QRF(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Vector<T, BW>&          tau,
                                SciPAL::Vector<T, BW>&          work = SciPAL::Vector<T, BW>(),
                                int&                            Lwork = 0)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // allocate space for tau
        tau.reinit(min_mn);

        Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // get work size
        cusolverStatus_t stat = cusolverDngeqrf_bufferSize(m, n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
        check_status(stat);

        // allocate work space
        work.reinit(Lwork);

        // do the QR decomp
        stat = cusolverDngeqrf(m, n, A.data(), A.leading_dim, tau.data(), work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
        check_status(stat);

        cudaFree(devInfo);
        return stat;
    }

    // Calculates the QR decomposition for a Matrix A (m x n) with A=Q*R for cublas matrices via cuSolver.
    // The content of A is destroyed during the process.
    // A is of matrix type s,d,c or z.
    template <typename T>
    static cusolverStatus_t QRF(SciPAL::Matrix<T, BW>&          A,
                                SciPAL::Matrix<T, BW>&          Q,
                                SciPAL::Matrix<T, BW>&          R)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        // allocate space for Q,R
        //Q.reinit(m, min_mn);
        //R.reinit(min_mn, n);
        Q.reinit(m, m);
        R.reinit(m, n);

        SciPAL::Vector<T, BW> tau(1);
        SciPAL::Vector<T, BW> work(1);
        int Lwork;

        cusolverStatus_t stat = QRF(A, tau, work, Lwork);

        //for (unsigned int i = 0; i < min_mn; i++)
        //    Q(i, i, 1.0);
        for (unsigned int i = 0; i < m; i++)
            Q(i, i, 1.0);

        // copy upper triangle part of A to R
        //for(unsigned int i = 0; i < n; i++)
        //    BW::copy(((i+1<min_mn)?i+1:min_mn), A.data() + i * A.leading_dim, A.stride, R.data() + i * R.leading_dim, R.stride);
        for(unsigned int i = 0; i < n; i++)
            BW::copy(((i+1<m)?i+1:m), A.data() + i * A.leading_dim, A.stride, R.data() + i * R.leading_dim, R.stride);

        cudaDeviceSynchronize();

        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // compute Q from the Householder vectors stored in tau and the (strict) lower triangle part of A
        // It seems, as if the parameters m and n are not the rows/columns of A but of C (or in this case Q). This is not what is written in the CUDA-documentation: http://docs.nvidia.com/cuda/cusolver/index.html#cuds-lt-t-gt-ormqr
        stat = cusolverDnormqr(CUBLAS_SIDE_LEFT, CUBLAS_OP_N, Q.n_rows(), Q.n_cols(), min_mn, A.data(), A.leading_dim, tau.data(), Q.data(), Q.leading_dim, work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
        check_status(stat);

        cudaFree(devInfo);
        return stat;
    }
    // Calculates the QR decomposition for a Matrix A (m x n) with A=Q*R for cublas matrices via cuSolver.
    // A is of matrix type s,d,c or z.
    template <typename T>
    static cusolverStatus_t QRF(const SciPAL::Matrix<T, BW>&    A,
                                SciPAL::Matrix<T, BW>&          Q,
                                SciPAL::Matrix<T, BW>&          R)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return QR(A_tmp, Q, R);
    }

    //---------------------------------------------------------
    // LDL functions:
    // Calculates the LDL decomposition for a Matrix A (n x n) with P*A*Pt=L*D*Lt or
    // P*A*Pt=Lt*D*L=U*D*Ut for cublas matrices via cuSolver.
    // It seems that the cusolver functions cusolverDn<X>sytrf do not work properly. Multiplying
    // the output leads to a results different to the input. Possilby there is a mistake in the
    // the wrapper function.
    // TODO: repair this wrapper function or wait until cusolverDn<X>sytrf works properly
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

        // allocate space for P,L,D
        P.reinit(n);
        L.reinit(n, n);
        D.reinit(n, n);

        // calculate QR factorization for A:
        int Lwork = 0;
        int *devInfo;
        cudaMalloc((void**)(&devInfo), sizeof(int));

        // get work size
        cusolverStatus_t stat = cusolverDnsytrf_bufferSize(n, A.data(), A.leading_dim, &Lwork);
        cudaDeviceSynchronize();
        check_status(stat);

        // allocate work space
        SciPAL::Vector<T, BW> work(Lwork);

        // do the LDL decomp
        stat = cusolverDnsytrf(uplo, n, A.data(), A.leading_dim, P.data(), work.data(), Lwork, devInfo);
        cudaDeviceSynchronize();
        check_status(stat);

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
    //Helper functions (move them to the Matrix class?):
    // Does not yet supports leading dimension nor stride
    template <typename T>
    static void lower_upper_triangle(const SciPAL::Matrix<T, BW>& A, SciPAL::Matrix<T, BW>& L, SciPAL::Matrix<T, BW>& U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);
        // set diagonal entries of L to 1.0
        // copy strict lower triangular part of A into L
        // copy upper triangular part of A into U
        for(unsigned int i = 0, ip1 = 1, im = 0; i < min_mn - 1; i++, ip1++, im += m)
        {
            unsigned int impip1 = im + ip1;
            L(i, i, 1.0);
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
    //SVD cusolverDn-wrapper:
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

    inline static void SVD_check_devInfo(cusolverStatus_t stat, int *devInfo)
    {
        // check error information (see cusolver doc.)
#ifdef DEBUG
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
#endif
    }

    //---------------------------------------------------------
    //LUD cusolverDn-wrapper:
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
    //QRF cusolverDn-wrapper:
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
    //LDL cusolverDn-wrapper:
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
