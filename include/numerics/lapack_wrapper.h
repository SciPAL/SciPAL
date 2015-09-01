#ifndef LAPACK_WRAPPER_HH
#define LAPACK_WRAPPER_HH

// This header includes wrapper for SVD, LUD and (partial) QRF.
// The methods for QRF who returns the Q and the R are not finished yet!

//STL includes
#include <string>

//BASE includes
#include <base/VTraits.h>
#include <base/ParallelArch.h>
#include <base/ArchTraits.h>

//...
#include <deal.II/base/timer.h>

// LAPACK
extern "C" {
//diagonalization of a tridiagonal matrix
void dstev_(char& jobz, int& size, double& d, double& e, double& z,
            int& ldz, double* work, int& info);

//SVD routines:
void sgesvd_(
        char            * JOBU,
        char            * JOBVT,
        const size_t    * M,
        const size_t    * N,
        float           * A,
        const size_t    * LDA,
        float           * S,
        float           * U,
        const size_t    * LDU,
        float           * VT,
        const size_t    * LDVT,
        float           * WORK,
        const size_t    * LWORK,
        int             * INFO);
void dgesvd_(
        char            * JOBU,
        char            * JOBVT,
        const size_t    * M,
        const size_t    * N,
        double          * A,
        const size_t    * LDA,
        double          * S,
        double          * U,
        const size_t    * LDU,
        double          * VT,
        const size_t    * LDVT,
        double          * WORK,
        const size_t    * LWORK,
        int             * INFO);
void cgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        float                       * S,
        SciPAL::CudaComplex<float>  * U,
        const size_t                * LDU,
        SciPAL::CudaComplex<float>  * VT,
        const size_t                * LDVT,
        SciPAL::CudaComplex<float>  * WORK,
        const size_t                * LWORK,
        float                       * RWORK,
        int                         * INFO);
void zgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        double                      * S,
        SciPAL::CudaComplex<double> * U,
        const size_t                * LDU,
        SciPAL::CudaComplex<double> * VT,
        const size_t                * LDVT,
        SciPAL::CudaComplex<double> * WORK,
        const size_t                * LWORK,
        double                      * RWORK,
        int                         * INFO);

//LUD routines:
void sgetrf_(
        const size_t    * M,
        const size_t    * N,
        float           * A,
        const size_t    * LDA,
        int             * IPIV,
        int             * INFO);
void dgetrf_(
        const size_t    * M,
        const size_t    * N,
        double          * A,
        const size_t    * LDA,
        int             * IPIV,
        int             * INFO);
void cgetrf_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        int                         * IPIV,
        int                         * INFO);
void zgetrf_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        int                         * IPIV,
        int                         * INFO);

//QRF routines:
void sgeqrf_(
        const size_t    * M,
        const size_t    * N,
        float           * A,
        const size_t    * LDA,
        float           * TAU,
        float           * WORK,
        const size_t    * LWORK,
        int             * INFO);
void dgeqrf_(
        const size_t    * M,
        const size_t    * N,
        double          * A,
        const size_t    * LDA,
        double          * TAU,
        double          * WORK,
        const size_t    * LWORK,
        int             * INFO);
void cgeqrf_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        SciPAL::CudaComplex<float>  * TAU,
        SciPAL::CudaComplex<float>  * WORK,
        const size_t                * LWORK,
        int                         * INFO);
void zgeqrf_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        SciPAL::CudaComplex<double> * TAU,
        SciPAL::CudaComplex<double> * WORK,
        const size_t                * LWORK,
        int                         * INFO);

//MQR routines:
void sormqr_(
        char            * SIDE,
        char            * TRANS,
        const size_t    * M,
        const size_t    * N,
        const size_t    * K,
        float           * A,
        const size_t    * LDA,
        float           * TAU,
        float           * C,
        const size_t    * LDC,
        float           * WORK,
        const size_t    * LWORK,
        int             * INFO);
void dormqr_(
        char            * SIDE,
        char            * TRANS,
        const size_t    * M,
        const size_t    * N,
        const size_t    * K,
        double          * A,
        const size_t    * LDA,
        double          * TAU,
        double          * C,
        const size_t    * LDC,
        double          * WORK,
        const size_t    * LWORK,
        int             * INFO);
void cunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        const size_t                * M,
        const size_t                * N,
        const size_t                * K,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        SciPAL::CudaComplex<float>  * TAU,
        SciPAL::CudaComplex<float>  * C,
        const size_t                * LDC,
        SciPAL::CudaComplex<float>  * WORK,
        const size_t                * LWORK,
        int                         * INFO);
void zunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        const size_t                * M,
        const size_t                * N,
        const size_t                * K,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        SciPAL::CudaComplex<double> * TAU,
        SciPAL::CudaComplex<double> * C,
        const size_t                * LDC,
        SciPAL::CudaComplex<double> * WORK,
        const size_t                * LWORK,
        int                         * INFO);
}

namespace SciPAL
{

struct lapack
{

    inline static void Init() {}
    inline static void Shutdown() {}
    // @sect4{Funktion: name}
    //!
    //! Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "LAPACK"; }

    //! Compile-time variable for identifying the side of PCIe bus this solver works on.
    static const ParallelArch arch = cpu;

    typedef typename archTraits<arch>::BlasType BW;

    //---------------------------------------------------------
    //SVD functions:
    template <typename T>
    static int SVD(SciPAL::Matrix<T, BW>&                                                   A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        //prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;
        int INFO = 0;
        char jobu[] = "S";
        char jobvt[] = "S";
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> RWORK(5*min_mn);        // Should be NULL if T is not a complex type.

        //get optimal workspace:
        gesvd(jobu, jobvt, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), &INFO);

        if(INFO == 0)
        {
            LWORK = size_t(abs(WORK(0)));
            WORK.reinit(LWORK);
        }
        else
        {
            check_status(INFO);
            return -1;
        }
        gesvd(jobu, jobvt, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), &INFO);
        return INFO;
    }
    template <typename T>
    static int SVD(const SciPAL::Matrix<T, BW>&                                             A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
        SciPAL::Matrix<T,BW> A_tmp(A);
        return SVD(A_tmp, U, S, Vt);
    }

    //---------------------------------------------------------
    //LUD functions:
    template <typename T>
    static int LUD(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<int, BW>&     P)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        P.reinit(min_mn);

        int INFO;
        getrf(&m, &n, A.data(), &(A.leading_dim), P.data(), &INFO);
        return INFO;
    }
    template <typename T>
    static int LUD(const SciPAL::Matrix<T,BW>&  A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return LUD(A_tmp, P, L ,U);
    }
    template <typename T>
    static int LUD(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m,n);

        L.reinit(m, min_mn);//   m,min_mn
        U.reinit(min_mn, n);//   min_mn,n

        int INFO = LUD(A, P);

        lower_upper_triangle(A, L, U);

        return INFO;
    }

    //---------------------------------------------------------
    //QRF functions:
    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>&           A,
                   SciPAL::Vector<T, BW>&           TAU)
    {
        const size_t m = A.n_rows();
        const size_t n = A.n_cols();
        const size_t min_mn = std::min(m, n);
        int INFO;
        //prepare call to LAPACK function:

        SciPAL::Vector<T, BW> WORK(1);
        TAU.reinit(min_mn);
        size_t LWORK = -1;
        //get optimal workspace:
        geqrf(&m, &n, A.data(), &(A.leading_dim), TAU.data(), WORK.data(), &LWORK, &INFO);
        if(INFO == 0)
        {
           LWORK = (size_t)abs(WORK(0));
           WORK.reinit(LWORK);
        }
        else
        {
          check_status(INFO);
          return -1;
        }
        geqrf(&m, &n, A.data(), &(A.leading_dim), TAU.data(), WORK.data(), &LWORK, &INFO);

        return INFO;
    }

    template <typename T>
    static int QRF(const SciPAL::Matrix<T, BW>&     A,
                   SciPAL::Matrix<T, BW>&           Output,
                   SciPAL::Vector<T, BW>&           TAU)
    {
        Output = A;
        return QRF(Output, TAU);
    }
    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>&           A,
                   SciPAL::Matrix<T, BW>&           Q,
                   SciPAL::Matrix<T, BW>&           R)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        SciPAL::Vector<T, BW> TAU(min_mn);

        //reinitialization
        Q.reinit(m, min_mn);
        R.reinit(min_mn, n);

        QRF(A, TAU);

        upper_triangle(A, R);

        for(size_t j = 0; j < min_mn; j++)
            Q(j, j, T(1));

        MQRF(A, TAU, Q);

        return 0;
    }
    template <typename T>
    static int QRF(const SciPAL::Matrix<T, BW>&     A,
                   SciPAL::Matrix<T, BW>&           Q,
                   SciPAL::Matrix<T, BW>&           R)
    {
        SciPAL::Matrix<T, BW> tmp;

        QRF(tmp, Q, R);

        return 0;
    }

    template <typename T>
    static int MQRF(const SciPAL::Matrix<T, BW>&    A,
                    const SciPAL::Vector<T, BW>&    TAU,
                    SciPAL::Matrix<T, BW>&          C)
    {
        const size_t m = C.n_rows_active(), n = C.n_cols_active();
        const size_t k = TAU.n_elements_active;
        int INFO;

        //prepare call to LAPACK function:
        char SIDE[] = "L";
        char TRANS[] = "N";

        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;

        //get optimal workspace:
        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &(A.leading_dim), TAU.data(), C.data(), &(C.leading_dim), WORK.data(), &LWORK, &INFO);

        if(INFO == 0)
        {
            LWORK = (size_t)abs(WORK(0));
            WORK.reinit(LWORK);
        }
        else
        {
            check_status(INFO);
            return -1;
        }

        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &(A.leading_dim), TAU.data(), C.data(), &(C.leading_dim), WORK.data(), &LWORK, &INFO);

        return INFO;
    }
    //---------------------------------------------------------
    // LDL functions:
    template <typename T>
    static int LDL(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<int, BW>&     /*P*/,
                   SciPAL::Matrix<T, BW>&       /*L*/,
                   SciPAL::Matrix<T, BW>&       /*D*/)
    {
        A.reinit(1, 1);
        // This method does nothing!
        return -1;
    }
    template <typename T>
    static int LDL(const SciPAL::Matrix<T, BW>& A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       D)
    {
        // This method does nothing!
        SciPAL::Matrix<T, BW> tmp(A);
        return LDL(tmp, P, L, D);
    }

private:
    inline static void check_status(const int & status)
    {
        //! The following is for tracking lapack usage
        //! std::cout << "Checking LAPACK status:" << std::endl;
#ifdef DEBUG
        std::string lapack_errors(" ");

        if (status != 0)
        {
            if(status == -1)
                lapack_errors += "Something went wrong.."; //TODO: expand error list

            if (lapack_errors == " ")
                lapack_errors = "unknown LAPACK error state";
#ifdef QT_NO_DEBUG
            AssertThrow(false, dealii::ExcMessage(lapack_errors.c_str() ) );
#else
            Assert(false, dealii::ExcMessage(lapack_errors.c_str() ) );
#endif
        }
#endif
    }

    //Helper functions (move them to the Matrix class?):
    template <typename T>
    static void lower_upper_triangle(const SciPAL::Matrix<T, BW>&               A,
                                     SciPAL::Matrix<T, BW>&                     L,
                                     SciPAL::Matrix<T, BW>&                     U)
    {
        // dimension of A
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m,n);
        // set diagonal entries of L to 1.0
        // copy strict lower triangular part of A into L
        // copy upper triangular part of A into U
        for(size_t i = 0; i < n; i++)
        {
            BW::copy(((i + 1 < m) ? i + 1 : m),
                     A.data() + i * A.leading_dim,
                     A.stride,
                     U.data() + i * U.leading_dim,
                     U.stride);
            if (n <= min_mn)
            {
                BW::copy(m - i + 1,
                         A.data() + i * A.leading_dim + i,
                         A.stride,
                         L.data() + i * U.leading_dim + i,
                         L.stride);
                L(i, i, T(1));
            }
        }
    }

    template <typename T>
    static void lower_upper_triangle_elementwise(const SciPAL::Matrix<T, BW>&   A,
                                                 SciPAL::Matrix<T, BW>&         L,
                                                 SciPAL::Matrix<T, BW>&         U)
    {
        // dimension of A
        const size_t m = A.n_rows(), n = A.n_cols();

        //get L and U from A, while the diagonal of the lower triangle matrix is 1
        for(size_t i = 0; i < m; i++)
        {
            for(size_t j = 0; j < n; j++)
            {
                if(i > j)
                    L(i, j, A(i,j));
                else
                    U(i, j, A(i,j));
                if(i == j)
                    L(i, i, (T)1);
            }
        }
    }

    template <typename T>
    static void upper_triangle(const SciPAL::Matrix<T, BW>&                     A,
                               SciPAL::Matrix<T, BW>&                           U)
    {
        // dimension of A
        const size_t m = A.n_rows(), n = A.n_cols();

        // copy upper triangular part of A into U.
        // Because we copy column-wise we need as much copy operation as we have columns.

        for(size_t i = 0; i < n; i++)
            BW::copy(((i + 1 < m) ? i + 1 : m),
                     A.data() + i * A.leading_dim,
                     A.stride,
                     U.data() + i * U.leading_dim,
                     U.stride);

    }

    template <typename T>
    static void lower_triangle(const SciPAL::Matrix<T, BW>&                     A,
                               SciPAL::Matrix<T, BW>&                           U)
    {
        // dimension of A
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        // copy lower triangular part of A into U (assume the diagonal elements are 1).
        // Because we copy column-wise we need as much copy operation as we have columns.

        for(size_t i = 0; i < min_mn; i++)
        {
            BW::copy(m - i + 1,
                     A.data() + i * A.leading_dim + i,
                     A.stride,
                     U.data() + i * U.leading_dim + i,
                     U.stride);
            U(i, i, T(1));
        }

    }

    template <typename T>
    static void upper_triangle_elementwise(const SciPAL::Matrix<T, BW>&         A,
                                           SciPAL::Matrix<T, BW>&               U)
    {
        // dimension of A
        const size_t m = A.n_rows(), n = A.n_cols();

        //get U from A:
        for(size_t i = 0; i < m; i++)
            for(size_t j = 0; j < n; j++)
                if(i <= j)
                    U(i, j, A(i, j));
    }


    //---------------------------------------------------------
    //SVD lapack-wrapper:
    inline static void gesvd(char * JOBU, char * JOBVT, const size_t * M, const size_t * N, float * A, const size_t * LDA, float * S, float * U, const size_t * LDU, float * VT, const size_t * LDVT, float * WORK, const size_t * LWORK, float * /*RWORK*/, int * INFO)
    {
        return sgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, const size_t * M, const size_t * N, double * A, const size_t * LDA, double * S, double * U, const size_t * LDU, double * VT, const size_t * LDVT, double * WORK, const size_t * LWORK, double * /*RWORK*/, int * INFO)
    {
        return dgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, const size_t * M, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, float * S, SciPAL::CudaComplex<float> * U, const size_t * LDU, SciPAL::CudaComplex<float> * VT, const size_t * LDVT, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, float * RWORK, int * INFO)
    {
        return cgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, const size_t * M, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, double * S, SciPAL::CudaComplex<double> * U, const size_t * LDU, SciPAL::CudaComplex<double> * VT, const size_t * LDVT, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, double * RWORK, int * INFO)
    {
        return zgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
    }

    //---------------------------------------------------------
    //LUD lapack-wrapper:
    inline static void getrf(const size_t * M, const size_t * N, float * A, const size_t * LDA, int * IPIV, int * INFO)
    {
        sgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(const size_t * M, const size_t * N, double * A, const size_t * LDA, int * IPIV, int * INFO)
    {
        dgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(const size_t * M, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, int * IPIV, int * INFO)
    {
        cgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(const size_t * M, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, int * IPIV, int * INFO)
    {
        zgetrf_(M, N, A, LDA, IPIV, INFO);
    }

    //---------------------------------------------------------
    //QRF lapack-wrapper:
    inline static void geqrf(const size_t * M, const size_t * N, float * A, const size_t * LDA, float * TAU, float * WORK, const size_t * LWORK, int * INFO)
    {
        sgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(const size_t * M, const size_t * N, double * A, const size_t * LDA, double * TAU, double * WORK, const size_t * LWORK, int * INFO)
    {
        dgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(const size_t * M, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, SciPAL::CudaComplex<float> * TAU, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, int * INFO)
    {
        cgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(const size_t * M, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, SciPAL::CudaComplex<double> * TAU, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, int * INFO)
    {
        zgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }

    inline static void xxmqr(char * SIDE, char * TRANS, const size_t * M, const size_t * N, const size_t * K, float * A, const size_t * LDA, float * TAU, float * C, const size_t * LDC, float * WORK, const size_t * LWORK, int * INFO)
    {
        sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, const size_t * M, const size_t * N, const size_t * K, double * A, const size_t * LDA, double * TAU, double * C, const size_t * LDC, double * WORK, const size_t * LWORK, int * INFO)
    {
        dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, const size_t * M, const size_t * N, const size_t * K, SciPAL::CudaComplex<float> * A, const size_t * LDA, SciPAL::CudaComplex<float> * TAU, SciPAL::CudaComplex<float> * C, const size_t * LDC, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, int * INFO)
    {
        cunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, const size_t * M, const size_t * N, const size_t * K, SciPAL::CudaComplex<double> * A, const size_t * LDA, SciPAL::CudaComplex<double> * TAU, SciPAL::CudaComplex<double> * C, const size_t * LDC, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, int * INFO)
    {
        zunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }

}; //END struct lapack
} //END namespace SciPAL
#endif // LAPACK_WRAPPER_HH
