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
        char    * JOBU,
        char    * JOBVT,
        int     * M,
        int     * N,
        float   * A,
        int     * LDA,
        float   * S,
        float   * U,
        int     * LDU,
        float   * VT,
        int     * LDVT,
        float   * WORK,
        int     * LWORK,
        int     * INFO);
void dgesvd_(
        char    * JOBU,
        char    * JOBVT,
        int     * M,
        int     * N,
        double  * A,
        int     * LDA,
        double  * S,
        double  * U,
        int     * LDU,
        double  * VT,
        int     * LDVT,
        double  * WORK,
        int     * LWORK,
        int     * INFO);
void cgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<float>  * A,
        int                         * LDA,
        float                       * S,
        SciPAL::CudaComplex<float>  * U,
        int                         * LDU,
        SciPAL::CudaComplex<float>  * VT,
        int                         * LDVT,
        SciPAL::CudaComplex<float>  * WORK,
        int                         * LWORK,
        float                       * RWORK,
        int                         * INFO);
void zgesvd_(
        char                        * JOBU,
        char                        * JOBVT,
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<double> * A,
        int                         * LDA,
        double                      * S,
        SciPAL::CudaComplex<double> * U,
        int                         * LDU,
        SciPAL::CudaComplex<double> * VT,
        int                         * LDVT,
        SciPAL::CudaComplex<double> * WORK,
        int                         * LWORK,
        double                      * RWORK,
        int                         * INFO);

//LUD routines:
void sgetrf_(
        int     * M,
        int     * N,
        float   * A,
        int     * LDA,
        int     * IPIV,
        int     * INFO);
void dgetrf_(
        int     * M,
        int     * N,
        double  * A,
        int     * LDA,
        int     * IPIV,
        int     * INFO);
void cgetrf_(
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<float>  * A,
        int                         * LDA,
        int                         * IPIV,
        int                         * INFO);
void zgetrf_(
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<double> * A,
        int                         * LDA,
        int                         * IPIV,
        int                         * INFO);

//QRF routines:
void sgeqrf_(
        int     * M,
        int     * N,
        float   * A,
        int     * LDA,
        float   * TAU,
        float   * WORK,
        int     * LWORK,
        int     * INFO);
void dgeqrf_(
        int     * M,
        int     * N,
        double  * A,
        int     * LDA,
        double  * TAU,
        double  * WORK,
        int     * LWORK,
        int     * INFO);
void cgeqrf_(
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<float>  * A,
        int                         * LDA,
        SciPAL::CudaComplex<float>  * TAU,
        SciPAL::CudaComplex<float>  * WORK,
        int                         * LWORK,
        int                         * INFO);
void zgeqrf_(
        int                         * M,
        int                         * N,
        SciPAL::CudaComplex<double> * A,
        int                         * LDA,
        SciPAL::CudaComplex<double> * TAU,
        SciPAL::CudaComplex<double> * WORK,
        int                         * LWORK,
        int                         * INFO);

//MQR routines:
void sormqr_(
        char    * SIDE,
        char    * TRANS,
        int     * M,
        int     * N,
        int     * K,
        float   * A,
        int     * LDA,
        float   * TAU,
        float   * C,
        int     * LDC,
        float   * WORK,
        int     * LWORK,
        int     * INFO);
void dormqr_(
        char    * SIDE,
        char    * TRANS,
        int     * M,
        int     * N,
        int     * K,
        double  * A,
        int     * LDA,
        double  * TAU,
        double  * C,
        int     * LDC,
        double  * WORK,
        int     * LWORK,
        int     * INFO);
void cunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        int                         * M,
        int                         * N,
        int                         * K,
        SciPAL::CudaComplex<float>  * A,
        int                         * LDA,
        SciPAL::CudaComplex<float>  * TAU,
        SciPAL::CudaComplex<float>  * C,
        int                         * LDC,
        SciPAL::CudaComplex<float>  * WORK,
        int                         * LWORK,
        int                         * INFO);
void zunmqr_(
        char                        * SIDE,
        char                        * TRANS,
        int                         * M,
        int                         * N,
        int                         * K,
        SciPAL::CudaComplex<double> * A,
        int                         * LDA,
        SciPAL::CudaComplex<double> * TAU,
        SciPAL::CudaComplex<double> * C,
        int                         * LDC,
        SciPAL::CudaComplex<double> * WORK,
        int                         * LWORK,
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
        int m = A.n_rows(), n = A.n_cols();
        int min_mn = std::min(m,n);

        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        //prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        int LWORK = -1;
        int INFO = 0;
        char jobu[] = "S";
        char jobvt[] = "S";
        SciPAL::Vector<typename SciPAL::VTraits<T, blas::arch>::NumberType, BW> RWORK(5*min_mn);        // Should be NULL if T is not a complex type.

        //get optimal workspace:
        gesvd(jobu, jobvt, &m, &n, A.data(), &m, \
                S.data(), U.data(), &m, \
                Vt.data(), &min_mn, WORK.data(),\
                &LWORK, RWORK.data(), &INFO);

        if(INFO == 0)
        {
            LWORK = (int)abs(WORK(0));
            WORK.reinit(LWORK);
        }
        else
        {
            check_status(INFO);
            return -1;
        }
        gesvd(jobu, jobvt, &m, &n, A.data(), &m, \
                S.data(), U.data(), &m, \
                Vt.data(), &min_mn, WORK.data(),\
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
        int m = A.n_rows(), n = A.n_cols();
        int min_mn = std::min(m,n);

        P.reinit(min_mn);

        int INFO;
        getrf(&m, &n, A.data(), &m, P.data(), &INFO);
        return INFO;
    }
    template <typename T>
    static int LUD(const SciPAL::Matrix<T,BW>&  A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        SciPAL::Matrix<T,blas> A_tmp(A);
        return LUD(A_tmp,P,L,U);
    }
    template <typename T>
    static int LUD(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        int m = A.n_rows(), n = A.n_cols();
        int min_mn = std::min(m,n);

        L.reinit(m, min_mn);//   m,min_mn
        U.reinit(min_mn, n);//   min_mn,n

        int INFO = LUD(A,P);

        lower_upper_triangle(A,L,U);

        return INFO;
    }

    //---------------------------------------------------------
    //QRF functions:
    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>&           A,
                   SciPAL::Vector<T, BW>&           TAU)
    {
        int m = A.n_rows();
        int n = A.n_cols();
        int min_mn = std::min(m, n);
        int INFO;
        //prepare call to LAPACK function:

        SciPAL::Vector<T, BW> WORK(1);
        TAU.reinit(min_mn);
        int LWORK = -1;
        //get optimal workspace:
        geqrf(&m, &n, A.data(), &m, TAU.data(), WORK.data(), &LWORK, &INFO);
        if(INFO == 0)
        {
           LWORK = (int)abs(WORK(0));
           WORK.reinit(LWORK);
        }
        else
        {
          check_status(INFO);
          return -1;
        }
        geqrf(&m, &n, A.data(), &m, TAU.data(), WORK.data(), &LWORK, &INFO);

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
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = (m<=n?m:n);

        SciPAL::Vector<T, blas> TAU(min_mn);

        //reinitialization
        Q.reinit(m, m);
        R.reinit(m, n);

        QRF(A, TAU);

        dealii::Timer testt;

        testt.start();

        upper_triangle(A, R);

        for(unsigned int j = 0; j < m; j++)
            Q(j, j, (T)1);

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
        int m = C.n_rows_active(), n = C.n_cols_active();
        int k = TAU.n_elements_active;
        int INFO;

        //prepare call to LAPACK function:
        char SIDE[] = "L";
        char TRANS[] = "N";

        SciPAL::Vector<T, BW> WORK(1);
        int LWORK = -1;

        //get optimal workspace:
        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &m, TAU.data(), C.data(), &m, WORK.data(), &LWORK, &INFO);

        if(INFO == 0)
        {
            LWORK = (int)abs(WORK(0));
            WORK.reinit(LWORK);
        }
        else
        {
            check_status(INFO);
            return -1;
        }

        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &m, TAU.data(), C.data(), &m, WORK.data(), &LWORK, &INFO);

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
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = std::min(m,n);
        // set diagonal entries of L to 1.0
        // copy strict lower triangular part of A into L
        // copy upper triangular part of A into U
        //   ip1 = i+1
        //   im = i*m
        //   impip1 = i*m+i+1
        for(unsigned int i=0, ip1=1, im=0; i<min_mn-1; i++, ip1++, im+=m)
        {
            unsigned int impip1 = im+ip1;
            L(i, i, 1.0);
            BW::copy(m-ip1, A.data()+impip1, 1, L.data()+impip1  , 1);
            BW::copy(  ip1, A.data()+im    , 1, U.data()+i*min_mn, 1);
        }
        L(min_mn-1, min_mn-1, 1.0);
        if(n >= m)
        {
            BW::copy(m*(n-m+1), A.data()+m*(m-1), 1, U.data()+m*(m-1), 1);
        }
        else
        {
            BW::copy(m-n, A.data()+m*(n-1)+n, 1, L.data()+m*(n-1)+n, 1);
            BW::copy(n, A.data()+m*(n-1), 1, U.data()+n*(n-1), 1);
        }
    }

    template <typename T>
    static void lower_upper_triangle_elementwise(const SciPAL::Matrix<T, BW>&   A,
                                                 SciPAL::Matrix<T, BW>&         L,
                                                 SciPAL::Matrix<T, BW>&         U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        //unsigned int min_mn = std::min(m,n);

        //get L and U from A:
        for(unsigned int i=0;i<m;i++)
        {
            for(unsigned int j=0;j<n;j++)
            {
                if(i>j) L(i,j, A(i,j));
                else    U(i,j, A(i,j));
                if(i==j) L(i,i, (T)1);
            }
        }
    }

    template <typename T>
    static void upper_triangle(const SciPAL::Matrix<T, BW>&               A,
                                     SciPAL::Matrix<T, BW>&               U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();
        unsigned int min_mn = std::min(m,n);
        // copy upper triangular part of A into U

        for(unsigned int i = 0; i < n; i++)
            BW::copy(((i+1<m)?i+1:m), A.data() + i * A.leading_dim, A.stride, U.data() + i * U.leading_dim, U.stride);

    }

    template <typename T>
    static void upper_triangle_elementwise(const SciPAL::Matrix<T, blas>&       A,
                                           SciPAL::Matrix<T, blas>&             U)
    {
        // dimension of A
        unsigned int m = A.n_rows(), n = A.n_cols();

        //get L and U from A:
        for(unsigned int i = 0; i < m; i++)
            for(unsigned int j = 0; j < n; j++)
                if(i <= j)
                    U(i, j, A(i, j));
    }


    //---------------------------------------------------------
    //SVD lapack-wrapper:
    inline static void gesvd(char * JOBU, char * JOBVT, int * M, int * N, float * A, int * LDA, float * S, float * U, int * LDU, float * VT, int * LDVT, float * WORK, int * LWORK, float * /*RWORK*/, int * INFO)
    {
        return sgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, int * M, int * N, double * A, int * LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK, double * /*RWORK*/, int * INFO)
    {
        return dgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, int * M, int * N, SciPAL::CudaComplex<float> * A, int * LDA, float * S, SciPAL::CudaComplex<float> * U, int * LDU, SciPAL::CudaComplex<float> * VT, int * LDVT, SciPAL::CudaComplex<float> * WORK, int * LWORK, float * RWORK, int * INFO)
    {
        return cgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
    }
    inline static void gesvd(char * JOBU, char * JOBVT, int * M, int * N, SciPAL::CudaComplex<double> * A, int * LDA, double * S, SciPAL::CudaComplex<double> * U, int * LDU, SciPAL::CudaComplex<double> * VT, int * LDVT, SciPAL::CudaComplex<double> * WORK, int * LWORK, double * RWORK, int * INFO)
    {
        return zgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
    }

    //---------------------------------------------------------
    //LUD lapack-wrapper:
    inline static void getrf(int * M, int * N, float * A, int * LDA, int * IPIV, int * INFO)
    {
        sgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(int * M, int * N, double * A, int * LDA, int * IPIV, int * INFO)
    {
        dgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(int * M, int * N, SciPAL::CudaComplex<float> * A, int * LDA, int * IPIV, int * INFO)
    {
        cgetrf_(M, N, A, LDA, IPIV, INFO);
    }
    inline static void getrf(int * M, int * N, SciPAL::CudaComplex<double> * A, int * LDA, int * IPIV, int * INFO)
    {
        zgetrf_(M, N, A, LDA, IPIV, INFO);
    }

    //---------------------------------------------------------
    //QRF lapack-wrapper:
    inline static void geqrf(int * M, int * N, float * A, int * LDA, float * TAU, float * WORK, int * LWORK, int * INFO)
    {
        sgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(int * M, int * N, double * A, int * LDA, double * TAU, double * WORK, int * LWORK, int * INFO)
    {
        dgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(int * M, int * N, SciPAL::CudaComplex<float> * A, int * LDA, SciPAL::CudaComplex<float> * TAU, SciPAL::CudaComplex<float> * WORK, int * LWORK, int * INFO)
    {
        cgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }
    inline static void geqrf(int * M, int * N, SciPAL::CudaComplex<double> * A, int * LDA, SciPAL::CudaComplex<double> * TAU, SciPAL::CudaComplex<double> * WORK, int * LWORK, int * INFO)
    {
        zgeqrf_(M, N, A, LDA, TAU, WORK, LWORK, INFO);
    }

    inline static void xxmqr(char * SIDE, char * TRANS, int * M, int * N, int * K, float * A, int * LDA, float * TAU, float * C, int * LDC, float * WORK, int * LWORK, int * INFO)
    {
        sormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, int * M, int * N, int * K, double * A, int * LDA, double * TAU, double * C, int * LDC, double * WORK, int * LWORK, int * INFO)
    {
        dormqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, int * M, int * N, int * K, SciPAL::CudaComplex<float> * A, int * LDA, SciPAL::CudaComplex<float> * TAU, SciPAL::CudaComplex<float> * C, int * LDC, SciPAL::CudaComplex<float> * WORK, int * LWORK, int * INFO)
    {
        cunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }
    inline static void xxmqr(char * SIDE, char * TRANS, int * M, int * N, int * K, SciPAL::CudaComplex<double> * A, int * LDA, SciPAL::CudaComplex<double> * TAU, SciPAL::CudaComplex<double> * C, int * LDC, SciPAL::CudaComplex<double> * WORK, int * LWORK, int * INFO)
    {
        zunmqr_(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO);
    }

}; //END struct lapack
} //END namespace SciPAL
#endif // LAPACK_WRAPPER_HH
