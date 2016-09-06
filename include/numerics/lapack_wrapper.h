//@sect3{File: lapack_wrapper.h}
// This header includes the lapack wrappers for SVD, LUD and QRF.

#ifndef LAPACK_WRAPPER_HH
#define LAPACK_WRAPPER_HH

//STL includes
#include <string>

//BASE includes
#include <lac/cublas_Vector.h>
#include <base/VTraits.h>
#include <base/ParallelArch.h>
#include <base/ArchTraits.h>

#include <lac/cublas_Vector.h>

#include <limits>

#define DEBUG_SVD

//Get externally defined LAPACK functions:
extern "C" {

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

//SVD divide and conquer routines:
void sgesdd_(
        char            * JOBZ,
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
        int             * IWORK,
        int             * INFO);
void dgesdd_(
        char            * JOBZ,
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
        int             * IWORK,
        int             * INFO);
void cgesdd_(
        char                        * JOBZ,
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
        int                         * IWORK,
        int                         * INFO);
void zgesdd_(
        char                        * JOBZ,
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
        int                         * IWORK,
        int                         * INFO);


//--------------------------------------------
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

//--------------------------------------------
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

//--------------------------------------------
//QP3 / RRQR routines:
void sgeqp3_(
        const size_t    * M,
        const size_t    * N,
        float           * A,
        const size_t    * LDA,
        int             * JPVT,
        float           * TAU,
        float           * WORK,
        const size_t    * LWORK,
        int             * INFO);
void dgeqp3_(
        const size_t    * M,
        const size_t    * N,
        double          * A,
        const size_t    * LDA,
        int             * JPVT,
        double          * TAU,
        double          * WORK,
        const size_t    * LWORK,
        int             * INFO);
void cgeqp3_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        int                         * JPVT,
        SciPAL::CudaComplex<float>  * TAU,
        SciPAL::CudaComplex<float>  * WORK,
        const size_t                * LWORK,
        float                       * RWORK,
        int                         * INFO);
void zgeqp3_(
        const size_t                * M,
        const size_t                * N,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        int                         * JPVT,
        SciPAL::CudaComplex<double> * TAU,
        SciPAL::CudaComplex<double> * WORK,
        const size_t                * LWORK,
        double                      * RWORK,
        int                         * INFO);

//--------------------------------------------
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

//--------------------------------------------
//SYEV/HEEV routines:
void ssyev_(
        char                        * JOBZ,
        char                        * UPLO,
        const size_t                * N,
        float                       * A,
        const size_t                * LDA,
        float                       * W,
        float                       * WORK,
        const size_t                * LWORK,
        int                         * INFO);
void dsyev_(
        char                        * JOBZ,
        char                        * UPLO,
        const size_t                * N,
        double                      * A,
        const size_t                * LDA,
        double                      * W,
        double                      * WORK,
        const size_t                * LWORK,
        int                         * INFO);
void cheev_(
        char                        * JOBZ,
        char                        * UPLO,
        const size_t                * N,
        SciPAL::CudaComplex<float>  * A,
        const size_t                * LDA,
        float                       * W,
        SciPAL::CudaComplex<float>  * WORK,
        const size_t                * LWORK,
        float                       * RWORK,
        int                         * INFO);
void zheev_(
        char                        * JOBZ,
        char                        * UPLO,
        const size_t                * N,
        SciPAL::CudaComplex<double> * A,
        const size_t                * LDA,
        double                      * W,
        SciPAL::CudaComplex<double> * WORK,
        const size_t                * LWORK,
        double                      * RWORK,
        int                         * INFO);

//--------------------------------------------
//STEV routines:
void sstev_(
        char                        * JOBZ,
        const size_t                * N,
        float                       * D,
        float                       * E,
        float                       * Z,
        const size_t                * LDZ,
        float                       * WORK,
        int                         * INFO);

void dstev_(
        char                        * JOBZ,
        const size_t                * N,
        double                      * D,
        double                      * E,
        double                      * Z,
        const size_t                * LDZ,
        double                      * WORK,
        int                         * INFO);

//--------------------------------------------
//LANGE norm routines:

float slange_(
        char                        * NORM,
        const size_t                * M,
        const size_t                * N,
        const float                 * A,
        const size_t                * LDA,
        float                       * WORK );


double dlange_(
        char                        * NORM,
        const size_t                * M,
        const size_t                * N,
        const double                * A,
        const size_t                * LDA,
        double                      * WORK );


float clange_(
        char                                * NORM,
        const size_t                        * M,
        const size_t                        * N,
        const SciPAL::CudaComplex<float>    * A,
        const size_t                        * LDA,
        float                               * WORK );

double zlange_(
        char                                * NORM,
        const size_t                        * M,
        const size_t                        * N,
        const SciPAL::CudaComplex<double>   * A,
        const size_t                        * LDA,
        double                              * WORK );

}

namespace SciPAL
{
// @sect4{struct lapack: LAPACK S,D,C,Z wrapper functions}
//
// We remove the type dependence of the functions provided by LAPACK by using C++'s polymorphism and templates.
// For the numerical purpose of the functions have a look at the LAPACK doc.
//
struct lapack
{
    // @sect5{Function: Init & Shutdown}
    // These functions aren't implemented since there is no need for a previous or common initialization for the used Lapack functions.
    inline static void Init() {}
    inline static void Shutdown() {}
    //---------------------------------------------------------
    // @sect5{Function: name (lapack)}
    //
    // Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "LAPACK"; }
    //---------------------------------------------------------

    // Compile-time variable for identifying the side of PCIe bus this solver works on:
    static const ParallelArch arch = cpu;

    typedef typename archTraits<arch>::BlasType BW;

    //---------------------------------------------------------
    // @sect5{SVD functions}
    // @sect6{Function: SVD}

    // Calculate singular-value-decomposition of $m \times n$ Matrix $A$ of type s, d, c, z with $A=U\cdot S\cdot Vt$.
    // \param A : input matrix, contents are overwritten during lapack function call.
    // \param U : output matrix of dimension $m \times min(m,n)$, contains left singular vectors.
    // \param S : output vector of length $min(m,n)$, contains singular values of real type.
    // \param Vt : output matrix of dimension $min(m,n) \times n$, contains right singular vectors.
    template <typename T>
    static int SVD(SciPAL::Matrix<T, BW>&                                                   A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD( Mat, Mat, Vec, Mat)\n";
#endif
#endif
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        //Resize output-matrices:
        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        //Prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;
        int INFO = 0;
	
        char jobu[] = "S";//The first min(m,n) columns/rows of U/V**T (the left/right singular vectors) are returned in the array U/Vt
        char jobvt[] = "S";//
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> RWORK(5*min_mn);        // Should be NULL if T is not a complex type.

        //Get optimal workspace:
        gesvd(jobu, jobvt, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), &INFO);

        //Resize workspace according to workspace query:
        if(INFO == 0)
        {
            LWORK = size_t(abs(WORK(0)));
            WORK.reinit(LWORK);
        }
        else //Throw ecxeption in case of error.
        {           
            check_status(INFO);
            return -1; //<-useless
        }

	//Call calculation routine:
        gesvd(jobu, jobvt, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), &INFO);
        return INFO;
    }

    // @sect6{Function: SVD (const A)}
    // Same as above, but keeps $A$ unchanged by creating a working copy $A_{tmp} = A$.
    template <typename T>
    static int SVD(const SciPAL::Matrix<T, BW>&                                             A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD( const Mat, Mat, Vec, Mat)\n";
#endif
#endif
        SciPAL::Matrix<T,BW> A_tmp(A);
        return SVD(A_tmp, U, S, Vt);
    }

    template <typename T>
    static int SVD(SciPAL::Matrix<T, BW>&                                                   A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD( Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;

        int result = SVD(A, U, S_tmp, Vt);

        // TODO: implement something like matix<diag> = vector
        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);
        for (unsigned int i = 0; i < S_tmp.n_elements_active; i++)
            S(i, i, S_tmp(i));

        return result;
    }

    template <typename T>
    static int SVD(const SciPAL::Matrix<T, BW>&                                             A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD( const Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        int result = SVD(A, U, S_tmp, Vt);
        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);
        for (unsigned int i = 0; i < S_tmp.n_elements_active; i++)
            S(i, i, S_tmp(i));
        return result;
    }

    template <typename T>
    static int SVD_truncate(const SciPAL::Matrix<T, BW>&                                    A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt,
                   typename SciPAL::VTraits<T, BW::arch>::NumberType MaxDiscWeight,
                   size_t MaxChi)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD_truncate( const Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        SciPAL::Matrix<T, BW> U_tmp,Vt_tmp;

        int result = SVD(A, U, S_tmp, Vt);

        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);


        // Calculate SV-Index at which either the demanded discarded weight is fullfilled
        // or the maximum dimension Chi is reached:
        typename SciPAL::VTraits<T, BW::arch>::NumberType DiscWeight = 0.0;
        unsigned int stop_index;

        for (unsigned int i = S_tmp.n_elements_active-1; i > 0 && DiscWeight < MaxDiscWeight && (S_tmp.n_elements_active-i)<=MaxChi ; i--)
        {
            DiscWeight += S_tmp(i)*S_tmp(i);
            stop_index = i;
        }

        S.reinit(stop_index, stop_index);
        for (unsigned int i = 0; i < stop_index; i++)
            S(i, i, S_tmp(i));
        SciPAL::SubMatrixView<T, BW> u_sub(U, 0, U.n_rows(), 0, stop_index);
        SciPAL::SubMatrixView<T, BW> vt_sub(Vt, 0, stop_index, 0, Vt.n_cols());
        U_tmp = u_sub;
        Vt_tmp = vt_sub;
        U = U_tmp;
        Vt = Vt_tmp;
        return result;
    }

    template <typename T>
    static int SVD_truncate(SciPAL::Matrix<T, BW>&                                           A,
                   SciPAL::Matrix<T, BW>&                                                    U,
                   SciPAL::Matrix<T, BW>&                                                    S,
                   SciPAL::Matrix<T, BW>&                                                   Vt,
                   typename SciPAL::VTraits<T, BW::arch>::NumberType             MaxDiscWeight,
                   size_t                                                               MaxChi)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SVD_truncate( Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        SciPAL::Matrix<T, BW> U_tmp,Vt_tmp;

        int result = SVD(A, U, S_tmp, Vt);

        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);


        // Calculate SV-Index at which either the demanded discarded weight is fullfilled
        // or the maximum dimension Chi is reached:
        typename SciPAL::VTraits<T, BW::arch>::NumberType DiscWeight = 0.0;
        size_t truncation_index = S_tmp.n_elements_active;

        if(MaxDiscWeight > 0)
        {
            while(DiscWeight < MaxDiscWeight)
            {
                truncation_index--;
                DiscWeight += S_tmp(truncation_index)*S_tmp(truncation_index);
            }
            truncation_index++;
        }

        truncation_index = std::min(truncation_index,MaxChi);

        S.reinit(truncation_index, truncation_index);
        for (unsigned int i = 0; i < truncation_index; i++)
            S(i, i, S_tmp(i));
        SciPAL::SubMatrixView<T, BW> u_sub(U, 0, U.n_rows(), 0, truncation_index);
        SciPAL::SubMatrixView<T, BW> vt_sub(Vt, 0, truncation_index, 0, Vt.n_cols());
        U_tmp = u_sub;
        Vt_tmp = vt_sub;
        U = U_tmp;
        Vt = Vt_tmp;
        return result;
    }


    template<typename T>
    static int RRQR_truncate(const SciPAL::Matrix<T, BW>&                       A,
                             SciPAL::Matrix<T, BW>&                             Q,
                             SciPAL::Matrix<T, BW>&                             R,
                             typename SciPAL::VTraits<T, BW::arch>::NumberType  MaxDiscWeight,
                             size_t                                             MaxChi
                             )
    {
        typedef typename SciPAL::VTraits<T, BW::arch>::NumberType NT;
#ifdef DEBUG
        SciPAL::Matrix<T, BW> A_tmp = A;
#endif
        SciPAL::Matrix<T, BW> R_tmp,Q_tmp;

        size_t min_mn = std::min(A.n_rows(),A.n_cols());

        // Allocate space for the results:
        SciPAL::Vector<int, BW> JPVT;

        // Execute the decomposition:
        int result = QP3(A, Q_tmp, R_tmp, JPVT);


        // Find largest sensible rank:

        NT R_max = abs(R_tmp(0,0));

        NT NT_precision = std::numeric_limits<NT>::epsilon();
        NT R_lower_limit = std::abs(R_max * NT_precision);

        //Find lowest non-zero diag-element
        size_t truncation_index=min_mn;

        while(abs(R_tmp(truncation_index-1,truncation_index-1))<R_lower_limit)
        {
            //std::cout << "R truncation_index:"<<truncation_index<<" \n";
            truncation_index--;

        };



#ifdef DEBUG

        std::cout << "R diag: \n";

        for(unsigned int j = 0; j < min_mn; j++)
            std::cout << R_tmp(j, j)<<"\t";
        std::cout << "\n";
        if(truncation_index < min_mn)
        {
            std::cout << "Lower limit chosen - truncation_index: "<<truncation_index<<", R(ti,ti) = "<<R_tmp(truncation_index,truncation_index)<<"\n";

            SciPAL::SubMatrixView<T,BW> R_sub(R_tmp,truncation_index,truncation_index);
            std::cout <<"Neglected ||R_2,2||: "<< R_sub.l2_norm()<<"\n";
        }
        else
        {
            std::cout <<"Matrix has full rank!\n";
        }
#endif

        NT sigma_cut = 0;
        NT ApproxDiscWeight = 0;
        if(MaxDiscWeight > 0)
        {
            // Is this neccesary since the rank-revealing truncation already stated, that R22 is empty?
            if(truncation_index < min_mn)
            {
                SciPAL::SubMatrixView<T,BW> R_sub(R_tmp,truncation_index,truncation_index);
                sigma_cut = R_sub.l2_norm();
                ApproxDiscWeight = std::pow(sigma_cut,2);
            }

            while(ApproxDiscWeight < MaxDiscWeight)
            {
                truncation_index--;
                SciPAL::RowVectorView<T,BW,const SciPAL::Matrix<T,BW>> R_row(R_tmp, truncation_index, truncation_index);
                sigma_cut += std::sqrt(std::pow(sigma_cut,2)+std::pow(R_row.l2_norm(),2));
                ApproxDiscWeight += std::pow(sigma_cut,2);

            }
            truncation_index++;
        }

        truncation_index = std::min(truncation_index,MaxChi);
#ifdef DEBUG
        if(truncation_index < min_mn)
        {
            std::cout << "Lower limit chosen - truncation_index: "<<truncation_index<<", R(ti,ti) = "<<R_tmp(truncation_index,truncation_index)<<"\n";

            SciPAL::SubMatrixView<T,BW> R_sub(R_tmp,truncation_index,truncation_index);
            std::cout <<"Neglected ||R_2,2||: "<< R_sub.l2_norm()<<"\n";
        }

        //SVD comparison:
        SciPAL::Vector<NT, BW> S_tmp;
        SciPAL::Matrix<T, BW> U_tmp,Vt_tmp;

        SVD(A_tmp, U_tmp, S_tmp, Vt_tmp);

        std::cout << "SV:\n";

        for (unsigned int i = 0; i <S_tmp.n_elements_active ; i++)
        {
            std::cout << S_tmp(i)<<"\t";
        }
        std::cout << "\n";
        std::cout << "DiscWeight:\n";
        S_tmp(S_tmp.n_elements_active-1, std::pow(S_tmp(S_tmp.n_elements_active-1),2));

        for (unsigned int i = S_tmp.n_elements_active-2; i>0  ; i--)
        {
            S_tmp(i,std::pow(S_tmp(i),2)+S_tmp(i+1));
        }
        S_tmp(0,std::pow(S_tmp(0),2)+S_tmp(1));

        for (unsigned int i = 0; i <S_tmp.n_elements_active ; i++)
        {
            std::cout << S_tmp(i)<<"\t";
        }
        std::cout << "\n";;
#endif

        // truncation

        SciPAL::SubMatrixView<T,BW> Q_relevant(Q_tmp, 0,Q_tmp.n_rows()  , 0 , truncation_index);
        SciPAL::SubMatrixView<T,BW> R_relevant(R_tmp, 0,truncation_index, 0 , R_tmp.n_cols());

        Q = Q_relevant;
        R = R_relevant;

        // Undo permutation on R
        R.permute_inv_col(JPVT);
        // Now, A is recomposed.
        return result;
    }
    //____________change to sdd begin


    // @sect6{Function: SDD}

    // Calculate singular-value-decomposition of $m \times n$ Matrix $A$ of type s, d, c, z with $A=U\cdot S\cdot Vt$ using lapacks gesdd.
    // \param A : input matrix, contents are overwritten during lapack function call.
    // \param U : output matrix of dimension $m \times min(m,n)$, contains left singular vectors.
    // \param S : output vector of length $min(m,n)$, contains singular values of real type.
    // \param Vt : output matrix of dimension $min(m,n) \times n$, contains right singular vectors.
    template <typename T, typename T2=int>
    static int SDD(SciPAL::Matrix<T, BW>&                                                   A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD( Mat, Mat, Vec, Mat)\n";
#endif
#endif
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);
        const size_t max_mn = std::max(m, n);

        //Resize output-matrices:
        U.reinit(m, min_mn);
        S.reinit(min_mn);
        Vt.reinit(min_mn, n);

        //Prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;
        SciPAL::Vector<T2, BW> IWORK(8*min_mn);

        int INFO = 0;

        char jobz[] = "S";//The first min(m,n) columns/rows of U/V**T (the left/right singular vectors) are returned in the array U/Vt
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> RWORK(std::max(5*min_mn*min_mn+5*min_mn, 2*max_mn*min_mn+2*min_mn*min_mn+min_mn));        // Should be NULL if T is not a complex type.

        //Get optimal workspace:
        gesdd(jobz, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), IWORK.data(), &INFO);

        //Resize workspace according to workspace query:
        if(INFO == 0)
        {
            LWORK = size_t(abs(WORK(0)));
            WORK.reinit(LWORK);
        }
        else //Throw ecxeption in case of error.
        {
            check_status(INFO);
            return -1; //<-useless
        }

    //Call calculation routine:
        gesdd(jobz, &m, &n, A.data(), &(A.leading_dim), \
                S.data(), U.data(), &(U.leading_dim), \
                Vt.data(), &(Vt.leading_dim), WORK.data(),\
                &LWORK, RWORK.data(), IWORK.data(), &INFO);
        return INFO;
    }

    // @sect6{Function: SDD (const A)}
    // Same as above, but keeps $A$ unchanged by creating a working copy $A_{tmp} = A$.
    template <typename T>
    static int SDD(const SciPAL::Matrix<T, BW>&                                             A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD( const Mat, Mat, Vec, Mat)\n";
#endif
#endif
        SciPAL::Matrix<T,BW> A_tmp(A);
        return SDD(A_tmp, U, S, Vt);
    }

    template <typename T>
    static int SDD(SciPAL::Matrix<T, BW>&                                                   A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD( Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;

        int result = SDD(A, U, S_tmp, Vt);

        // TODO: implement something like matix<diag> = vector
        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);
        for (unsigned int i = 0; i < S_tmp.n_elements_active; i++)
            S(i, i, S_tmp(i));

        return result;
    }

    template <typename T>
    static int SDD(const SciPAL::Matrix<T, BW>&                                             A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD( const Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        int result = SDD(A, U, S_tmp, Vt);
        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);
        for (unsigned int i = 0; i < S_tmp.n_elements_active; i++)
            S(i, i, S_tmp(i));
        return result;
    }

    template <typename T>
    static int SDD_truncate(const SciPAL::Matrix<T, BW>&                                    A,
                   SciPAL::Matrix<T, BW>&                                                   U,
                   SciPAL::Matrix<T, BW>&                                                   S,
                   SciPAL::Matrix<T, BW>&                                                   Vt,
                   typename SciPAL::VTraits<T, BW::arch>::NumberType MaxDiscWeight,
                   size_t MaxChi)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD_truncate( const Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        SciPAL::Matrix<T, BW> U_tmp,Vt_tmp;

        int result = SDD(A, U, S_tmp, Vt);

        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);


        // Calculate SV-Index at which either the demanded discarded weight is fullfilled
        // or the maximum dimension Chi is reached:
        typename SciPAL::VTraits<T, BW::arch>::NumberType DiscWeight = 0.0;
        unsigned int stop_index;

        for (unsigned int i = S_tmp.n_elements_active-1; i > 0 && DiscWeight < MaxDiscWeight && (S_tmp.n_elements_active-i)<=MaxChi ; i--)
        {
            DiscWeight += S_tmp(i)*S_tmp(i);
            stop_index = i;
        }

        S.reinit(stop_index, stop_index);
        for (unsigned int i = 0; i < stop_index; i++)
            S(i, i, S_tmp(i));
        SciPAL::SubMatrixView<T, BW> u_sub(U, 0, U.n_rows(), 0, stop_index);
        SciPAL::SubMatrixView<T, BW> vt_sub(Vt, 0, stop_index, 0, Vt.n_cols());
        U_tmp = u_sub;
        Vt_tmp = vt_sub;
        U = U_tmp;
        Vt = Vt_tmp;
        return result;
    }

    template <typename T>
    static int SDD_truncate(SciPAL::Matrix<T, BW>&                                           A,
                   SciPAL::Matrix<T, BW>&                                                    U,
                   SciPAL::Matrix<T, BW>&                                                    S,
                   SciPAL::Matrix<T, BW>&                                                   Vt,
                   typename SciPAL::VTraits<T, BW::arch>::NumberType             MaxDiscWeight,
                   size_t                                                               MaxChi)
    {
#ifdef DEBUG
#ifdef DEBUG_SVD
        std::cerr <<"Using SDD_truncate( Mat, Mat, Mat, Mat)\n";
#endif
#endif
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> S_tmp;
        SciPAL::Matrix<T, BW> U_tmp,Vt_tmp;

        int result = SDD(A, U, S_tmp, Vt);

        S.reinit(S_tmp.n_elements_active, S_tmp.n_elements_active);


        // Calculate SV-Index at which either the demanded discarded weight is fullfilled
        // or the maximum dimension Chi is reached:
        typename SciPAL::VTraits<T, BW::arch>::NumberType DiscWeight = 0.0;
        size_t truncation_index = S_tmp.n_elements_active;

        if(MaxDiscWeight > 0)
        {
            while(DiscWeight < MaxDiscWeight)
            {
                truncation_index--;
                DiscWeight += S_tmp(truncation_index)*S_tmp(truncation_index);
            }
            truncation_index++;
        }

        truncation_index = std::min(truncation_index,MaxChi);

        S.reinit(truncation_index, truncation_index);
        for (unsigned int i = 0; i < truncation_index; i++)
            S(i, i, S_tmp(i));
        SciPAL::SubMatrixView<T, BW> u_sub(U, 0, U.n_rows(), 0, truncation_index);
        SciPAL::SubMatrixView<T, BW> vt_sub(Vt, 0, truncation_index, 0, Vt.n_cols());
        U_tmp = u_sub;
        Vt_tmp = vt_sub;
        U = U_tmp;
        Vt = Vt_tmp;
        return result;
    }

    //____________change to sdd end

    //---------------------------------------------------------
    // @sect5{LUD functions}

    // @sect6{Function: LUD}
    // Calculate LU-decomposition of $m \times n$ Matrix $A$ of type s, d, c, z with $A=P\cdot L\cdot U$.
    //
    // $P$ is a permutation matrix.
    // $L$ is lower triangular (lower trapezoidal if $m > n$) with unit diagonal elements.
    // $U$ is upper triangular (upper trapezoidal if $m < n$).
    // \param A : input/output matrix, contains $L$ (except for its unit diagonal) and $U$ after function call.
    // \param P : output vector, describes permutation matrix as pivot indices; for $1 \leq i \leq min(m,n)$, row $i$ of the matrix was interchanged with row $P(i)$.
    template <typename T, typename T2=int>
    static int LUD(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<T2, BW>&     P)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

	//Resize output:
        P.reinit(min_mn);

        int INFO;

	//Call LAPACK routine:
        getrf(&m, &n, A.data(), &(A.leading_dim), P.data(), &INFO);
        return INFO;
    }

    // @sect6{Function: LUD (L,U seperate)}
    // Same as above but explicitly copies the content of $A$ into matrices $L$ and $U$.
    // \param A : input/output matrix, contains $L$ (except for its unit diagonal) and $U$ after function call.
    // \param P : output vector, describes permutation matrix as pivot indices; for $1 \leq i \leq min(m,n)$, row $i$ of the matrix was interchanged with row $P(i)$.
    // \param L : output matrix; lower triangular (lower trapezoidal if $m > n$) of dimension $m \times min(m,n)$ with unit diagonal elements.
    // \param U : output matrix; upper triangular (upper trapezoidal if $m < n$) of dimension $min(m,n) \times n$.
    template <typename T, typename T2=int>
    static int LUD(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<T2, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m,n);

	//Resize output:
        L.reinit(m, min_mn);
        U.reinit(min_mn, n);

	//Call LU-decomposition:
        int INFO = LUD(A, P);

	//Copy contents of A into L and U:
        lower_upper_triangle(A, L, U);

        return INFO;
    }

    // @sect6{Function: LUD (const A, L/U seperate)}
    // Same as above, but keeps $A$ unchanged by creating a working copy $A_{tmp} = A$.
    template <typename T, typename T2=int>
    static int LUD(const SciPAL::Matrix<T,BW>&  A,
                   SciPAL::Vector<T2, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       U)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        return LUD(A_tmp, P, L ,U);
    }



     //---------------------------------------------------------
   // @sect5{QRF functions}

    // @sect6{Function: QRF}
    // Calculate QR-factorization of $m \times n$ Matrix $A$ of type s, d, c, z with $A=Q\cdot R$.

    // $Q$ orthogonal matrix of dimension $m \times m$.
    // $R$ is upper triangular matrix of dimension $m \times n$ (upper trapezoidal if $m < n$).

    // Lapack doc: 'The matrix $Q$ is represented as a product of elementary reflectors $Q = H(1) H(2) ... H(k)$ , where $k = min(m,n)$.
    // Each $H(i)$ has the form $H(i) = I - TAU \cdot v \cdot v^T$ where $TAU$ is a real scalar, and $v$ is a real vector with $v(1:i-1) = 0$ and $v(i) = 1; v(i+1:m)$ is stored on exit in $A(i+1:m,i)$, and $TAU$ in $TAU(i)$.'

    // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
    // \param TAU : output vector of length $min(m,n)$ contains the scalar factors of the elementary reflector.

    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>&           A,
                   SciPAL::Vector<T, BW>&           TAU)
    {
        const size_t m = A.n_rows();
        const size_t n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        //Prepare call to LAPACK function:
        int INFO;
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;

        //Resize output:
        TAU.reinit(min_mn);

        //Get optimal workspace:
        geqrf(&m, &n, A.data(), &(A.leading_dim), TAU.data(), WORK.data(), &LWORK, &INFO);

        //Resize workspace according to workspace query:
        if(INFO == 0)
        {
            LWORK = (size_t)abs(WORK(0));
            WORK.reinit(LWORK);
        }
        else//Throw ecxeption in case of error.
        {
            check_status(INFO);
            return -1; //<-useless
        }
        //Call calculation routine:
        geqrf(&m, &n, A.data(), &(A.leading_dim), TAU.data(), WORK.data(), &LWORK, &INFO);

        return INFO;
    }

    // @sect6{Function: QRF (const A)}
    // Same as above, but keeps $A$ unchanged by using a working copy $A_{tmp} = A$.
    template <typename T>
    static int QRF(const SciPAL::Matrix<T, BW>&     A,
                   SciPAL::Matrix<T, BW>&           A_tmp,
                   SciPAL::Vector<T, BW>&           TAU)
    {
        A_tmp = A;
        return QRF(A_tmp, TAU);
    }

    // @sect6{Function: QRF (Q,R seperate)}
    // Same as above but explicitly processec the output content of $A$ into matrices $Q$ and $R$.
    // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
    // \param Q : output vector, orthogonal matrix of dimension $m \times min(m,n)$.
    // \param R : output matrix; upper triangular (upper trapezoidal if $m < n$) of dimension $min(m,n) \times n$.
    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>&           A,
                   SciPAL::Matrix<T, BW>&           Q,
                   SciPAL::Matrix<T, BW>&           R)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        SciPAL::Vector<T, BW> TAU(min_mn);

        //Resize output:
        Q.reinit(m, min_mn);
        R.reinit(min_mn, n);

        //Call QR-factorizaion:
        QRF(A, TAU);

        //Copy upper triangle/trapezoidal into R:
        upper_triangle(A, R);

        //Calculate Q by applying the product of elementary reflectors:
        //Make input matrix Q equal to trapezoidal unit matrix
        for(size_t j = 0; j < min_mn; j++)
            Q(j, j, T(1));
        //and call Lapack function to evaluate the product:
        MQRF(A, TAU, Q);

        return 0;
    }

    // @sect6{Function: QRF (const A, Q/R seperate)}
    // Same as above, but keeps $A$ unchanged by using a working copy $A_{tmp} = A$
    template <typename T>
    static int QRF(const SciPAL::Matrix<T, BW>&     A,
                   SciPAL::Matrix<T, BW>&           Q,
                   SciPAL::Matrix<T, BW>&           R)
    {
        SciPAL::Matrix<T, BW> A_tmp(A);
        QRF(A_tmp, Q, R);
        return 0;
    }

    // @sect6{Function: QRF (Q, R)}
    // Same as above, but using Q as input. Might be worth a look for optimization!
    template <typename T>
    static int QRF(SciPAL::Matrix<T, BW>& Q,
                   SciPAL::Matrix<T, BW>& R)
    {
        SciPAL::Matrix<T, BW> A;
        A = Q;

        return QRF(A, Q, R);
    }

    //---------------------------------------------------------
  // @sect5{QP3 functions}

   // @sect6{Function: QP3}
   // Calculate QR-factorization of $m \times n$ Matrix $A$ of type s, d, c, z with $A\cdot P=Q\cdot R$.

   // $P$ permutation matrix of dimensions $n \times m$.

   // $Q$ orthogonal matrix of dimension $m \times m$.
   // $R$ is upper triangular matrix of dimension $m \times n$ (upper trapezoidal if $m < n$).

   // Lapack doc: 'The matrix $Q$ is represented as a product of elementary reflectors $Q = H(1) H(2) ... H(k)$ , where $k = min(m,n)$.
   // Each $H(i)$ has the form $H(i) = I - TAU \cdot v \cdot v^T$ where $TAU$ is a real scalar, and $v$ is a real vector with $v(1:i-1) = 0$ and $v(i) = 1; v(i+1:m)$ is stored on exit in $A(i+1:m,i)$, and $TAU$ in $TAU(i)$.'

   // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
   // \param TAU : output vector of length $min(m,n)$ contains the scalar factors of the elementary reflector.

   template <typename T, typename T2=int>
   static int QP3(SciPAL::Matrix<T, BW>&           A,
                  SciPAL::Vector<T, BW>&           TAU,
                  SciPAL::Vector<T2,BW>&       JPVT)
   {
       const size_t m = A.n_rows();
       const size_t n = A.n_cols();
       const size_t min_mn = std::min(m, n);

       //Prepare call to LAPACK function:
       int INFO;
       SciPAL::Vector<T, BW> WORK(1);
       size_t LWORK = -1;

       //Resize output:
       TAU.reinit(min_mn);

       //Resize output:
       JPVT.reinit(n);

       //Get optimal workspace:
       geqp3(&m, &n, A.data(), &(A.leading_dim),JPVT.data(), TAU.data(), WORK.data(), &LWORK, &INFO);

       //Resize workspace according to workspace query:
       if(INFO == 0)
       {
           LWORK = (size_t)abs(WORK(0));
           WORK.reinit(LWORK);
       }
       else//Throw ecxeption in case of error.
       {
           check_status(INFO);
           return -1; //<-useless
       }
       //Call calculation routine:
       geqp3(&m, &n, A.data(), &(A.leading_dim),JPVT.data() , TAU.data(), WORK.data(), &LWORK, &INFO);

       return INFO;
   }

   // @sect6{Function: QP3 (const A)}
   // Same as above, but keeps $A$ unchanged by using a working copy $A_{tmp} = A$.
   template <typename T, typename T2=int>
   static int QP3(const SciPAL::Matrix<T, BW>&     A,
                  SciPAL::Matrix<T, BW>&           A_tmp,
                  SciPAL::Vector<T, BW>&           TAU,
                  SciPAL::Vector<T2,BW>&       JPVT)
   {
       A_tmp = A;
       return QP3(A_tmp, TAU, JPVT);
   }

   // @sect6{Function: QP3 (Q,R seperate)}
   // Same as above but explicitly processec the output content of $A$ into matrices $Q$ and $R$.
   // \param A : input/output matrix, contains vectors $v$ in its strict lower triangular part and $R$ after function call.
   // \param Q : output vector, orthogonal matrix of dimension $m \times min(m,n)$.
   // \param R : output matrix; upper triangular (upper trapezoidal if $m < n$) of dimension $min(m,n) \times n$.
   template <typename T, typename T2=int>
   static int QP3(SciPAL::Matrix<T, BW>&           A,
                  SciPAL::Matrix<T, BW>&           Q,
                  SciPAL::Matrix<T, BW>&           R,
                  SciPAL::Vector<T2,BW>&       JPVT)
   {
       const size_t m = A.n_rows(), n = A.n_cols();
       const size_t min_mn = std::min(m, n);

       SciPAL::Vector<T, BW> TAU(min_mn);
       //Resize output:
       Q.reinit(m, min_mn);
       R.reinit(min_mn, n);
       JPVT.reinit(n);

       //Call QR-factorizaion:
       QP3(A, TAU, JPVT);

       //Copy upper triangle/trapezoidal into R:
       upper_triangle(A, R);

       //Calculate Q by applying the product of elementary reflectors:
       //Make input matrix Q equal to trapezoidal unit matrix
       for(size_t j = 0; j < min_mn; j++)
           Q(j, j, T(1));
       //and call Lapack function to evaluate the product:
       MQRF(A, TAU, Q);

       return 0;
   }

   // @sect6{Function: QP3 (const A, Q/R seperate)}
   // Same as above, but keeps $A$ unchanged by using a working copy $A_{tmp} = A$
   template <typename T, typename T2=int>
   static int QP3(const SciPAL::Matrix<T, BW>&     A,
                  SciPAL::Matrix<T, BW>&           Q,
                  SciPAL::Matrix<T, BW>&           R,
                  SciPAL::Vector<T2,BW>&       JPVT)
   {
       SciPAL::Matrix<T, BW> A_tmp(A);
       QP3(A_tmp, Q, R, JPVT);
       return 0;
   }

   // @sect6{Function: QP3 (Q, R)}
   // Same as above, but using Q as input. Might be worth a look for optimization!
   template <typename T, typename T2=int>
   static int QP3(SciPAL::Matrix<T, BW>&        Q,
                  SciPAL::Matrix<T, BW>&        R,
                  SciPAL::Vector<T2,BW>&    JPVT)
   {
       SciPAL::Matrix<T, BW> A;
       A = Q;

       return QP3(A, Q, R, JPVT);
   }

  // @sect6{Function: MQR}
  // Overwrites the matrix $C$ of dimension  with $Q\cdot C$ where $Q$ is a real orthogonal matrix defined as the product of $k$ elementary reflectors $Q = H(1) H(2) ... H(k)$ as returned by xGEQRF.
  // \param A : input matrix, contains results of xgeqrf function call / values of elementary reflectors.
  // \param TAU : input vector, contains results of xgeqrf function call.
  // \param C : output matrix; contains results of product $Q\cdot C$.
    template <typename T>
    static int MQRF(const SciPAL::Matrix<T, BW>&    A,
                    const SciPAL::Vector<T, BW>&    TAU,
                    SciPAL::Matrix<T, BW>&          C)
    {
        const size_t m = C.n_rows_active(), n = C.n_cols_active();
        const size_t k = TAU.n_elements_active;
        int INFO;

        // Prepare call to LAPACK function:
        char SIDE[] = "L"; // Apply Q from the left.
        char TRANS[] = "N";//
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;

        //Get optimal workspace:
        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &(A.leading_dim), TAU.data(), C.data(), &(C.leading_dim), WORK.data(), &LWORK, &INFO);

	//Resize workspace according to workspace query:
        if(INFO == 0)
        {
           LWORK = (size_t)abs(WORK(0));
           WORK.reinit(LWORK);
        }
        else //Throw ecxeption in case of error.
        {
            check_status(INFO);
            return -1; //<-useless
        }
	//Call calculation routine:
        xxmqr(SIDE, TRANS, &m, &n, &k, A.data(), &(A.leading_dim), TAU.data(), C.data(), &(C.leading_dim), WORK.data(), &LWORK, &INFO);

        return INFO;
    }


    //---------------------------------------------------------
   // @sect5{LDL functions}
   // @sect6{Function: LDL}
   // Not yet implemented.
    template <typename T>
    static int LDL(SciPAL::Matrix<T, BW>&       A,
                   SciPAL::Vector<int, BW>&     /*P*/,
                   SciPAL::Matrix<T, BW>&       /*L*/,
                   SciPAL::Matrix<T, BW>&       /*D*/)
    {
        // This method does nothing!
        A.reinit(1, 1);
        return -1;
    }

   // @sect6{Function: LDL (const A)}
   // Same as above, but keeps $A$ unchanged by using a working copy $A_{tmp} = A$.
   // Not yet implemented.
    template <typename T>
    static int LDL(const SciPAL::Matrix<T, BW>& A,
                   SciPAL::Vector<int, BW>&     P,
                   SciPAL::Matrix<T, BW>&       L,
                   SciPAL::Matrix<T, BW>&       D)
    {
        // This method does nothing!
        SciPAL::Matrix<T, BW> A_tmp(A);
        return LDL(A_tmp, P, L, D);
    }

    //---------------------------------------------------------
    // @sect5{Eigen functions}
    // @sect6{Function: Eigen}
    // \param A : Input/output matrix; on entry: symmetric/hermitsh matrix. On exit: contains the eigenvectors of A.
    // \param EV : Contains the real eigenvalues in ascending order
    template <typename T>
    static int Eigen(SciPAL::Matrix<T, BW>&                                                     A,
                     SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     EV)
    {
        Assert(A.n_rows() == A.n_cols(), dealii::ExcMessage("Dimension missmatch."));

        char jobz[] = "V";

        char uplo[] = "U";

        const size_t n = A.n_rows();

        EV.reinit(n);

//        std::cout << "INPUT MATRIX (before xxev):" << std::endl;
//        A.print();
//        std::cout << "CHECK THE EIGENVALUES (before xxev)!!!" << std::endl;
//        EV.print();

        //Prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;
        int INFO = 0;

        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> RWORK(3*n-2);        // Should be NULL if T is not a complex type.

        //Get optimal workspace:
        xxev(jobz, uplo, &n, A.data(), &(A.leading_dim), \
             EV.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);


        //Resize workspace according to workspace query:
        if(INFO == 0)
        {
            LWORK = size_t(abs(WORK(0)));
            WORK.reinit(LWORK);
        }
        else //Throw ecxeption in case of error.
        {
            check_status(INFO);
        }

        //Call calculation routine:
        xxev(jobz, uplo, &n, A.data(), &(A.leading_dim), \
             EV.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);
//        std::cout << "INPUT MATRIX (after xxev):" << std::endl;
//        A.print();
//        std::cout << "CHECK THE EIGENVALUES (after xxev)!!!" << std::endl;
//        EV.print();
        return INFO;
    }

    // @sect6{Function: Eigen}
    // \param A : Input matrix: symmetric/hermitsh matrix. As it is const, no eigenvectors are calculated.
    // \param EV : Contains the real eigenvalues in ascending order
    template <typename T>
    static int Eigen(const SciPAL::Matrix<T, BW>&                                               A,
                     SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     EV)
    {
        Assert(A.n_rows() == A.n_cols(), dealii::ExcMessage("Dimension missmatch."));

        SciPAL::Matrix<T, BW> A_tmp(A);

        char jobz[] = "N";

        char uplo[] = "U";

        const size_t n = A.n_rows();

        //Prepare call to LAPACK function:
        SciPAL::Vector<T, BW> WORK(1);
        size_t LWORK = -1;
        int INFO = 0;

        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> RWORK(3*n-2);        // Should be NULL if T is not a complex type.

        //Get optimal workspace:
        xxev(jobz, uplo, &n, A_tmp.data(), &(A_tmp.leading_dim), \
             EV.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);


        //Resize workspace according to workspace query:
        if(INFO == 0)
        {
            LWORK = size_t(abs(WORK(0)));
            WORK.reinit(LWORK);
        }
        else //Throw ecxeption in case of error.
        {
            check_status(INFO);
        }

        //Call calculation routine:
        xxev(jobz, uplo, &n, A_tmp.data(), &(A_tmp.leading_dim), \
             EV.data(), WORK.data(), &LWORK, RWORK.data(), &INFO);
        return INFO;
    }

    // @sect6{Function: Eigen}
    // \param D : Input/Output vector: n diagonal elements of the tridiagonal matrix A. On exit: Eigenvalues in ascending order.
    // \param E : Input/Output vector: (n-1) sub-diagonal elements of the tridiagonal matrix A, stored in (1, N-1). On Exit: Content get destroyed.
    // \param Z : Output matrix: Eigenvectors of the matrix A
    template <typename T>
    static int Eigen(SciPAL::Vector<T, BW>&     D,
                     SciPAL::Vector<T, BW>&     E,
                     SciPAL::Matrix<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     Z)
//    static int Eigen(SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     D,
//                     SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     E,
//                     SciPAL::Matrix<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     Z)
    {
//        if (D.n_elements_active > 1)
//            Assert(D.n_elements_active == (E.n_elements_active-1), dealii::ExcMessage("Vectors have wrong dimension."));

        char jobz[] = "V";
        const size_t n = D.n_elements_active;
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> WORK(n==1?1:(2*n)-2);
        int INFO = 0;

        Z.reinit(n, n);

        stev(jobz, &n, D.data(), E.data(), Z.data(), &(Z.leading_dim), WORK.data(), &INFO);

        return INFO;
    }

    // @sect6{Function: Eigen}
    // \param D : Input/Output vector: n diagonal elements of the tridiagonal matrix A. On exit: Eigenvalues in ascending order.
    // \param E : Input/Output vector: (n-1) sub-diagonal elements of the tridiagonal matrix A, stored in (1, N-1). On Exit: Content get destroyed.
    template <typename T>
    static int Eigen(SciPAL::Vector<T, BW>&     D,
                     SciPAL::Vector<T, BW>&     E)
//    static int Eigen(SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     D,
//                     SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>&     E)
    {
        Assert(D.n_elements_active == E.n_elements_active-1, dealii::ExcMessage("Vectors have wrong dimension."))

        char jobz[] = "N";
        const size_t n = D.n_elements_active;
        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> WORK(1);
        int INFO = 0;

        SciPAL::Matrix<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> Z(1, 1);

        stev(jobz, &n, D.data(), E.data(), Z.data(), &(Z.leading_dim), WORK.data(), &INFO);

        return INFO;
    }

    //---------------------------------------------------------

    // @sect6{Function: norm2}
    // \param M : Input matrix
    // Calculates the 2-norm of submatrix M
    template <typename T,typename BW>
    static typename SciPAL::VTraits<T, BW::arch>::NumberType norm2(const SciPAL::SubMatrixView<T,BW>& M)
    {
        char norm[] = "F"; //Frobenius aka 2-norm

        SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW> WORK(1);

        const size_t m = M.n_rows();
        const size_t n = M.n_cols();

        const size_t ldm = M.leading_dim();

        return lange(norm, &m, &n, M.data(), &ldm, WORK.data());
    }

private:

   // @sect5{Function: check_status}
   // Check if status parameter is equal to zero; otherwise throw exception.
   // Detailed error information is not included yet since for almost every lapack function, the given info parameter has got a different meaning in its corresponding context.
    inline static void check_status(const int & status)
    {
        std::string lapack_errors(" ");

        if (status != 0)
        {
            if(status == -1)
                lapack_errors += "Something went wrong.."; //TODO: expand error list

            if (lapack_errors == " ")
                lapack_errors = "unknown LAPACK error state.";
#ifdef QT_NO_DEBUG
            AssertThrow(false, dealii::ExcMessage(lapack_errors.c_str() ) );
#else
            Assert(false, dealii::ExcMessage(lapack_errors.c_str() ) );
#endif
        }
    }

    //---------------------------------------------------------
   // @sect5{Copy functions}
   // @sect6{Function: lower_upper_triangle}
   //
   // Helper function for the LU-decomposition. Copies lower and upper triangular (trapezoidal) parts of $A$ into $L$ and $U$ using blas copy.
   // \param A : input Matrix of dimension $m \times n$.
   // \param L : output Matrix of dimension $m \times n$, contains strict lower triangular (trapezoidal) values of $A$ and unit diagonal.
   // \param U : output Matrix of dimension $m \times n$, contains upper triangular (trapezoidal) values of $A$.
    template <typename T>
    static void lower_upper_triangle(const SciPAL::Matrix<T, BW>&               A,
                                     SciPAL::Matrix<T, BW>&                     L,
                                     SciPAL::Matrix<T, BW>&                     U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m,n);

        for(size_t i = 0; i < n; i++)
        {
	    // Copy upper triangular part of A into U:
            BW::copy(((i + 1 < m) ? i + 1 : m),
                     A.data() + i * A.leading_dim,
                     A.stride,
                     U.data() + i * U.leading_dim,
                     U.stride);
            if (i < min_mn)
            {
		// Copy strict lower triangular part of A into L:
                BW::copy(m - (i + 1),
                         A.data() + i * A.leading_dim + i + 1,
                         A.stride,
                         L.data() + i * L.leading_dim + i + 1,
                         L.stride);
                L(i, i, T(1)); // Set diagonal entries of L to 1.0.
            }
        }
    }

   // @sect6{Function: lower_upper_triangle_elementwise}
   //
   // Helper function for the LU-decomposition. Copies lower and upper triangular (trapezoidal) parts of $A$ into $L$ and $U$ using element wise copy. It's not recommended to use this function other than for debugging purposes!
   // \param A : input Matrix of dimension $m \times n$.
   // \param L : output Matrix of dimension $m \times n$, contains strict lower triangular (trapezoidal) values of $A$ and unit diagonal.
   // \param U : output Matrix of dimension $m \times n$, contains upper triangular (trapezoidal) values of $A$.
    template <typename T>
    static void lower_upper_triangle_elementwise(const SciPAL::Matrix<T, BW>&   A,
                                                 SciPAL::Matrix<T, BW>&         L,
                                                 SciPAL::Matrix<T, BW>&         U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        for(size_t i = 0; i < m; i++)
        {
            for(size_t j = 0; j < n; j++)
            {
                if(i > j)
                    L(i, j, A(i,j)); 	// Copy strict lower triangular part of A into L.
                else
                    U(i, j, A(i,j));    // Copy upper triangular part of A into U.
                if(i == j)        	// Set diagonal entries of L to 1.0.
                    L(i, i, (T)1);
            }
        }
    }
   // @sect6{Function: upper_triangle}
   //
   // Helper function. Copies upper triangular (trapezoidal) parts of $A$ into $U$ using blas copy.
   // \param A : input Matrix of dimension $m \times n$.
   // \param U : output Matrix of dimension $m \times n$, contains upper triangular (trapezoidal) values of $A$.
    template <typename T>
    static void upper_triangle(const SciPAL::Matrix<T, BW>&                     A,
                               SciPAL::Matrix<T, BW>&                           U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();

	// Copy upper triangular part of A into U:
        for(size_t i = 0; i < n; i++)
            BW::copy(((i + 1 < m) ? i + 1 : m),
                     A.data() + i * A.leading_dim,
                     A.stride,
                     U.data() + i * U.leading_dim,
                     U.stride);
    }

   // @sect6{Function: lower_triangle}
   //
   // Helper function. Copies lower triangular (trapezoidal) parts of $A$ into $L$ using blas copy.
   // \param A : input Matrix of dimension $m \times n$.
   // \param L : output Matrix of dimension $m \times n$, contains strict lower triangular (trapezoidal) values of $A$ and unit diagonal.
    template <typename T>
    static void lower_triangle(const SciPAL::Matrix<T, BW>&                     A,
                               SciPAL::Matrix<T, BW>&                           L)
    {
        const size_t m = A.n_rows(), n = A.n_cols();
        const size_t min_mn = std::min(m, n);

        // Copy lower triangular part of A into L (assume the diagonal elements are 1).
        // Because we copy column-wise we need as much copy operation as we have columns.

        for(size_t i = 0; i < min_mn; i++)
        {
            BW::copy(m - i - 1,
                     A.data() + i * A.leading_dim + i + 1,
                     A.stride,
                     L.data() + i * L.leading_dim + i + 1,
                     L.stride);
            L(i, i, T(1)); // Set diagonal entries of L to 1.0.
        }

    }
   // @sect6{Function: upper_triangle_elementwise}
   //
   // Copies lower and upper triangular (trapezoidal) parts of $A$ into $L$ and $U$ using element wise copy. It's not recommended to use this function other than for debugging purposes!
   // \param A : input Matrix of dimension $m \times n$.
   // \param U : output Matrix of dimension $m \times n$, contains upper triangular (trapezoidal) values of $A$.
    template <typename T>
    static void upper_triangle_elementwise(const SciPAL::Matrix<T, BW>&         A,
                                           SciPAL::Matrix<T, BW>&               U)
    {
        const size_t m = A.n_rows(), n = A.n_cols();

	// Copy upper triangular part of A into U:
        for(size_t i = 0; i < m; i++)
            for(size_t j = 0; j < n; j++)
                if(i <= j)
                    U(i, j, A(i, j));
    }


    //---------------------------------------------------------
   // @sect5{Lapack function wrappers}
   // @sect6{Wrapper: SVD lapack-wrapper}
   // Wraps the general matrix svd lapack function for the four types s, d, c and z.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/d8/d49/sgesvd_8f.html
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

    // @sect6{Wrapper: SVD divide and conquer lapack-wrapper}
    // Wraps the general matrix sdd lapack function for the four types s, d, c and z.
     // Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga76f797b6a9e278ad7b21aae2b4a55d76.html#ga76f797b6a9e278ad7b21aae2b4a55d76
     inline static void gesdd(char * JOBZ, const size_t * M, const size_t * N, float * A, const size_t * LDA, float * S, float * U, const size_t * LDU, float * VT, const size_t * LDVT, float * WORK, const size_t * LWORK, float * /*RWORK*/, int * IWORK, int * INFO)
     {
         return sgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);
     }

     inline static void gesdd(char * JOBZ, const size_t * M, const size_t * N, double * A, const size_t * LDA, double * S, double * U, const size_t * LDU, double * VT, const size_t * LDVT, double * WORK, const size_t * LWORK, double * /*RWORK*/, int * IWORK, int * INFO)
     {
         return dgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);
     }

     inline static void gesdd(char * JOBZ, const size_t * M, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, float * S, SciPAL::CudaComplex<float> * U, const size_t * LDU, SciPAL::CudaComplex<float> * VT, const size_t * LDVT, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, float * RWORK, int * IWORK, int * INFO)
     {
         return cgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);
     }

     inline static void gesdd(char * JOBZ, const size_t * M, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, double * S, SciPAL::CudaComplex<double> * U, const size_t * LDU, SciPAL::CudaComplex<double> * VT, const size_t * LDVT, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, double * RWORK, int * IWORK, int * INFO)
     {
         return zgesdd_(JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);
     }

   // @sect6{Wrapper: LUD lapack-wrapper}
   // Wraps the general matrix LU decomposition lapack function for the four types s, d, c and z.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/de/de2/sgetrf_8f.html
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

   // @sect6{Wrapper: QRF lapack-wrapper}
   // Wraps the general matrix QR factorization lapack function for the four types s, d, c and z.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/df/d97/sgeqrf_8f.html
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

    // @sect6{Wrapper: QP3 lapack-wrapper}
    // Wraps the rank revealing matrix QR factorization lapack function for the four types s, d, c and z.
     // Lapack doc: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
     inline static void geqp3(const size_t * M, const size_t * N, float * A, const size_t * LDA, int * JPVT, float * TAU, float * WORK, const size_t * LWORK, int * INFO)
     {
         sgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO);
     }

     inline static void geqp3(const size_t * M, const size_t * N, double * A, const size_t * LDA, int * JPVT, double * TAU, double * WORK, const size_t * LWORK, int * INFO)
     {
         dgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO);
     }

     template <typename NT=float, typename BW=blas>
     inline static void geqp3(const size_t * M, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, int * JPVT, SciPAL::CudaComplex<float> * TAU, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, int * INFO)
     {
         SciPAL::Vector<NT, BW> RWORK(2*(*N));
         cgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK.data(), INFO);
     }

     template <typename NT=double, typename BW=blas>
     inline static void geqp3(const size_t * M, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, int * JPVT, SciPAL::CudaComplex<double> * TAU, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, int * INFO)
     {
         SciPAL::Vector<NT, BW> RWORK(2*(*N));
         zgeqp3_(M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK.data(), INFO);
     }

   // @sect6{Wrapper: MQR lapack-wrapper}
   // Wraps the lapack function of the product of elementary reflectors from the QR factorization for the four types s, d, c and z.
   // Lapack doc: http://www.netlib.org/lapack/explore-html/d0/d98/sormqr_8f.html
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

    // @sect6{Wrapper: EV lapack-wrapper}
    // Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the four types s, d, c and z.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html
    inline static void xxev(char * JOBZ, char * UPLO, const size_t * N, float * A, const size_t * LDA, float * W, float * WORK, const size_t * LWORK, float * /*RWORK*/, int * INFO)
    {
        return ssyev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
    }

    inline static void xxev(char * JOBZ, char * UPLO, const size_t * N, double * A, const size_t * LDA, double * W, double * WORK, const size_t * LWORK, double * /*RWORK*/, int * INFO)
    {
        return dsyev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO);
    }

    inline static void xxev(char * JOBZ, char * UPLO, const size_t * N, SciPAL::CudaComplex<float> * A, const size_t * LDA, float * W, SciPAL::CudaComplex<float> * WORK, const size_t * LWORK, float * RWORK, int * INFO)
    {
        return cheev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO);
    }

    inline static void xxev(char * JOBZ, char * UPLO, const size_t * N, SciPAL::CudaComplex<double> * A, const size_t * LDA, double * W, SciPAL::CudaComplex<double> * WORK, const size_t * LWORK, double * RWORK, int * INFO)
    {
        return zheev_(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO);
    }

    // @sect6{Wrapper: EV lapack-wrapper}
    // Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the two real types s and d.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html
    inline static void stev(char * JOBZ, const size_t * N, float * D, float * E, float * Z, const size_t * LDZ, float * WORK, int * INFO)
    {
        return sstev_(JOBZ, N, D, E, Z, LDZ, WORK, INFO);
    }

    inline static void stev(char * JOBZ, const size_t * N, double * D, double * E, double * Z, const size_t * LDZ, double * WORK, int * INFO)
    {
        return dstev_(JOBZ, N, D, E, Z, LDZ, WORK, INFO);
    }

    // @sect6{Wrapper: Matrix norm lapack-wrapper}
    // Wraps the lapack function for computing the eigenvalues and, optionally, the left and/or right eigenvectors for the two real types s and d.
    // Lapack doc: http://www.netlib.org/lapack/explore-html/d1/d74/group__eigen_s_y.html and http://www.netlib.org/lapack/explore-html/df/db2/cheev_8f.html

    inline static float lange(char * NORM, const size_t * M, const size_t * N, const float * A, const size_t * LDA, float * WORK)
    {
        return slange_( NORM, M, N, A, LDA, WORK );
    }

    inline static double lange(char * NORM, const size_t * M, const size_t * N, const double * A, const size_t * LDA, double * WORK)
    {
        return dlange_( NORM, M, N, A, LDA, WORK );
    }

    inline static float lange(char * NORM, const size_t * M, const size_t * N, const SciPAL::CudaComplex<float> * A, const size_t * LDA, float * WORK)
    {
        return clange_( NORM, M, N, A, LDA, WORK );
    }

    inline static double lange(char * NORM, const size_t * M, const size_t * N, const SciPAL::CudaComplex<double> * A, const size_t * LDA, double * WORK)
    {
        return zlange_( NORM, M, N, A, LDA, WORK );
    }
}; //END struct lapack
} //END namespace SciPAL
#endif // LAPACK_WRAPPER_HH
