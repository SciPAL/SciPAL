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


#ifndef BLAS_WRAPPER_HH
#define BLAS_WRAPPER_HH

#ifdef __APPLE__
//! On Apple systems BLAS is part of the Accelerate framework.
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif



//! deal.II components
#include <deal.II/base/exceptions.h>

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <vector_types.h> //! für float2, double2 Datentyp
#include <base/ParallelArch.h>

// @sect3{struct blas: cblas S,D,C,Z wrapper functions}
//!
//! We remove the type dependence of the function of CBLAS by using C++'s polymorphism and templates.
//!
//! Naming convention: The cblas_ prefix and the single character code for the number is omitted.
//! For instance: cblas_sgemv becomes gemv.
//! For the numerical purpose of the functions have a look at the BLAS doc.
struct blas {

    //! Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "cpu-blas"; }

    //! Compile-time variable for identifying the side of PCIe bus this BLAS works on.
    static const ParallelArch arch = cpu;

 #ifdef dfkgskdd

    // @sect4{Funktion: check_status}
    //!
    //! Wertet das cublasStatus Argument aus und wirft im Falle
    //! eines Fehlers eine Exception aus.
    static void  check_status(const cublasStatus & status)
    {

        std::string cublas_errors(" ");

        if (status != CUBLAS_STATUS_SUCCESS)
        {
            if (status == CUBLAS_STATUS_NOT_INITIALIZED)
                cublas_errors += "cublas not initialized ";
            if (status == CUBLAS_STATUS_MAPPING_ERROR)
                cublas_errors +="mapping error ";
            if (status == CUBLAS_STATUS_INVALID_VALUE)
                cublas_errors +="invalid value ";
            if (status == CUBLAS_STATUS_ALLOC_FAILED)
                cublas_errors +="allocation failed ";
            if (status == CUBLAS_STATUS_ARCH_MISMATCH)
                cublas_errors +="architecture mismatch ";
            if (status == CUBLAS_STATUS_EXECUTION_FAILED)
                cublas_errors +="execution failed ";
            if (status == CUBLAS_STATUS_INTERNAL_ERROR)
                cublas_errors +="cublas internal error ";

            if (cublas_errors == " ")
                cublas_errors = "unknown cublas error state";

#ifdef QT_NO_DEBUG
       AssertThrow(false, dealii::ExcMessage(cublas_errors.c_str() ) );
#else
       AssertThrow(false, dealii::ExcMessage(cublas_errors.c_str() ) );
#endif
        }

    }
#endif


    //! Initialize cblas. In this particular case this function does nothing.
    static void  Init() {

        std::cout << "ATLAS Init" << std::endl;
    }


    //! Shutdown cblas. In this particular case this function does nothing.
    static void Shutdown() {

        std::cout << "ATLAS Down" << std::endl;
    }


    //! Basic memory management. Encapsulates allocation, deallocation and resizing the array.
    template<typename T>
    class Data {

         //! Pointer to the first element of the array.
        T * __data;

         //! Number of elements in the array.
        size_t __n_el;

    public:
        //! This attribute can be used to determine an optimal
        //! value for the leading dimension by (re)allocating memory
        //! in multiples of it.
        //! For CUDA 32 is a reasonable value because this corresponds to the number of
        //! threads in a CUDA warp. For floats this is also the number of
        //! entries in a cache line. For the CPU this should be reasonable as well.
        static const int leading_dim_multiplier = 32;

        //! Default constructor. Sets up nothing.
        Data() : __data(0)
        {}

        //! Construct an array of length @p n.
        //! @param n : Number of elements to allocate.
        Data(size_t n)
            : __data(0), __n_el(0)
        {
            resize(n);
        }

        //! Resize the array to length @p n. All previous data is erased.
        //! @param n : New number of elements.
        void resize(size_t n)
        {
#ifdef QT_NO_DEBUG
            AssertThrow(n > 0,
                        dealii::ExcMessage("allocation of 0 elements not allowed"));
#else
            Assert(n > 0,
                   dealii::ExcMessage("allocation of 0 elements not allowed"));
#endif

            if (__n_el == 0)
            {
                __data = new T[n];
                __n_el = n;
            }
            else {
            if (__n_el != n)
                delete __data;
                __data = new T[n];
                __n_el = n;
            }

            if(__n_el == n)
                for(size_t ii = 0; ii <  __n_el; ii++)
                    __data[ii] = T();

        };

        ~Data()
        {
            if (__n_el > 0)
            {
                delete __data;
                __data = 0;
            }
        }

        //! Read-Write access to pointer to the first element of the array.
        T * data() { return __data; }

        //! Read access to pointer to the first element of the array.
        const T * data() const { return __data; }
    };


public:
    // @sect4{Funktion: SetMatrix}
    //!
    //! @param rows : Anzahl Zeilen.
    //! @param cols : Anzahl Spalten.
    //! @param A : Quellmatrix A.
    //! @param lda : leading dimension von A.
    //! @param B : Zielmatrix B.
    //! @param ldb : leading dimension von B.
    template<typename T,  typename T2>
    static void SetMatrix(int rows, int cols, const T2 *const &A,
                   int lda, T *&B, int ldb)
    {
        //! cublasStatus status = cublasSetMatrix(rows, cols, sizeof(T),
        //!                                      A, lda, B, ldb);

        copy(rows*cols, (A), 1, reinterpret_cast< T2*>(B), 1);

        //! check_status(status);
    }

    // @sect4{Funktion: GetMatrix}
    //!
    //! @param rows : Anzahl Zeilen.
    //! @param cols : Anzahl Spalten.
    //! @param A : Quellmatrix A.
    //! @param lda : leading dimension von A.
    //! @param B : Zielmatrix B.
    //! @param ldb : leading dimension von B.
    template<typename T>
    static void GetMatrix(int rows, int cols, const T * const &A,
                   int lda, T *&B, int ldb)
    {
        //! cublasStatus status = cublasGetMatrix(rows, cols, sizeof(T),
        //!                                      A, lda, B, ldb);

        //!

        for (int c = 0; c < cols; c++)
        {
            copy(rows, A+c*lda, 1, B+c*ldb, 1);


        }

        //! check_status(status);
    }

    // @sect4{Funktion: SetVector}
    //!
    //! @param n_el : Anzahl Elemente.
    //! @param src : Quellvektor src.
    //! @param inc_src : Speicher Abstand zwischen Elemente in Vector src.
    //! @param dst : Zielvektor dst.
    //! @param inc_dst : Speicher Abstand zwischen Elemente in Vector dst.
    template<typename T>
    static void SetVector(int n_el, const T * const src, int inc_src, T *dst, int inc_dst)
    {
        //! cublasStatus status = cublasSetVector(n_el, sizeof(T),
        //!                                      src, inc_src, dst, inc_dst);

        //! check_status(status);

        copy(n_el, src, inc_src, dst, inc_dst);

    }

    // @sect4{Funktion: GetVector}
    //!
    //! @param n_el : Anzahl Elemente.
    //! @param A: Quellvektor A.
    //! @param inc_src : Speicher Abstand zwischen Elemente in Vector A.
    //! @param B : Zielvektor B.
    //! @param inc_dst : Speicher Abstand zwischen Elemente in Vector B.
    template<typename T>
    static void GetVector(int n_el, const T * const &A, int inc_src, T *&B, int inc_dst)
    {
        //! cublasStatus status = cublasGetVector(n_el, sizeof(T),
        //!                                      A, inc_src, B, inc_dst);

        copy (n_el, A, inc_src, B, inc_dst);

        //! check_status(status);
    }

    // @sect4{Funktion: asum}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    static float
    asum (int n, const float *x,
          int incx)
    {
       float sum = cblas_sasum(n, x, incx);

       return sum;
    }

    static double
    asum (int n, const double *x,
          int incx)
    {
       double sum = cblas_dasum(n, x, incx);

       return sum;
    }

    static float
    asum (int n, const std::complex<float> *x,
          int incx)
    {
       float sum = cblas_scasum(n, reinterpret_cast<const float*>(x), incx);

       return sum;
    }

    static double
    asum (int n, const std::complex<double> *x,
          int incx)
    {
       double sum = cblas_dzasum(n, reinterpret_cast<const double*>(x), incx);

       return sum;
    }

    static float
    asum (int n, const float2 *x,
          int incx)
    {
       float sum = cblas_scasum(n, reinterpret_cast<const float*>(x), incx);

       return sum;
    }

    static double
    asum (int n, const double2 *x,
          int incx)
    {
       double sum = cblas_dzasum(n, reinterpret_cast<const double*>(x), incx);

       return sum;
    }

    // @sect4{Funktion: axpy}
    //!
    //! @param n : Anzahl Elemente.
    //! @param alpha: Skalar.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor (y= alpha*x+y).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.

    static void
    axpy (int n, float alpha, const float *x,
          int incx, float *y, int incy)
    {
       cblas_saxpy(n, alpha, x, incx, y, incy);

       //! cublas_status status = cublasGetError();

       //! check_status(status);
    }

    static void
    axpy (int n, double alpha, const double *x,
          int incx, double *y, int incy)
    {
       cblas_daxpy(n, alpha, x, incx, y, incy);

       //! cublasStatus status = cublasGetError();

       //! check_status(status);
    }

    static void
    axpy (int n, std::complex<double> alpha, const std::complex<double> *x,
          int incx, std::complex<double> *y, int incy)
    {
       cblas_zaxpy(n, reinterpret_cast<const double*>(&alpha), reinterpret_cast<const double*>(x), incx, reinterpret_cast<double*>(y), incy);
    }

    static void
    axpy (int n, std::complex<float> alpha, const std::complex<float> *x,
          int incx, std::complex<float> *y, int incy)
    {
       cblas_caxpy(n, reinterpret_cast<const float*>(&alpha), reinterpret_cast<const float*>(x), incx, reinterpret_cast<float*>(y), incy);
    }

    static void
    axpy (int n, double2 alpha, const double2 *x,
          int incx, double2 *y, int incy)
    {
       cblas_zaxpy(n, reinterpret_cast<const double*>(&alpha), reinterpret_cast<const double*>(x), incx, reinterpret_cast<double*>(y), incy);
    }

    static void
    axpy (int n, float2 alpha, const float2 *x,
          int incx, float2 *y, int incy)
    {
       cblas_caxpy(n, reinterpret_cast<const float*>(&alpha), reinterpret_cast<const float*>(x), incx, reinterpret_cast<float*>(y), incy);
    }

    // @sect4{Funktion: copy}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor (y= alpha*x+y).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.

    static void
    copy(int n, const float *x, int incx, float *y, int incy)
    {
        cblas_scopy(n, x, incx, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void
    copy(int n, const double *x, int incx, double *y, int incy)
    {
        cblas_dcopy(n, x, incx, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void
    copy(int n, const std::complex<double> *x, int incx, std::complex<double> *y, int incy)
    {
        cblas_zcopy(n, reinterpret_cast<const double*>(x), incx, reinterpret_cast<double*>(y), incy);
    }

    static void
    copy(int n, const std::complex<float> *x, int incx, std::complex<float> *y, int incy)
    {
        cblas_ccopy(n, reinterpret_cast<const float*>(x), incx, reinterpret_cast<float*>(y), incy);
    }

    static void
    copy(int n, const double2 *x, int incx, double2 *y, int incy)
    {
        cblas_zcopy(n, reinterpret_cast<const double*>(x), incx, reinterpret_cast<double*>(y), incy);
    }

    static void
    copy(int n, const float2 *x, int incx, float2 *y, int incy)
    {
        cblas_ccopy(n, reinterpret_cast<const float*>(x), incx, reinterpret_cast<float*>(y), incy);
    }

    // @sect4{Funktion: scal}
    //!
    //! @param n : Anzahl Elemente.
    //! @param alpha: Skalar.
    //! @param x : Quellvektor x (auch output).
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.


    static void
    scal (int n, float alpha, float *x, int incx)
    {
        cblas_sscal (n, alpha, x, incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void
    scal (int n, double alpha, double *x, int incx)
    {
        cblas_dscal (n, alpha, x, incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void
    scal (int n, std::complex<float> alpha, std::complex<float> *x, int incx)
    {
       cblas_cscal(n, reinterpret_cast<const float*>(&alpha), reinterpret_cast<float*>(x), incx);
    }

    static void
    scal (int n, std::complex<double>  alpha, std::complex<double> *x, int incx)
    {
       cblas_zscal(n, reinterpret_cast<const double*>(&alpha), reinterpret_cast<double*>(x), incx);
    }

    static void
    scal (int n, float2 alpha, float2 *x, int incx)
    {
       cblas_cscal(n, reinterpret_cast<const float*>(&alpha), reinterpret_cast<float*>(x), incx);
    }

    static void
    scal (int n, double2  alpha, double2 *x, int incx)
    {
       cblas_zscal(n, reinterpret_cast<const double*>(&alpha), reinterpret_cast<double*>(x), incx);
    }

    static void
    scal (int n, float alpha, float2 *x, int incx)
    {
        float2 a;
        a.x = alpha;
        a.y = 0;
       cblas_cscal(n, reinterpret_cast<const float*>(&a), reinterpret_cast<float*>(x), incx);
    }

    static void
    scal (int n, double  alpha, double2 *x, int incx)
    {
        double2 a;
        a.x = alpha;
        a.y = 0;
        cblas_zscal (n, reinterpret_cast<const double*>(&a), reinterpret_cast<double*>(x), incx);
    }



    // @sect4{Funktion: gemv}
    //!
    //! @param trans : gibt an ob A transponiert ist oder nicht. Sei trans = 'N' oder 'n' so ist op(A)= A, sei trans = 'T', 't','C' oder 'c' so ist op(A)= trans(A)
    //! @param m : Anzahl Zeilen in Matrix A.
    //! @param n : Anzahl Spalten in Matrix A.
    //! @param alpha: Skalar fuer A.
    //! @param A : Matrix A
    //! @param lda : leading dimension von A.
    //! @param x : Vektor mit der laenge von mindestens (1+(n-1)*abs(incx)) falls trans = 'N' oder 'n', sonst mindestens der laenge (1+(m-1)*abs(incx)).
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param beta : Skalar fuer Vektor y.
    //! @param y : Vektor mit der laenge von mindestens (1+(n-1)*abs(incy)) falls trans = 'N' oder 'n', sonst mindestens der laenge (1+(m-1)*abs(incy)).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.

    static void gemv (char trans, int m, int n, float alpha,
                      const float * const A, int lda,
                      const float * const x,  int incx, float beta,
                      float *y, int incy)
    {
        CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
        cblas_sgemv (CblasColMajor, tr, m, n, alpha, A, lda, x, incx, beta, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void gemv (char trans, int m, int n, double alpha,
                      const double * const A, int lda,
                      const double * const x,  int incx, double beta,
                      double *y, int incy)
    {
        CBLAS_TRANSPOSE tr = ( ( (trans == 't') || (trans == 'T') ) ? CblasTrans : CblasNoTrans );
        cblas_dgemv (CblasColMajor, tr, m, n, alpha, A, lda, x, incx, beta, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    // @sect4{Funktion: ger}
    //!
    //! @param m : Anzahl Zeilen in Matrix A.
    //! @param n : Anzahl Spalten in Matrix A.
    //! @param alpha: Skalar fuer x*trans(y).
    //! @param x : Vektor mit der laenge von mindestens der laenge (1+(m-1)*abs(incx)).
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Vektor mit der laenge von mindestens (1+(n-1)*abs(incy)).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.
    //! @param A : Matrix A
    //! @param lda : leading dimension von A.

    static void
    ger(int m, int n, float alpha, const float *x,
        int incx, const float *y, int incy, float *A,
        int lda)
    {
         cblas_sger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);

         //! cublasStatus status = cublasGetError();

         //! check_status(status);
    }

    static void
    ger(int m, int n, double alpha, const double *x,
        int incx, const double *y, int incy, double *A,
        int lda)
    {
         cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);

         //! cublasStatus status = cublasGetError();

         //! check_status(status);
    }



    // @sect4{Funktion: gemm}
    //!
    //! @param transa : gibt an ob A transponiert ist oder nicht. Sei transa = 'N' oder 'n' so ist op(A)= A, sei transa = 'T' oder 't' so ist op(A)= trans(A), sei transa = 'C' oder 'c' so ist op(A)=adjoint(A)
    //! @param transb : gibt an ob B transponiert ist oder nicht. Sei transb = 'N' oder 'n' so ist op(B)= A, sei transb = 'T' oder 't' so ist op(B)= trans(B), sei transb = 'C' oder 'c' so ist op(B)=adjoint(B)
    //! @param m : Anzahl Zeilen in Matrix A und Matrix C.
    //! @param n : Anzahl Spalten in Matrix B und Matrix C.
    //! @param k : Anzahl Spalten in Matrix A und Zeilen in Matrix B.
    //! @param alpha: Skalar fuer op(A)*op(B).
    //! @param A : Matrix A
    //! @param lda : leading dimension von A.
    //! @param B : Matrix B.
    //! @param ldb : leading dimension von B.
    //! @param beta : Skalar fuer Matrix C.
    //! @param C : Matrix C.
    //! @param ldc : leading dimension von C.
    static void gemm(char transa, char transb, int m, int n, int k, float alpha,
                     const float * const A, int lda, const float * const B, int ldb,
                     float beta, float * C, int ldc)
    {
        CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : CblasNoTrans );
        CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : CblasNoTrans );

        cblas_sgemm(CblasColMajor,
                    tr_a, tr_b,
                    m, n, k,
                    alpha,
                    A, lda,
                    B, ldb,
                    beta,
                    C, ldc);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }


    static void gemm(char transa, char transb, int m, int n, int k, double alpha,
                     const double * const A, int lda, const double * const B, int ldb,
                     double beta, double * C, int ldc)
    {
        CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : CblasNoTrans );
        CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : CblasNoTrans );

        cblas_dgemm(CblasColMajor, tr_a, tr_b, m, n, k, alpha,
                    A, lda, B, ldb,
                    beta, C, ldc);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }

    static void gemm(char transa, char transb, int m, int n, int k, const float2 alpha,
                     const float2 * const A, int lda, const float2 * const B, int ldb,
                     const float2 beta, float2 * C, int ldc)
    {
        CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasConjTrans : CblasNoTrans );
        CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasConjTrans : CblasNoTrans );

        cblas_cgemm(CblasColMajor,
                    tr_a, tr_b,
                    m, n, k,
                    reinterpret_cast<const float*>(&alpha),
                    reinterpret_cast<const float*>(A), lda,
                    reinterpret_cast<const float*>(B), ldb,
                    reinterpret_cast<const float*>(&beta),
                    reinterpret_cast<float*>(C), ldc);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }


    static void gemm(char transa, char transb, int m, int n, int k, const double2 alpha,
                     const double2 * const A, int lda, const double2 * const B, int ldb,
                     const double2 beta, double2 * C, int ldc)
    {
        CBLAS_TRANSPOSE tr_a = ( ( (transa == 't') || (transa == 'T') ) ? CblasTrans : ( (transa == 'c') || (transa == 'C') ) ? CblasConjTrans : CblasNoTrans );
        CBLAS_TRANSPOSE tr_b = ( ( (transb == 't') || (transb == 'T') ) ? CblasTrans : ( (transb == 'c') || (transb == 'C') ) ? CblasConjTrans : CblasNoTrans );

        cblas_zgemm(CblasColMajor,
                    tr_a, tr_b,
                    m, n, k,
                    reinterpret_cast<const double*>(&alpha),
                    reinterpret_cast<const double*>(A), lda,
                    reinterpret_cast<const double*>(B), ldb,
                    reinterpret_cast<const double*>(&beta),
                    reinterpret_cast<double*>(C), ldc);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);
    }


    // @sect4{Funktion: nrm2}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.

    static float
    nrm2(int n, const float *x, int incx)
    {
        float result = cblas_snrm2 (n, x, incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static double
    nrm2(int n, const double *x, int incx)
    {
        double result = cblas_dnrm2 (n, x, incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static float
    nrm2(int n, const float2 *x, int incx)
    {
        float result = cblas_scnrm2 (n, reinterpret_cast<const float*>(x), incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static double
    nrm2(int n, const double2 *x, int incx)
    {
        double result = cblas_dznrm2 (n, reinterpret_cast<const double*>(x), incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static float
            nrm2(int n, const std::complex<float> *x, int incx)
    {
        float result = cblas_scnrm2 (n, reinterpret_cast<const float*>(x), incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static double
            nrm2(int n, const std::complex<double> *x, int incx)
    {
        double result = cblas_dznrm2 (n, reinterpret_cast<const double*>(x), incx);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    // @sect4{Funktion: dot}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor.
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.

    static float  dot(int n, const float *x, int incx, const float *y, int incy)
    {
        float result = cblas_sdot(n, x, incx, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

    static double dot(int n, const double *x,
                      int incx, const double *y, int incy)
    {
        double result = cblas_ddot(n, x, incx, y, incy);

        //! cublasStatus status = cublasGetError();

        //! check_status(status);

        return result;
    }

}; //! struct CBW END




#endif //! BLAS_WRAPPER_HH
