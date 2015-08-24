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


#ifndef CUBLAS_WRAPPER_HH
#define CUBLAS_WRAPPER_HH

//STL includes
#include <string>
#include <complex>
#include <iostream>
#include <typeinfo>

//include deal.II exceptions
#include <deal.II/base/exceptions.h>

//SciPAL includes
#include<base/ParallelArch.h>
#include<base/CUDA_error_check.h>

//CUDA includes
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

//! A global handle object. this is needed later by all cublas functions
static cublasHandle_t handle;

// @sect3{struct cublas: cublas S,D,C,Z wrapper functions}
//!
//! We remove the type dependence of the function of CBLAS by using C++'s polymorphism and templates.
//!
//! Naming convention: The cblas_ prefix and the single character code for the number is omitted.
//! For instance: cblas_sgemv becomes gemv.
//! For the numerical purpose of the functions have a look at the (CU)BLAS doc.
//!
struct cublas {

    typedef cublas blas_wrapper_type;
    // @sect4{Funktion: name}
    //!
    //! Return a string identifier. This is useful in the output of performance tests.
    static inline std::string name () { return "cublas"; }

    //! Compile-time variable for identifying the side of PCIe bus this BLAS works on.
    static const ParallelArch arch = gpu_cuda;

    //! TODO: set und get-Routinen verpacken

    // @sect4{Funktion: check_status}
    //!
    //! Evaluate the @p cublasStatus argument and throws
    //! an exception in case an error occurs.
    static void  check_status(const cublasStatus_t & status)
    {
        //! The following is for tracking cublas usage
        //! std::cout << "Checking cublas status" << std::endl;
#ifdef DEBUG
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
       Assert(false, dealii::ExcMessage(cublas_errors.c_str() ) );
#endif
        }
#endif
    }

    static void  check_status(const cudaError_t & status)
    {
        //! The following is for tracking cublas usage
        //! std::cout << "Checking cublas status" << std::endl;
#ifdef DEBUG
        std::string cuda_errors(" ");

        if (status != cudaSuccess)
        {
            cuda_errors += cudaGetErrorString(status);


#ifdef QT_NO_DEBUG
       AssertThrow(false, dealii::ExcMessage(cuda_errors.c_str() ) );
#else
       Assert(false, dealii::ExcMessage(cuda_errors.c_str() ) );
#endif
        }
#endif
    }




    // @sect4{Funktion: Init}
    //!
    //! Initialize cublas.
    //! If it fails an exception is thrown. To get a more readable output
    //! we use deal.II's Assert and AssertThrow macros. ONe could be even more
    //! precise by implementing a taylor-made exception when then could be tracked with C++'s
    //! try-catch clause. But that seems to be a bit of an overkill.
    static void  Init() {

        cublasStatus_t s = cublasCreate(&handle);

#ifdef DEBUG
        if (s == CUBLAS_STATUS_SUCCESS)
            std::cout << "cublas init succeeded" << std::endl;
#endif

#ifndef DEBUG
        AssertThrow(s == CUBLAS_STATUS_SUCCESS,
                    dealii::ExcMessage("cublas init failed"));
#else
        Assert(s == CUBLAS_STATUS_SUCCESS,
               dealii::ExcMessage("cublas init failed"));
#endif
    }


    // @sect4{Funktion: Shutdown}
    //!
    //! Shut down cublas. If it fails an exception is thrown.
    static void Shutdown() {

        cublasStatus_t s = cublasDestroy(handle);

#ifdef DEBUG
        if (s == CUBLAS_STATUS_SUCCESS)
            std::cout << "cublas shutdown succeeded" << std::endl;
#endif

#ifdef QT_NO_DEBUG
        AssertThrow(s == CUBLAS_STATUS_SUCCESS,
                    dealii::ExcMessage("cublas shutdown failed"));
#else
        Assert(s == CUBLAS_STATUS_SUCCESS,
                    dealii::ExcMessage("cublas shutdown failed"));
#endif

    }

    //! This structure encapsulates the basic memory management.
    template<typename T>
    struct Data {

        //! This attribute can be used to determine an optimal
        //! value for the leading dimension by (re)allocating memory
        //! in multiples of it.
        //! Its current value is 32 which corresponds to the number of
        //! threads in a CUDA warp. For floats this is also the number of
        //! entries in a cache line.
        static const int leading_dim_multiplier = 32;

        Data() : dev_ptr(0) {}

        Data(size_t n_rows, size_t n_cols = 1)
        {
            alloc(n_rows, n_cols);
        }

        ~Data() {  free_dev_ptr(); }

        T * data() { return dev_ptr; }

        const T * data() const { return dev_ptr; }

        void resize(size_t n_rows, size_t n_cols = 1)
        {
            free_dev_ptr();
            alloc(n_rows, n_cols);
        }

   /*     size_t leading_dim() const {
            size_t n_bytes = sizeof(T);

                if( pitch_in_bytes % n_bytes != 0 )
                {
                    std::cout << "something is wrong with the pitch pitch_in_bytes % n_bytes ="<< pitch_in_bytes % n_bytes<< std::endl;
                return 0;
                }
                else
                {
                    size_t result = pitch_in_bytes;//  / n_bytes;
                    return result;
                }
        }
        */

    protected:
        void swap(Data<T>& other)
        {
            std::swap(dev_ptr, other.dev_ptr);
        }

    private:
        void free_dev_ptr()
        {
            if (dev_ptr == 0) return;

             cudaError_t status = cudaFree(dev_ptr);
             check_status(status);

            dev_ptr = 0;
        }


        void alloc(size_t rows, size_t cols = 1)
        {
            //cols and rows switched on purpose to meet cublas alignment requirements!!!!!!!1
//            cudaError_t status = cudaMallocPitch((void**)&dev_ptr, &pitch_in_bytes,
//                                             rows, cols);
            //cols and rows switched on purpose to meet cublas alignment requirements!!!!!!!1

            cudaError_t status = cudaMalloc( (void**)&dev_ptr, rows*cols*sizeof(T) );
            check_status(status);

            // set everything to zero
            status = cudaMemset( dev_ptr, 0, rows*cols*sizeof(T) );
//            status = cudaMemset2D(dev_ptr, pitch_in_bytes,
//                                  0, cols*sizeof(T), rows*sizeof(T));
            // According to the
            // <a href="http://developer.download.nvidia.com/compute/cuda/4_2/rel/toolkit/docs/online/sync_async.html#memset_sync_async_behavior">online documentation</a>
            // cudaMemset is asynchronous.
            // Thus:
            cudaThreadSynchronize();
            // ans only then we check the error status.
            check_status(status);

            // Initially, we wanted to use cublas allocation routine. However, it seemed
            // to suffer from 32 bit issues as no more then 512MB could alocated with it.
            // cublasStatus_t status = cublasAlloc( n, sizeof(T), (void**)&dev_ptr);
            // check_status(status);
        }

        T * dev_ptr;
        size_t pitch_in_bytes;
    };


private:


public:



    // @sect4{Funktion: SetMatrix}
    //!
    //! @param rows : Anzahl Zeilen.
    //! @param cols : Anzahl Spalten.
    //! @param A : Quellmatrix A.
    //! @param lda : leading dimension von A.
    //! @param B : Zielmatrix B.
    //! @param ldb : leading dimension von B.
    template<typename T, typename T2>
    static void SetMatrix(int rows, int cols, const T2 *const A,
                   int lda, T *B, int ldb)
    {
        if(sizeof(T) != sizeof(T2))
            std::cout<<"You are trying to copy matrices with different Types: T="
                     << typeid(T).name() << ", T2=" << typeid(T2).name() ;

        cublasStatus_t status = cublasSetMatrix(rows, cols, sizeof(T),
                                                A, lda,
                                                B, ldb);
//        cudaError_t status = cudaMemcpy2D( B, ldb,
//                               A, lda,
//                               cols, rows, cudaMemcpyHostToDevice);

//        gpuErrchk( cudaPeekAtLastError() );
//        gpuErrchk( cudaDeviceSynchronize() );
        check_status(status);
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
    static void GetMatrix(int rows, int cols, const T * const A,
                   int lda, T *B, int ldb)
    {
        cublasStatus_t status = cublasGetMatrix(rows, cols, sizeof(T),
                                              A, lda, B, ldb);

        check_status(status);
    }

    // @sect4{Funktion: SetVector}
    //!
    //! @param n_el : Anzahl Elemente.
    //! @param src : Quellvektor src.
    //! @param inc_src : Speicher Abstand zwischen Elemente in Vector src.
    //! @param dst : Zielvektor dst.
    //! @param inc_dst : Speicher Abstand zwischen Elemente in Vector dst.
    template<typename T>
    static void SetVector(int n_el, const T * const src, int inc_src,
                          T *dst, int inc_dst)
    {
        cublasStatus_t status = cublasSetVector(n_el, sizeof(T),
                                              src, inc_src, dst, inc_dst);

        check_status(status);

    }

    // @sect4{Funktion: GetVector}
    //!
    //! @param n_el : Anzahl Elemente.
    //! @param A: Quellvektor A.
    //! @param inc_src : Speicher Abstand zwischen Elemente in Vector A.
    //! @param B : Zielvektor B.
    //! @param inc_dst : Speicher Abstand zwischen Elemente in Vector B.
    template<typename T>
    static void GetVector(int n_el, const T * const src, int inc_src,
                          T *dst, int inc_dst)
    {
        cublasStatus_t status = cublasGetVector(n_el, sizeof(T),
                                              src, inc_src, dst, inc_dst);

        check_status(status);
    }

    // @sect4{Function: Memset}
    //!
    //! @param n_el : Anzahl Elemente.
    //! @param A: Quellvektor A.
    //! @param inc_src : Speicher Abstand zwischen Elemente in Vector A.
    //! @param B : Zielvektor B.
    //! @param inc_dst : Speicher Abstand zwischen Elemente in Vector B.

    template<typename T>
    static void Memset(int n_el, T* ptr, int val)
    {

//      cublasStatus_t status ;
//      status =
              cudaMemset( ptr, val, n_el*sizeof(T) );
//         check_status(status);
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
       float sum;
       cublasStatus_t status = cublasSasum(handle, n, x, incx, &sum);
       check_status(status);
       return sum;
    }

    static double
    asum (int n, const double *x,
          int incx)
    {
        double sum;
        cublasStatus_t status = cublasDasum(handle, n, x, incx, &sum);
        check_status(status);
        return sum;
    }

    static float
    asum (int n, const cuComplex *x,
          int incx)
    {

        float sum;
        cublasStatus_t status = cublasScasum(handle, n, x, incx, &sum);
        check_status(status);
        return sum;

    }

    static double
    asum (int n, const cuDoubleComplex *x,
          int incx)
    {
        double sum;
        cublasStatus_t status = cublasDzasum(handle, n, x, incx, &sum);
        check_status(status);
        return sum;

    }

    // @sect4{Funktion: amin}
    //!
    //! This function returns the index of entry with the minimal value.
    //!
    //! @param n : number of elements.
    //! @param x : source vector x.
    //! @param incx : stride between consecutive elements of x.

    static int
    amin (int n, const float *x,
          int incx)
    {
       int result;
       cublasStatus_t status = cublasIsamin(handle, n, x, incx, &result);
       check_status(status);
       return result - 1;
    }

    static int
    amin (int n, const double *x,
          int incx)
    {
        int result;
        cublasStatus_t status = cublasIdamin(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;
    }

    static int
    amin (int n, const cuComplex *x,
          int incx)
    {

        int result;
        cublasStatus_t status = cublasIcamin(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;

    }

    static int
    amin (int n, const cuDoubleComplex *x,
          int incx)
    {
        int result;
        cublasStatus_t status = cublasIzamin(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;

    }


    // @sect4{Funktion: amax}
    //!
    //! @param n : number of elements.
    //! @param x : source vector x.
    //! @param incx : stride between consecutive elements of x.

    static int
    amax (int n, const float *x,
          int incx)
    {
       int result;
       cublasStatus_t status = cublasIsamax(handle, n, x, incx, &result);
       check_status(status);
       return result - 1;
    }

    static int
    amax (int n, const double *x,
          int incx)
    {
        int result;
        cublasStatus_t status = cublasIdamax(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;
    }

    static int
    amax (int n, const cuComplex *x,
          int incx)
    {

        int result;
        cublasStatus_t status = cublasIcamax(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;

    }

    static int
    amax (int n, const cuDoubleComplex *x,
          int incx)
    {
        int result;
        cublasStatus_t status = cublasIzamax(handle, n, x, incx, &result);
        check_status(status);
        return result - 1;

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
       cublasStatus_t status = cublasSaxpy(handle, n, &alpha, x, incx, y, incy);
       check_status(status);
    }

    static void
    axpy (int n, double alpha, const double *x,
          int incx, double *y, int incy)
    {

        cublasStatus_t status = cublasDaxpy(handle, n, &alpha, x, incx, y, incy);
        check_status(status);

    }

    static void
    axpy (int n, cuComplex alpha, const cuComplex *x,
          int incx, cuComplex *y, int incy)
    {
        cublasStatus_t status = cublasCaxpy(handle, n, &alpha, x, incx, y, incy);
        check_status(status);
    }

    static void
    axpy (int n, cuDoubleComplex alpha, const cuDoubleComplex *x,
          int incx, cuDoubleComplex *y, int incy)
    {
        cublasStatus_t status = cublasZaxpy(handle, n, &alpha, x, incx, y, incy);
        check_status(status);
    }

    // @sect4{Funktion: copy}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor (y= alpha*x+y).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.


    /*
    template <typename T>
    struct fct_ptrs;

    template<>
    struct<float> {

        typedef (*cublasScopy) ptr_t;

    };


    template<>
    struct<double> {

        typedef double T;

        typedef (*cublasDcopy)(int, T*, int, T*, int) copy_ptr_t;
        typedef (*cublasDscal) scal_ptr_t;

    };



    static void
    template<typename T>
    copy(int n, const T *x, int incx, T *y, int incy)
    {
        cublasHandle_t h;
        cublasStatus_t status = cublasCreate(&h);
        check_status(status);


        (*(typename fct_ptrs<T>::ptr_t))(h, n, x, incx, y, incy);

        check_status(status);
        status = cublasDestroy(h);
        check_status(status);

    }


    static void
    copy_impl(int n, const float *x, int incx, float *y, int incy)
    {
        cublasStatus_t status = cublasScopy(handle, n, x, incx, y, incy);
        check_status(status);
    }
     */

    static void
    copy(int n, const float *x, int incx, float *y, int incy)
    {
        cublasStatus_t status = cublasScopy(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    copy(int n, const double *x, int incx, double *y, int incy)
    {
        cublasStatus_t status;
        status = cublasDcopy(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    copy(int n, const cuComplex *x, int incx, cuComplex *y, int incy)
    {
        cublasStatus_t status = cublasCcopy(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    copy(int n, const cuDoubleComplex *x, int incx, cuDoubleComplex *y, int incy)
    {
        cublasStatus_t status = cublasZcopy(handle, n, x, incx, y, incy);
        check_status(status);
    }

    // @sect4{Funktion: swap}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor (y= alpha*x+y).
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.

    static void
    swap(int n, float *x, int incx, float *y, int incy)
    {
        cublasStatus_t status = cublasSswap(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    swap(int n, double *x, int incx, double *y, int incy)
    {
        cublasStatus_t status;
        status = cublasDswap(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    swap(int n, cuComplex *x, int incx, cuComplex *y, int incy)
    {
        cublasStatus_t status = cublasCswap(handle, n, x, incx, y, incy);
        check_status(status);
    }

    static void
    swap(int n, cuDoubleComplex *x, int incx, cuDoubleComplex *y, int incy)
    {
        cublasStatus_t status = cublasZswap(handle, n, x, incx, y, incy);
        check_status(status);
    }

    // @sect4{Funktion:}
    //!
    //! @param n : Anzahl Elemente.
    //! @param alpha: Skalar.
    //! @param x : Quellvektor x (auch output).
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    static void
    scal (int n, float alpha, float *x, int incx)
    {
        cublasStatus_t status = cublasSscal (handle, n, &alpha, x, incx);
        check_status(status);
    }

    static void
    scal (int n, double alpha, double *x, int incx)
    {
        cublasStatus_t status = cublasDscal (handle, n, &alpha, x, incx);
        check_status(status);
    }


    static void
            scal (int n, cuComplex alpha, cuComplex *x, int incx)
    {
        cublasStatus_t status = cublasCscal (handle, n, &alpha, x, incx);
        check_status(status);
    }



    static void
    scal (int n, cuDoubleComplex alpha, cuDoubleComplex *x, int incx)
    {
        cublasStatus_t status = cublasZscal (handle, n, &alpha, x, incx);
        check_status(status);
    }



    static void
    scal (int n, std::complex<float> alpha, cuComplex *x, int incx)
    {
        cuComplex a;
        a.x = alpha.real();
        a.y = alpha.imag();

        cublasStatus_t status = cublasCscal (handle, n, &a, x, incx);
        check_status(status);
    }



    static void
    scal (int n, std::complex<double> alpha, cuDoubleComplex *x, int incx)
    {
        cuDoubleComplex a;
        a.x = alpha.real();
        a.y = alpha.imag();

        cublasStatus_t status = cublasZscal (handle, n, &a, x, incx);
        check_status(status);
    }

    // @sect4{Funktion: geam}
    //! C = α op ( A ) + β op ( B )
    //! op(A) = A 	if  transa == CUBLAS_OP_N
    //!  transpose(A) if  transa == CUBLAS_OP_T
    //!  hermitian(A) = if  transa == CUBLAS_OP_C
    //! ldx = leading dimension

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          float  *alpha,
          const float *A, int lda,
          const float *beta,
          const float *B, int ldb,
          float *C, int ldc )
    {
    cublasStatus_t status = cublasSgeam(handle,
                                        opA,  opB,
                                        rows,  cols,
                                        alpha, A, lda,
                                        beta, B, ldb,
                                        C, ldc);
    check_status(status);
    }

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          const double  *alpha,
          const double *A, int lda,
          const double *beta,
          const double *B, int ldb,
          double *C, int ldc )
    {
    cublasStatus_t status = cublasDgeam(handle,
                                        opA,  opB,
                                        rows,  cols,
                                        alpha, A, lda,
                                        beta, B, ldb,
                                        C, ldc);
    check_status(status);
    }

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          const cuComplex  *alpha,
          const cuComplex *A, int lda,
          const cuComplex *beta,
          const cuComplex *B, int ldb,
          cuComplex *C, int ldc )
    {
    cublasStatus_t status = cublasCgeam(handle,
                                        opA,  opB,
                                        rows,  cols,
                                        alpha, A, lda,
                                        beta, B, ldb,
                                        C, ldc);
    check_status(status);
    }

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          const cuDoubleComplex  *alpha,
          const cuDoubleComplex *A, int lda,
          const cuDoubleComplex *beta,
          const cuDoubleComplex *B, int ldb,
          cuDoubleComplex *C, int ldc )
    {
    cublasStatus_t status = cublasZgeam(handle,
                                        opA,  opB,
                                        rows,  cols,
                                        alpha, A, lda,
                                       beta, B, ldb,
                                        C, ldc);
    check_status(status);
    }

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          const std::complex<float>  *alpha,
          const cuComplex *A, int lda,
          const std::complex<float> *beta,
          const cuComplex *B, int ldb,
          cuComplex *C, int ldc )
    {
        cuComplex _alpha;
        _alpha.x = (*alpha).real();
        _alpha.y = (*alpha).imag();

        cuComplex _beta;
        _beta.x = (*beta).real();
        _beta.y = (*beta).imag();

        cublasStatus_t status = cublasCgeam(handle,
                                            opA,  opB,
                                            rows,  cols,
                                            &_alpha, A, lda,
                                            &_beta, B, ldb,
                                            C, ldc);
        check_status(status);
    }

    static void
    geam (cublasOperation_t opA, cublasOperation_t opB, int rows, int cols,
          const std::complex<double>  *alpha,
          const cuDoubleComplex *A, int lda,
          const std::complex<double> *beta,
          const cuDoubleComplex *B, int ldb,
          cuDoubleComplex *C, int ldc )
    {
        cuDoubleComplex _alpha;
        _alpha.x = (*alpha).real();
        _alpha.y = (*alpha).imag();

        cuDoubleComplex _beta;
        _beta.x = (*beta).real();
        _beta.y = (*beta).imag();

        cublasStatus_t status = cublasZgeam(handle,
                                            opA,  opB,
                                            rows,  cols,
                                            &_alpha, A, lda,
                                            &_beta, B, ldb,
                                            C, ldc);
        cudaThreadSynchronize();
        check_status(status);
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

        cublasOperation_t transaction;

        if(trans=='N' || trans=='n')
            transaction = CUBLAS_OP_N;

        if(trans=='T' || trans=='t')
            transaction = CUBLAS_OP_T;

        if(trans=='C' || trans=='c')
            transaction = CUBLAS_OP_C;

        cublasStatus_t status = cublasSgemv (handle, transaction, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
        check_status(status);
    }

    static void gemv (char trans, int m, int n, double alpha,
                      const double * const A, int lda,
                      const double * const x,  int incx, double beta,
                      double *y, int incy)
    {
        cublasOperation_t transaction;

        if(trans=='N' || trans=='n')
            transaction = CUBLAS_OP_N;

        if(trans=='T' || trans=='t')
            transaction = CUBLAS_OP_T;

        if(trans=='C' || trans=='c')
            transaction = CUBLAS_OP_C;

        cublasStatus_t status = cublasDgemv (handle, transaction, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
        check_status(status);
    }

    static void gemv (char trans, int m, int n, cuComplex alpha,
                      const cuComplex * const A, int lda,
                      const cuComplex * const x,  int incx, cuComplex beta,
                      cuComplex *y, int incy)
    {
        cublasOperation_t transaction;

        if(trans=='N' || trans=='n')
            transaction = CUBLAS_OP_N;

        if(trans=='T' || trans=='t')
            transaction = CUBLAS_OP_T;

        if(trans=='C' || trans=='c')
            transaction = CUBLAS_OP_C;

        cublasStatus_t status = cublasCgemv (handle, transaction, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
        check_status(status);
    }

    static void gemv (char trans, int m, int n, cuDoubleComplex alpha,
                      const cuDoubleComplex * const A, int lda,
                      const cuDoubleComplex * const x,  int incx, cuDoubleComplex beta,
                      cuDoubleComplex *y, int incy)
    {
        cublasOperation_t transaction;

        if(trans=='N' || trans=='n')
            transaction = CUBLAS_OP_N;

        if(trans=='T' || trans=='t')
            transaction = CUBLAS_OP_T;

        if(trans=='C' || trans=='c')
            transaction = CUBLAS_OP_C;

        cublasStatus_t status = cublasZgemv (handle, transaction, m, n, &alpha, A, lda, x, incx, &beta, y, incy);
        check_status(status);
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
         cublasStatus_t status = cublasSger(handle, m, n, &alpha, x, incx, y, incy, A, lda);
         check_status(status);
    }

    static void
    ger(int m, int n, double alpha, const double *x,
        int incx, const double *y, int incy, double *A,
        int lda)
    {
        cublasStatus_t status = cublasDger(handle, m, n, &alpha, x, incx, y, incy, A, lda);
        check_status(status);
    }



    // @sect4{Funktion: gemm}
    //!
    //! @param transa : gibt an ob A transponiert ist oder nicht. Sei transa = 'N' oder 'n' so ist op(A)= A, sei trans = 'T', 't','C' oder 'c' so ist op(A)= trans(A)
    //! @param transb : gibt an ob B transponiert ist oder nicht. Sei transb = 'N' oder 'n' so ist op(B)= A, sei trans = 'T', 't','C' oder 'c' so ist op(B)= trans(B)
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
    static void gemm(char transa, char transb,
                     int m, int n, int k, float alpha,
                     const float * const A, int lda,
                     const float * const B, int ldb,
                     float beta, float * C, int ldc)
    {

        cublasOperation_t transactionA, transactionB;

        if(transa=='N' || transa=='n')
            transactionA = CUBLAS_OP_N;

        if(transa=='T' || transa=='t')
            transactionA = CUBLAS_OP_T;

        if(transa=='C' || transa=='c')
            transactionA = CUBLAS_OP_C;


        if(transb=='N' || transb=='n')
            transactionB = CUBLAS_OP_N;

        if(transb=='T' || transb=='t')
            transactionB = CUBLAS_OP_T;

        if(transb=='C' || transb=='c')
            transactionB = CUBLAS_OP_C;

        cublasStatus_t status = cublasSgemm(handle, transactionA, transactionB, m, n, k, &alpha,
                    A, lda, B, ldb,
                    &beta, C, ldc);

        check_status(status);
    }


    static void gemm(char transa, char transb, int m, int n, int k, double alpha,
                     const double * const A, int lda, const double * const B, int ldb,
                     double beta, double * C, int ldc)
    {
        cublasOperation_t transactionA, transactionB;

        if(transa=='N' || transa=='n')
            transactionA = CUBLAS_OP_N;

        if(transa=='T' || transa=='t')
            transactionA = CUBLAS_OP_T;

        if(transa=='C' || transa=='c')
            transactionA = CUBLAS_OP_C;


        if(transb=='N' || transb=='n')
            transactionB = CUBLAS_OP_N;

        if(transb=='T' || transb=='t')
            transactionB = CUBLAS_OP_T;

        if(transb=='C' || transb=='c')
            transactionB = CUBLAS_OP_C;

        cublasStatus_t status = cublasDgemm(handle,
                                            transactionA,
                                            transactionB,
                                            m, n, k, &alpha,
                                            A, lda, B, ldb,
                                            &beta, C, ldc);


        check_status(status);
    }

    static void gemm(char transa, char transb, int m, int n, int k, cuComplex alpha,
                         const cuComplex * const A, int lda, const cuComplex* const B, int ldb,
                         cuComplex beta, cuComplex * C, int ldc)
        {
            cublasOperation_t transactionA, transactionB;

            if(transa=='N' || transa=='n')
                transactionA = CUBLAS_OP_N;

            if(transa=='T' || transa=='t')
                transactionA = CUBLAS_OP_T;

            if(transa=='C' || transa=='c')
                transactionA = CUBLAS_OP_C;


            if(transb=='N' || transb=='n')
                transactionB = CUBLAS_OP_N;

            if(transb=='T' || transb=='t')
                transactionB = CUBLAS_OP_T;

            if(transb=='C' || transb=='c')
                transactionB = CUBLAS_OP_C;

            cublasStatus_t status = cublasCgemm(handle, transactionA, transactionB, m, n, k, &alpha,
                        A, lda, B, ldb,
                        &beta, C, ldc);

            check_status(status);
        }

    static void gemm(char transa, char transb, int m, int n, int k, cuDoubleComplex alpha,
                         const cuDoubleComplex * const A, int lda, const cuDoubleComplex* const B, int ldb,
                         cuDoubleComplex beta, cuDoubleComplex * C, int ldc)
        {
            cublasOperation_t transactionA, transactionB;

            if(transa=='N' || transa=='n')
                transactionA = CUBLAS_OP_N;

            if(transa=='T' || transa=='t')
                transactionA = CUBLAS_OP_T;

            if(transa=='C' || transa=='c')
                transactionA = CUBLAS_OP_C;


            if(transb=='N' || transb=='n')
                transactionB = CUBLAS_OP_N;

            if(transb=='T' || transb=='t')
                transactionB = CUBLAS_OP_T;

            if(transb=='C' || transb=='c')
                transactionB = CUBLAS_OP_C;

            cublasStatus_t status = cublasZgemm(handle, transactionA, transactionB, m, n, k, &alpha,
                        A, lda, B, ldb,
                        &beta, C, ldc);

            check_status(status);
        }


    // @sect4{Funktion: dgmm}
    //! This function performs the matrix-matrix multiplication
    //! C = A × diag(x) if  mode == CUBLAS_SIDE_RIGHT
    //! C = diag(x) × A if  mode == CUBLAS_SIDE_LEFT
    //! @param m : number of rows of matrix A and C.
    //! @param n : number of columns of matrix A and C.
    //! @param A : <type> array of dimensions lda x n with lda>=max(1,m)
    //! @param lda : leading dimension of two-dimensional array used to store the matrix A.
    //! @param x : one-dimensional <type> array of size | i n c | × m if mode ==
    //!            CUBLAS_SIDE_LEFT and | i n c | × n if mode == CUBLAS_SIDE_RIGHT
    //! @param incx : stride of one-dimensional array x.
    //! @param C : <type> array of dimensions ldc x n with ldc>=max(1,m).
    //! @param ldc : leading dimension of a two-dimensional array used to store the matrix C.


    static void dgmm(cublasSideMode_t mode,
                     int m, int n,
                     const float *A, int lda,
                     const float *x, int incx,
                     float *C, int ldc)
    {

        cublasStatus_t status = cublasSdgmm( handle, mode,
                                   m,  n,
                                   A,  lda,
                                   x, incx,
                                   C, ldc);

        check_status(status);
    }

    static void dgmm(cublasSideMode_t mode,
                     int m, int n,
                     const double *A, int lda,
                     const double *x, int incx,
                     double *C, int ldc)
    {

        cublasStatus_t status = cublasDdgmm( handle, mode,
                                   m,  n,
                                   A,  lda,
                                   x, incx,
                                   C, ldc);

        check_status(status);
    }


    static void dgmm(cublasSideMode_t mode,
                     int m, int n,
                     const float2 *A, int lda,
                     const float2 *x, int incx,
                     float2 *C, int ldc)
    {

        cublasStatus_t status = cublasCdgmm( handle, mode,
                                   m,  n,
                                   A,  lda,
                                   x, incx,
                                   C, ldc);

        check_status(status);
    }

    static void dgmm(cublasSideMode_t mode,
                     int m, int n,
                     const double2 *A, int lda,
                     const double2 *x, int incx,
                     double2 *C, int ldc)
    {

        cublasStatus_t status = cublasZdgmm( handle, mode,
                                   m,  n,
                                   A,  lda,
                                   x, incx,
                                   C, ldc);

        check_status(status);
    }


    // @sect4{Funktion: nrm2}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.

    static float
    nrm2(int n, const float *x, int incx)
    {
        float result;
        cublasStatus_t status = cublasSnrm2 (handle, n, x, incx, &result);
        check_status(status);

        return result;
    }

    static double
    nrm2(int n, const double *x, int incx)
    {
        double result;
        cublasStatus_t status = cublasDnrm2 (handle, n, x, incx, &result);
        check_status(status);

        return result;
    }

    static float
    nrm2(int n, const cuComplex *x, int incx)
    {
        cuComplex result = dotc(n, x, incx, x, incx);
        result.x = sqrt(result.x);
        result.y = 0.;

        return result.x;

    }



    static double
    nrm2(int n, const cuDoubleComplex *x, int incx)
    {
        cuDoubleComplex result = dotc(n, x, incx, x, incx);
        result.x = sqrt(result.x);
        result.y = 0.;

        return result.x;

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
        float result;
        cublasStatus_t status = cublasSdot(handle, n, x, incx, y, incy, &result);
        check_status(status);

        return result;
    }

    static double dot(int n, const double *x,
                      int incx, const double *y, int incy)
    {
        double result;
        cublasStatus_t status = cublasDdot(handle, n, x, incx, y, incy, &result);
        check_status(status);

        return result;
    }

    static cuComplex dot(int n, const cuComplex *x,
                      int incx, const cuComplex *y, int incy)
    {
        return dotc(n, x, incx, y, incy);
    }

    static cuDoubleComplex dot(int n, const cuDoubleComplex *x,
                      int incx, const cuDoubleComplex *y, int incy)
    {
        return dotc(n, x, incx, y, incy);
    }

    // @sect4{Funktion: dotu}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor.
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.
    //! Funktioniert auch fuer cufftComplex, da dieses per typedef
    //! cuComlex cufftComplex definiert ist.
    static cuComplex dotu(int n, const cuComplex *x,
                      int incx, const cuComplex *y, int incy)
    {
        cuComplex result;
        cublasStatus_t status = cublasCdotu(handle, n, x, incx, y, incy, &result);
        check_status(status);

        std::cout << __FUNCTION__ << " : " << " c*c^* = " << result.x << ", " << result.y << std::endl;

        return result;
    }


    static cuDoubleComplex dotu(int n, const cuDoubleComplex *x,
                      int incx, const cuDoubleComplex *y, int incy)
    {
        cuDoubleComplex result;
        cublasStatus_t status = cublasZdotu(handle, n, x, incx, y, incy, &result);
        check_status(status);

        return result;
    }

    // @sect4{Funktion: dotc}
    //!
    //! @param n : Anzahl Elemente.
    //! @param x : Quellvektor x.
    //! @param incx : Speicher Abstand zwischen Elemente in Vector x.
    //! @param y : Zielvektor.
    //! @param incy: Speicher Abstand zwischen Elemente in Vector y.
    //! Funktioniert auch fuer cufftComplex, da dieses per typedef
    //! cuComlex cufftComplex definiert ist.
    static cuComplex dotc(int n, const cuComplex *x,
                      int incx, const cuComplex *y, int incy)
    {
        cuComplex result;
        cublasStatus_t status = cublasCdotc(handle, n, x, incx, y, incy, &result);
        check_status(status);

        //!std::cout << __FUNCTION__ << " : " << " c*c^* = " << result.x << ", " << result.y << std::endl;

        return result;
    }


    static cuDoubleComplex dotc(int n, const cuDoubleComplex *x,
                      int incx, const cuDoubleComplex *y, int incy)
    {
        cuDoubleComplex result;
        cublasStatus_t status  = cublasZdotc(handle, n, x, incx, y, incy, &result);
        check_status(status);

        return result;
    }


    // @sect4{Funktion: trsm}
    //!
    static inline void trsm(char side, char uplo, char transa, char diag,  int m, int n,
                     float alpha,
                     const float * A, int lda, float * B, int ldb)
    {

        cublasSideMode_t sideMode;
        cublasFillMode_t uploMode;
        cublasOperation_t transaction;
        cublasDiagType_t diagMode;

        if(side=='L' || side=='l')
            sideMode = CUBLAS_SIDE_LEFT;

        if(side=='R' || side=='r')
            sideMode = CUBLAS_SIDE_RIGHT;


        if(uplo=='U' || uplo=='u')
            uploMode = CUBLAS_FILL_MODE_UPPER;

        if(uplo=='L' || uplo=='l')
            uploMode = CUBLAS_FILL_MODE_LOWER;


        if(transa=='N' || transa=='n')
            transaction = CUBLAS_OP_N;

        if(transa=='T' || transa=='t')
            transaction = CUBLAS_OP_T;

        if(transa=='C' || transa=='c')
            transaction = CUBLAS_OP_C;


        if(diag=='U' || diag=='u')
            diagMode = CUBLAS_DIAG_UNIT;

        if(diag=='N' || diag=='n')
            diagMode = CUBLAS_DIAG_NON_UNIT;

        cublasStatus_t status = cublasStrsm(handle, sideMode, uploMode, transaction, diagMode, m, n, &alpha,
                    A, lda, B, ldb);

        check_status(status);
    }

    static inline void trsm(char side, char uplo, char transa, char diag,  int m, int n,
                     double alpha,
                     const double * A, int lda, double * B, int ldb)
    {
        cublasSideMode_t sideMode;
        cublasFillMode_t uploMode;
        cublasOperation_t transaction;
        cublasDiagType_t diagMode;

        if(side=='L' || side=='l')
            sideMode = CUBLAS_SIDE_LEFT;

        if(side=='R' || side=='r')
            sideMode = CUBLAS_SIDE_RIGHT;


        if(uplo=='U' || uplo=='u')
            uploMode = CUBLAS_FILL_MODE_UPPER;

        if(uplo=='L' || uplo=='l')
            uploMode = CUBLAS_FILL_MODE_LOWER;


        if(transa=='N' || transa=='n')
            transaction = CUBLAS_OP_N;

        if(transa=='T' || transa=='t')
            transaction = CUBLAS_OP_T;

        if(transa=='C' || transa=='c')
            transaction = CUBLAS_OP_C;


        if(diag=='U' || diag=='u')
            diagMode = CUBLAS_DIAG_UNIT;

        if(diag=='N' || diag=='n')
            diagMode = CUBLAS_DIAG_NON_UNIT;

        cublasStatus_t status = cublasDtrsm(handle, sideMode, uploMode, transaction, diagMode, m, n, &alpha,
                    A, lda, B, ldb);

        check_status(status);
    }

}; //! struct CW END




#endif //! CUBLAS_WRAPPER_HH
