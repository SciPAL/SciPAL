// To outfox QTCreator's syntax highlighting and especially
// nvcc we put all cuda-related code into files with names
// ending on .cu.c and include them here.
// Then, in the project file we only have to take of one source
// file. This reduces the amount of maintenance.
#include "cuda_utils.cu.c"

#include <base/PrecisionTraits.h>

#include <thrust/functional.h>
// #include <thrust/device_ptr.h>
// #include <thrust/execution_policy.h>

#include <cuComplex.h>

//#ifndef CUDAHELPER_H
//#define CUDAHELPER_H

//#include<stdlib.h>

////@sect4{CUDA helper functions}
////CUDA helper functions  for error handling

////@sect5{checkCudaErrors}
////@brief print the proper CUDA error strings if a CUDA host call returns an error
//#define checkCudaErrors(err)    __checkCudaErrors (err, __FILE__, __LINE__)

//inline void __checkCudaErrors( cudaError err, const char *file, const int line ) {
//    if( cudaSuccess != err) {
//        printf("%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString(err));
//        std::abort();
//    }
//}

////@sect5{getLastCudaError}
////@brief print the proper error string of a previous kernel call when calling cudaGetLastError
//#define getLastCudaError(msg)   __getLastCudaError (msg, __FILE__, __LINE__)

//inline void __getLastCudaError( const char *errorMessage, const char *file, const int line ) {
//    cudaError_t err = cudaGetLastError();
//    if( cudaSuccess != err) {
//        printf("%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n", file, line, errorMessage, (int)err, cudaGetErrorString(err));
//        std::abort();
//    }
//}
//#endif // CUDAHELPER_H




//@sect5{Struct: __abs}
//@brief Unary functions to calculate the absolute value for an element
template <typename T>
struct __abs : public thrust::unary_function<T,T> {
    // Ensure that we won't compile any un-specialized types
    __device__
    T operator()(T x) {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};


template <>
struct __abs <float> : public thrust::unary_function<float,float> {
    __device__
    float operator()(float x) {
        return fabsf(x);
    }
};

// @sect5{Function: element_norm_product}
//@brief Elementwise multiplication of two arrays with complex values and division by the number of elements in each array.
//@param arr1 first input array, output will be calculated inplace here
//@param arr2 second input array
//@param width width of the array
//@param height height of the array

typedef PrecisionTraits<double, gpu_cuda> CType;

//Double specialization
//template <>
//struct __element_norm_product
//        // <double>
//{
//    cuDoubleComplex normalization;

//    __element_norm_product(double _n) {
//        normalization = make_cuDoubleComplex(1.0/_n, 0);
//    }

//    __host__ __device__
//    cuDoubleComplex operator()(const cuDoubleComplex& x, const cuDoubleComplex& y) const {
//        return cuCmul(cuCmul(x, y), normalization);
//    }
//};

// template<typename T>
void /*step35::Kernels<T>::*/
element_norm_product(CType* arr1,
        CType* arr2,
        int width, int height, int depth)
{
    int size=depth*width*(height/2+1);
    double normalization=width*height*depth;
   //  __element_norm_product op(normalization);

   // thrust::device_ptr<CType> ptr1 = thrust::device_pointer_cast(arr1);
   // thrust::device_ptr<CType> ptr2 = thrust::device_pointer_cast(arr2);
   // thrust::transform(thrust::device, ptr1, ptr1 + size, ptr2, ptr1, op);
    cudaDeviceSynchronize();
  //  getLastCudaError("element_norm_product<<<>>> execution failed\n");
}

#include "cuda_kernel_step-28.cu.c"
