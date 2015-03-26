// Header-files of CUDA utility-library
//Header for referencing CUDA based part
#include <lac/scipal_kernels_wrapper.cu.h>

//UnaryExpressions

template<typename T, typename X, typename op>
__host__ __device__
void __apply_element(T *d_dst,
                     const SciPAL::DevUnaryExpr<X, op> & Ax,
                     int i)
{
    // typedef typename PrecisionTraits<T, arch>::NumberType NumberType;
    // typedef typename Complex<NumberType, arch>::ComplexType Complex;
    //hier vielleicht noch auspacken der Ax innereien um doppeltes lesen auf glob mem
    //zu verhindern
    //d_dst[i] = Ax(i);
    typedef SciPAL:: ExprTree<T> e;
    d_dst[i] = e::eval(Ax, i);
}

// @sect4{Kernel: __apply}
//
//This kernel

template <typename T, typename X, typename op>
__global__
void
__apply(T *d_dst, const SciPAL::DevUnaryExpr<X, op> Ax, int size)
{
    //Calculate the thread ID. The thread ID determines which pixel is calculated.
    int x = blockDim.x*blockIdx.x+threadIdx.x;

    //Prevents kernel to calculate something outside the image vector.
    if(x<size)
     __apply_element<T, X, op>(d_dst, Ax, x);

}

//Wrapper functions
template<typename T>
template <typename X, typename op>
void
ImplCUDA<T>::apply(SciPAL::ShapeData<T> &d_dst,
                  const SciPAL::DevUnaryExpr<X, op> & Ax)
{
#if __CUDA_ARCH__ < 200
    int threads_per_block = 512;
#else
    int threads_per_block = 1024;
#endif
    int size = d_dst.n_rows * d_dst.n_cols;
    int blocks = (size + threads_per_block - 1) / threads_per_block;

    __apply<T, X><<<blocks, threads_per_block>>>(d_dst.data_ptr, Ax, size);
    cudaDeviceSynchronize();

}

//CPU specialization:

template <typename T>
template <typename X, typename op>
void
ImplOpenMP<T>::apply(SciPAL::ShapeData<T> &d_dst,
                  const SciPAL::DevUnaryExpr<X, op> &Ax)
{
    #pragma omp parallel for
    for(int i = 0; i < d_dst.n_rows * d_dst.n_cols; i++)
     __apply_element<T, X, op>(d_dst.data_ptr, Ax, i);
}
/////////////////////////////////////////
//Binary Function
template <typename T, typename L, typename op, typename R>
__device__ __host__ __forceinline__
void __apply_element(T *d_dst,
                     const SciPAL::DevBinaryExpr<L, op, R> & Ax,
                     int i)
{

    //hier vielleicht noch auspacken der Ax innereien um doppeltes lesen auf glob mem
    //zu verhindern
    // By construction expressions provide element-wise access by operator[].
    T tmp;
    tmp = Ax[i];
    d_dst[i] = tmp;
}

// @sect4{Kernel: __apply}
//
//This kernel

template <typename T, typename L, typename op, typename R>
__global__
void
__apply(T *d_dst,
        const ::SciPAL::DevBinaryExpr<L, op, R> Ax, int size)
{
    //Calculate the thread ID. The thread ID determines which pixel is calculated.
    int x = blockDim.x*blockIdx.x+threadIdx.x;

    //Prevents kernel to calculate something outside the image vector.
//    if(x<size)
     __apply_element<T, L, op, R>(d_dst, Ax, x);

}

//Wrapper functions
template< typename T>
template <typename L, typename op, typename R>
void
ImplCUDA<T>::apply(SciPAL::ShapeData<T> & d_dst,
                  const SciPAL::DevBinaryExpr<L, op, R> & Ax)
{
#if __CUDA_ARCH__ < 200
    int threads_per_block = 512;
#else
    int threads_per_block = 1024;
#endif

    int size = d_dst.n_rows * d_dst.n_cols;
    int blocks = (size + threads_per_block - 1) / threads_per_block;

    __apply<T, L, op, R><<<blocks, threads_per_block>>>
                        (d_dst.data_ptr, Ax, size);
    cudaDeviceSynchronize();
}

//CPU specialization:

template<typename T>
template <typename L, typename op, typename R>
void
ImplOpenMP<T>::apply(SciPAL::ShapeData<T> & d_dst,
                  const SciPAL::DevBinaryExpr<L, op, R> &Ax)
{
    #pragma omp parallel for
    for(int i = 0; i < d_dst.n_rows * d_dst.n_cols; i++)
     __apply_element<T, L, op, R>(d_dst.data_ptr, Ax, i);

}



//////////////////////////////////////////////////////////////////////////////////////////

namespace bw_types {
template <typename T, typename BW> class Vector;
}
struct blas;
struct cublas;

#include <lac/instantiations.h>
