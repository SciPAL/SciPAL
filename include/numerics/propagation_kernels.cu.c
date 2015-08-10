// Header-files of CUDA utility-library
//Header for referencing CUDA based part
//! \file
#include "propagation_kernels_wrapper.cu.h"
#include <iostream>
#include <base/ParallelArch.h>
#include <lac/Shape.h>
#include <typeinfo>

#define TILE_DIM (32)
#define BLOCK_ROWS (8)

using namespace SciPAL;

//////////////////////////////////////////////////////////////////////////////////////////

// @sect4{Kernel: Host + Devce: __transpose2_element}
//is for aike
//The template allows to test the algorithm with T=float2 or double2 numbers.
template <typename T, ParallelArch arch>
__host__ __device__
void
__transpose2_element(T *d_devPtr, int width, int height, int index)
{

    typedef typename PrecisionTraits<T, arch>::NumberType NumberType;
    typedef typename SciPAL::CudaComplex<NumberType> Complex;

    int x;
    int y;
    x = index % (width);
    y = index / (width);
    if( x < y)
    {
        Complex val(d_devPtr[index]);
        int new_index = x * width + y;
        Complex tmp(d_devPtr[new_index]);
        d_devPtr[new_index] = toNumberType2(val);
        d_devPtr[index] = toNumberType2(tmp);
        //        printf("i am working: x=%d, y=%d, index=%d, new_index=%d, val=(%f,%f), tmp=(%f,%f) "
        //               "my result: d_devPtr[new_index] = (%f, %f), d_devPtr[index]=(%f, %f)\n",
        //               x ,y, index, new_index, val.real(), val.imag(), tmp.real(), tmp.imag(),
        //               d_devPtr[new_index].x, d_devPtr[new_index].y, d_devPtr[index].x, d_devPtr[index].y);

    }
    else
    {
        //     printf("i am lazy: x=%d, y=%d\n", x, y);
    }

}

// No bank-conflict transpose
// Same as transposeCoalesced except the first tile dimension is padded
// to avoid shared memory bank conflicts.



//const int TILE_DIM =32;
//const int BLOCK_ROWS = 8;


template<typename T>
__global__  void
transposeNoBankConflicts_inplace
(T *data_ptr)
{
    typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;
    __shared__ T tile[TILE_DIM][TILE_DIM+1];

    if( blockIdx.x <= blockIdx.y)
    {
        int x = blockIdx.x * TILE_DIM + threadIdx.x;
        int y = blockIdx.y * TILE_DIM + threadIdx.y;
        int width = gridDim.x * TILE_DIM;

        //store date from source patch in shared mem
        for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
            tile[threadIdx.y+j][threadIdx.x] = data_ptr[(y+j)*width + x];

        __syncthreads();

        x = blockIdx.y * TILE_DIM + threadIdx.x;  // transpose block offset
        y = blockIdx.x * TILE_DIM + threadIdx.y;

        if (typeid(T) == typeid(double2))
        {
        unsigned long long int tmp2;
        ulonglong2 * pointer =
                reinterpret_cast<ulonglong2*>(&data_ptr[0]);

        ulonglong2 * tile_pointer =
                reinterpret_cast<ulonglong2*>(&tile[0][0]);

        // write contents of shared mem to target tile, while doing that
        // we exchange contents of shared mem with values from target tile
        for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
        {

            tmp2 = atomicExch(&(pointer[(y+j)*width + x].x),
                    tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].x);

            tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].x = (tmp2);

            tmp2 = atomicExch(&(pointer[(y+j)*width + x].y),
                    tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].y);

            tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].y = (tmp2);

        }
        }
        else if (sizeof(T) == 8 /*64 bit data types*/ )
        {
        unsigned long long int tmp2;
        ulonglong1 * pointer =
                reinterpret_cast<ulonglong1*>(&data_ptr[0]);

        ulonglong1 * tile_pointer =
                reinterpret_cast<ulonglong1*>(&tile[0][0]);

        // write contents of shared mem to target tile, while doing that
        // we exchange contents of shared mem with values from target tile
        for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
        {

            tmp2 = atomicExch(&(pointer[(y+j)*width + x].x),
                    tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].x);

            tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j].x = (tmp2);

        }
        }
        else if (sizeof(T) == 4 /*32 bit data types*/ )
        {
        unsigned long long int tmp2;
        int * pointer =
                reinterpret_cast<int*>(&data_ptr[0]);

        int * tile_pointer =
                reinterpret_cast<int*>(&tile[0][0]);

        // write contents of shared mem to target tile, while doing that
        // we exchange contents of shared mem with values from target tile
        for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
        {

            tmp2 = atomicExch(&(pointer[(y+j)*width + x]),
                    tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j]);

            tile_pointer[threadIdx.x * (TILE_DIM+1) + threadIdx.y + j] = (tmp2);

        }
        }


        __syncthreads();

        //contents of shared mem already transposed, so we can go back to “normal”
        //x, y coordinates
        x = blockIdx.x * TILE_DIM + threadIdx.x;
        y = blockIdx.y * TILE_DIM + threadIdx.y;

        //write contents of shared mem to source tile.
        for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
            data_ptr[(y+j)*width + x] = tile[threadIdx.y+j][threadIdx.x];

        __syncthreads();

    }
}

// @sect4{Kernel: __transpose2}
//
//This kernel generates rhe quadratic phase factors for Fresnel propagation.

template <typename T>
__global__
void
__transpose2(T *d_devPtr, int width,
             int height, int size)
{
    //Calculate the thread ID. The thread ID determines which pixel is calculated.
    int x = blockDim.x*blockIdx.x+threadIdx.x;

    //Prevents kernel to calculate something outside the image vector.
    if(x<size)
        __transpose2_element<T,gpu_cuda>(d_devPtr, width, height, x);

}

//Wrapper functions
template<typename T>
template< typename dummy>
template <typename U>
inline void
PropagationKernelsImpl<T>::Impl<gpu_cuda, dummy>
::transpose2(U *d_devPtr, int width,
             int height, int size )
{
    //#if __CUDA_ARCH__ < 200
    //    int threads_per_block = 512;
    //#else
    //    int threads_per_block = 1024;
    //#endif
    //    int threads_per_block = 512;
    //    int blocks = (size + threads_per_block - 1) / threads_per_block;

    dim3 dimGrid(width/TILE_DIM, height/TILE_DIM, 1);
    dim3 dimBlock(TILE_DIM, BLOCK_ROWS, 1);

    transposeNoBankConflicts_inplace<T><<<dimGrid, dimBlock>>>(d_devPtr);

    //    __transpose2<T><<<blocks,threads_per_block>>>(d_devPtr, width, height, size);
    cudaThreadSynchronize();

}
//CPU specialization:

template<typename T>
template< typename dummy>
template <typename U>
inline void
PropagationKernelsImpl<T>::Impl<cpu, dummy>::transpose2(U *d_devPtr, int width,
                                                        int height, int size )
{
#pragma omp parallel for
    for(int i = 0;i < size;i++)
    {__transpose2_element<T, cpu>(d_devPtr, width, height, i); }

}

/////////////////////////////////////////////////////////////////////////////////

// @sect4{Kernel: __merge_complex}
//
//This kernel generates rhe quadratic phase factors for Fresnel propagation.

template <typename T>
__global__
void
__merge_cplx(typename PrecisionTraits<T, gpu_cuda>::NumberType *re,
             typename PrecisionTraits<T, gpu_cuda>::NumberType *im,
             T *d_devPtr)
{

    {
        typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;
        __shared__ NumberType retile[TILE_DIM][TILE_DIM+1];
        __shared__ NumberType imtile[TILE_DIM][TILE_DIM+1];

        //        if( blockIdx.x <= blockIdx.y)
        {
            int x = blockIdx.x * TILE_DIM + threadIdx.x;
            int y = blockIdx.y * TILE_DIM + threadIdx.y;
            int width = gridDim.x * TILE_DIM;

            //store data from source patch in shared mem
            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                retile[threadIdx.y+j][threadIdx.x] = re[(y+j)*width + x];

            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                imtile[threadIdx.y+j][threadIdx.x] = im[(y+j)*width + x];

            __syncthreads();


            //write contents of shared mem to target.
            typedef typename SciPAL::CudaComplex<NumberType> Complex;
            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                d_devPtr[(y+j)*width + x] =
                        toNumberType2(Complex(retile[threadIdx.y+j][threadIdx.x],
                        imtile[threadIdx.y+j][threadIdx.x]));

            __syncthreads();

        }
    }


}

//Wrapper functions
template<typename T>
template< typename dummy>
inline void
PropagationKernelsImpl<T>::Impl<gpu_cuda, dummy>::merge_cplx(rShape &re, rShape &im,
                                                             cShape &d_devPtr,
                                                             int width, int height, int size )
{
    //#if __CUDA_ARCH__ < 200
    //    int threads_per_block = 512;
    //#else
    //    int threads_per_block = 1024;
    //#endif
    //    int blocks = (size + threads_per_block - 1) / threads_per_block;

    dim3 dimGrid(width/TILE_DIM, height/TILE_DIM, 1);
    dim3 dimBlock(TILE_DIM, BLOCK_ROWS, 1);

    __merge_cplx<T><<<dimGrid, dimBlock>>>(re.data_ptr, im.data_ptr, d_devPtr.data_ptr);
    cudaThreadSynchronize();

}
//CPU specialization:

template<typename T>
template< typename dummy>
inline void
PropagationKernelsImpl<T>::Impl<cpu, dummy>::merge_cplx(rShape &re, rShape &im,
                                                        cShape &d_devPtr,
                                                        int width, int height, int size )
{
    std::cerr<< "not implemented for cpu" << std::endl;

}
/////////////////////////////////////////////////////////////////////////////////

// @sect4{Kernel: __separate_complex}
//
//This kernel generates rhe quadratic phase factors for Fresnel propagation.

template <typename T>
__global__
void
__separate_cplx(typename PrecisionTraits<T, gpu_cuda>::NumberType *re,
                typename PrecisionTraits<T, gpu_cuda>::NumberType *im,
                T *d_devPtr)
{

    {
        typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;
        typedef typename PrecisionTraits<T, gpu_cuda>::ComplexType ComplexType;
        typedef typename SciPAL::CudaComplex<NumberType> Complex;

        __shared__ ComplexType tile[TILE_DIM][TILE_DIM+1];

        {
            int x = blockIdx.x * TILE_DIM + threadIdx.x;
            int y = blockIdx.y * TILE_DIM + threadIdx.y;
            int width = gridDim.x * TILE_DIM;

            //read from complex vector
            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                tile[threadIdx.y+j][threadIdx.x] = (d_devPtr[(y+j)*width + x]);

            __syncthreads();


            //store data from source patch in shared mem
            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                re[(y+j)*width + x] = tile[threadIdx.y+j][threadIdx.x].x;

            for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
                im[(y+j)*width + x] = tile[threadIdx.y+j][threadIdx.x].y;

            __syncthreads();

        }
    }
}

//Wrapper functions
template<typename T>
template< typename dummy>
inline void
PropagationKernelsImpl<T>::Impl<gpu_cuda, dummy>::separate_cplx(rShape &re, rShape &im,
                                                                cShape &d_devPtr,
                                                                int width, int height, int size )
{
    //#if __CUDA_ARCH__ < 200
    //    int threads_per_block = 512;
    //#else
    //    int threads_per_block = 1024;
    //#endif
    //    int blocks = (size + threads_per_block - 1) / threads_per_block;

    dim3 dimGrid(width/TILE_DIM, height/TILE_DIM, 1);
    dim3 dimBlock(TILE_DIM, BLOCK_ROWS, 1);

    __separate_cplx<T><<<dimGrid, dimBlock>>>(re.data_ptr, im.data_ptr, d_devPtr.data_ptr);
    cudaThreadSynchronize();

}
//CPU specialization:

template<typename T>
template< typename dummy>
inline void
PropagationKernelsImpl<T>::Impl<cpu, dummy>::separate_cplx(rShape &re, rShape &im,
                                                           cShape &d_devPtr,
                                                           int width, int height, int size )
{
    std::cerr<< "not implemented for cpu" << std::endl;

}


//////////////////////////////////////////////////////////////////////////////////////////
//This has to be at the end of file, because all functions have to be declared and their bodies defined before the
//class can be explictly instantiated by the compiler.
//The template parameter is in this case float2/double2, because the Kernels work on elements
//of these types.
template class PropagationKernels<float2,gpu_cuda>;
template class PropagationKernels<double2,gpu_cuda>;

template class PropagationKernelsImpl<float2>::Impl<gpu_cuda,float2>;
template class PropagationKernelsImpl<double2>::Impl<gpu_cuda,double2>;
template void PropagationKernelsImpl<double2>::Impl<(ParallelArch)1, double2>::transpose2<SciPAL::CudaComplex<double> >(SciPAL::CudaComplex<double>*, int, int, int);

template class PropagationKernels<float2,cpu>;
template class PropagationKernels<double2,cpu>;

template class PropagationKernelsImpl<float2>::Impl<cpu,float2>;
template class PropagationKernelsImpl<double2>::Impl<cpu,double2>;

//template void PropagationKernelsImpl<double2>::Impl<(ParallelArch)0, double2>::transpose<pcomplex>(SciPAL::Shape<pcomplex>&, int, int, int);
//template void PropagationKernelsImpl<float2>::Impl<cpu, float2>::amplitude_adaption_r<MagProj>(float2 *, float2 const*, int );
//template void PropagationKernelsImpl<double2>::Impl<cpu, double2>::amplitude_adaption_r<MagProj>(double2*,double2 const*, int );

