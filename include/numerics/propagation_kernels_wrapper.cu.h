#ifndef PROPAGATION_KERNELS_CU_H
#define PROPAGATION_KERNELS_CU_H
//CUDA specific include
//! \file
#include <omp.h>
#include <base/CudaComplex.h>
#include <base/PrecisionTraits.h>
#include <lac/Shape.h>

using namespace SciPAL;

// @sect3{Struct: PropagationDataStruct}
template<typename NumberType> //, ParallelArch arch>
struct
PropagationDataStruct{

    PropagationDataStruct( NumberType _delta_z,
                           NumberType _k, //2π/λ
                           NumberType _pixel_size,
                           int _width, int _height):
        delta_z(_delta_z), k(_k), k2(k *k), pixel_size(_pixel_size), width(_width), height(_height),
        delta_x2(4*M_PI*M_PI / (pixel_size * pixel_size* width * width)),
        delta_y2(4*M_PI*M_PI / (pixel_size * pixel_size* height * height)),
        x_0(_width/2), y_0(height/2), size(width*height)
    {
//        if(arch == gpu_cuda)
        {
//        cudaMemcpyToSymbol(ddelta_x2, &delta_x2, sizeof(double));
//        cudaMemcpyToSymbol(ddelta_y2, &delta_y2, sizeof(double));
//        cudaMemcpyToSymbol(dwidth, &width, sizeof(int));
//        cudaMemcpyToSymbol(dheight, &height, sizeof(int));
//        cudaMemcpyToSymbol(dx0, &x_0, sizeof(int));
//        cudaMemcpyToSymbol(dy0, &y_0, sizeof(int));
//        cudaMemcpyToSymbol(ddelta_z, &delta_z, sizeof(double));
        }
//        else
//        {
//            ddelta_x2 = delta_x2;
//        }

    }

    NumberType delta_z;
    NumberType k;
    NumberType k2; //squared 2π/λ
    NumberType pixel_size;
    int width; int height;
    NumberType delta_x2;
    NumberType delta_y2;
    int x_0;
    int y_0;
    int size;
};

// @sect3{Class: PropagationKernelsImpl}
template<typename T>
struct PropagationKernelsImpl {

// @sect3{Class: Impl}
template <ParallelArch arch, typename dummy > class Impl;

template <typename dummy> class Impl<gpu_cuda, dummy> {
public:

//        typedef typename PrecisionTraits<T, gpu_cuda>::ComplexType NumberType2;
        typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;
        typedef typename SciPAL::CudaComplex<NumberType> NumberType2;
        typedef typename SciPAL::Shape<NumberType2> cShape;
        typedef typename SciPAL::Shape<NumberType> rShape;

        void generate_k(cShape &d_devPtr, NumberType delta_z,
                        NumberType k, NumberType pixel_size, int width,
                        int height, int size);

        void fourier_factor(rShape &sinfac,rShape &cosfac,
                            NumberType delta_z,NumberType k, NumberType pixel_size,
                            int width, int height, int size);

        void transpose(NumberType2* d_devPtr, int width, int height, int size );

        template <typename U>
        void transpose2(U *d_devPtr, int width, int height, int size );

        void separate_cplx(rShape &re, rShape &im,
                           cShape &d_devPtr,
                           int width, int height, int size );

        void merge_cplx(rShape &re, rShape &im,
                           cShape &d_devPtr,
                           int width, int height, int size );


private:
};

template <typename dummy> class Impl<cpu, dummy> {
    public:

        typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;
        typedef typename SciPAL::CudaComplex<NumberType> NumberType2;
        typedef typename SciPAL::Shape<NumberType2> cShape;
        typedef typename SciPAL::Shape<NumberType> rShape;

        void generate_k(cShape &d_devPtr, NumberType delta_z,
                        NumberType k, NumberType pixel_size, int width,
                        int height, int size);

        void fourier_factor(rShape &sinfac, rShape &cosfac,
                            NumberType delta_z,NumberType k, NumberType pixel_size,
                            int width, int height, int size);

        void transpose(NumberType2* d_devPtr, int width, int height, int size );

        template <typename U>
        void transpose2(U *d_devPtr, int width, int height, int size );

        void separate_cplx(rShape &re, rShape &im,
                           cShape &d_devPtr,
                           int width, int height, int size );

        void merge_cplx(rShape &re, rShape &im,
                           cShape &d_devPtr,
                           int width, int height, int size );

private:
    };

};

// @sect3{Class: PropagationKernels}
// Primary template.
template<typename T, ParallelArch arch>
struct PropagationKernels
        : PropagationKernelsImpl<T>::template Impl<arch, T>
{};

template<typename T>
struct PropagationKernels<T, gpu_cuda>
        : PropagationKernelsImpl<T>::template Impl<gpu_cuda, T>
{
    PropagationKernels(int /*num_omp_threads*/)
    {}
};

template<typename T>
struct PropagationKernels<T,cpu>
        : PropagationKernelsImpl<T>::template Impl<cpu, T>
{
    PropagationKernels(int num_omp_threads)
    {
        omp_set_num_threads(num_omp_threads);
    }
};

#endif // PROPAGATION_KERNELS_CU_H
