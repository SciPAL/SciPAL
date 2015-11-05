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


#ifndef CUDA_FFT_H
#define CUDA_FFT_H

//including the headers in order, forces an overwrite of the standard complex format by the fftw format
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
//STL
#include <complex>
#include <iostream>
//FFTs: FFTW and cuFFT
#include <fftw3.h>
#include <cufft.h>
//deal.II
#include <deal.II/base/parameter_handler.h>
//SciPAL
#include <base/VTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/Shape.h>
#include <numerics/propagation_kernels_wrapper.cu.h>
//Qt Includes
#include <QDir>
#include <QString>

enum TransformDirection {FORWARD, BACKWARD};
namespace SciPAL {


// @sect3{struct: TransformType}
// The struct TransformType chosses to type of tranformation e.g. complex to complex.
// The type of source/target datatype is derived from this Argument. Details is FFTTraits.

template<typename T> struct TransformType;

template<>
struct TransformType<float>
{
    static const cufftType_t FFTType_C2C = CUFFT_C2C;
    static const cufftType_t FFTType_C2R = CUFFT_C2R;
    static const cufftType_t FFTType_R2C = CUFFT_R2C;
};

template<>
struct TransformType<double>
{
    static const cufftType_t FFTType_C2C = CUFFT_Z2Z;
    static const cufftType_t FFTType_C2R = CUFFT_Z2D;
    static const cufftType_t FFTType_R2C = CUFFT_D2Z;
};


//@sect3{struct: FFTTraits}
//This struct holds information of the source/destination datatype for different
//combinations of computing hardware and TransformType.
//Additionally every instantiation holds the appropriate execution function which is
//encapsulated in the exec-function.
template<cufftType_t cufft_type, ParallelArch arch> struct FFTTraits;

//FFTTraits specializations for CUDAFFT.
template<>
struct FFTTraits<CUFFT_C2C, gpu_cuda>
{
    typedef CudaComplex<float> T_dst;
    typedef CudaComplex<float> T_src;

    static cufftResult exec(cufftHandle &plan, cufftComplex *dst, cufftComplex *src, int direction)
    {
        return cufftExecC2C(plan, src, dst, direction);
    }
};

template<>
struct FFTTraits<CUFFT_R2C, gpu_cuda>
{
    typedef CudaComplex<float>    T_dst;
    typedef float T_src;

    static cufftResult exec(cufftHandle &plan, cufftComplex *dst, cufftReal *src, int /*direction*/)
    {
        return cufftExecR2C(plan, src, dst);
    }
};

template<>
struct FFTTraits<CUFFT_C2R, gpu_cuda>
{
    typedef float T_dst;
    typedef CudaComplex<float> T_src;

    static cufftResult exec(cufftHandle &plan, cufftReal *dst, cufftComplex *src, int /*direction*/)
    {
        return cufftExecC2R(plan, src, dst);
    }
};

template<>
struct FFTTraits<CUFFT_Z2Z, gpu_cuda>
{
    typedef CudaComplex<double> T_dst;
    typedef CudaComplex<double> T_src;

    static cufftResult exec(cufftHandle &plan, cufftDoubleComplex *dst, cufftDoubleComplex *src, int direction)
    {
        return cufftExecZ2Z(plan, src, dst, direction);
    }
};

template<>
struct FFTTraits<CUFFT_D2Z, gpu_cuda>
{
    typedef CudaComplex<double> T_dst;
    typedef double T_src;

    static cufftResult exec(cufftHandle &plan, cufftDoubleComplex *dst, cufftDoubleReal *src, int /*direction*/)
    {
        return cufftExecD2Z(plan, src, dst);
    }
};

template<>
struct FFTTraits<CUFFT_Z2D, gpu_cuda>
{
    typedef double T_dst;
    typedef CudaComplex<double> T_src;

    static cufftResult exec(cufftHandle &plan, cufftDoubleReal *dst, cufftDoubleComplex *src, int /*direction*/)
    {
        return cufftExecZ2D(plan, src, dst);
    }
};

//FFTTraits specializations for FFTW on the CPU.
template<>
struct FFTTraits<CUFFT_C2C, cpu>
{
    typedef CudaComplex<float> T_dst;
    typedef CudaComplex<float> T_src;
};

template<>
struct FFTTraits<CUFFT_R2C, cpu>
{
    typedef float     T_dst;
    typedef CudaComplex<float> T_src;
};

template<>
struct FFTTraits<CUFFT_C2R, cpu>
{
    typedef float T_dst;
    typedef CudaComplex<float> T_src;
};

template<>
struct FFTTraits<CUFFT_Z2Z, cpu>
{
    typedef CudaComplex<double> T_dst;
    typedef CudaComplex<double> T_src;
};

template<>
struct FFTTraits<CUFFT_D2Z, cpu>
{
    typedef CudaComplex<double> T_dst;
    typedef double T_src;
};

template<>
struct FFTTraits<CUFFT_Z2D, cpu>
{
    typedef double T_dst;
    typedef CudaComplex<double> T_src;
};

// @sect3{class: CUDAFFTbase}
// This base class holds some fundamental informations needed for the FFT
// e.g. dimension, datatypes.
template <typename T,int d, cufftType_t cufft_type, ParallelArch arch>
class CUDAFFTbase {

public:

    static const int         dim     = d;
    static const cufftType_t variant = cufft_type;

    typedef typename FFTTraits<cufft_type, arch>::T_dst T_dst;
    typedef typename FFTTraits<cufft_type, arch>::T_src T_src;


    typedef SciPAL::Shape<T_dst, typename archTraits<arch>::BlasType, matrix> FFTDataDST;
    typedef SciPAL::Shape<T_src, typename archTraits<arch>::BlasType, matrix> FFTDataSRC;



    void operator() (FFTDataDST &dst, FFTDataSRC &src);

    void check_folder()
    {
        QDir dir(QString("fftw_wisdom"));
        if(!dir.exists())
        {
            dir.mkpath(".");
            std::cout<<"creating folder"<<std::endl;
        }

    }

private:
};

// @sect3{struct: CUDAFFT}
// This struct inherits from CUDAFFTbase this is necessary in order to specialize all the
// needed template combinations correctly.
//Primary template of CUDAFFT.
template <typename T, int dim, cufftType_t cufft_type, ParallelArch arch>
struct CUDAFFT
        : SciPAL::CUDAFFTbase<T, dim, cufft_type, arch> {};

// These are the specializations for the CPU version. Since the FFTW commands differ for
// float and doble precision, we have to create a float and double specialization separately.
template <int dim, cufftType_t cufft_type>
struct CUDAFFT<float, dim, cufft_type, cpu>
        : SciPAL::CUDAFFTbase<float, dim, cufft_type, cpu>
{
    typedef typename CUDAFFTbase<float, dim, cufft_type, cpu>::FFTDataDST FFTDataDST;
    typedef typename CUDAFFTbase<float, dim, cufft_type, cpu>::FFTDataSRC FFTDataSRC;

    //creates fftw handle object CPU
    fftwf_plan plan_cpu_f_forward;
    fftwf_plan plan_cpu_f_backward;
    int width, height;
    FFTDataDST &dst;
    FFTDataSRC &src;
    bool initialized;
    TransformDirection old_direction;

    CUDAFFT(int _height, int _width, FFTDataDST &_dst, FFTDataSRC &_src): width(_width),
        height(_height), dst(_dst), src(_src), initialized(false)
    {
        //initializiation of multithreaded fftw
        dealii::ParameterHandler prm;
        prm.declare_entry("FFTW-Threads", "4",
                          dealii::Patterns::Integer(),
                          "sets the number of Threads, the FFTW uses");
        prm.read_input("fftw-param.prm");
        int num_threads = prm.get_integer("FFTW-Threads");
        int error = fftwf_init_threads();
#ifdef DEBUG
        std::cout<<"fftw says(0 bad, everything else ok):"<<error<<std::endl;
#endif
        fftwf_plan_with_nthreads(num_threads);
    }


    ~CUDAFFT()
    {
        fftwf_destroy_plan(this->plan_cpu_f_forward);
        fftwf_destroy_plan(this->plan_cpu_f_backward);

    }
    void reinit(bool inverse)
    {
        switch(dim)
        {
        case 1:
            // TODO: does not cover c2r r2c etc...
            this->plan_cpu_f_forward = fftwf_plan_dft_1d(width,
                                                 reinterpret_cast<fftwf_complex*>(src.data_ptr),
                                                 reinterpret_cast<fftwf_complex*>(dst.data_ptr),
                                                 (inverse ? FFTW_BACKWARD : FFTW_FORWARD),
                                                 FFTW_ESTIMATE);

            break;
        case 2:
            this->plan_cpu_f_forward = fftwf_plan_dft_2d( width, height,
                                                  reinterpret_cast<fftwf_complex*>(src.data_ptr),
                                                  reinterpret_cast<fftwf_complex*>(dst.data_ptr),
                                                  (inverse ? FFTW_BACKWARD : FFTW_FORWARD),
                                                  FFTW_ESTIMATE);

            break;
        case 3:
            /*     plan_cpu = fftw_plan_dft_3d(int n0, int n1, int n2,
                    fftw_complex *in, fftw_complex *out,
                    int sign, unsigned flags);*/
            // cufftPlan3d(&plan, width ,depth, height , CUFFT_);
        default:
            break;
        }
    }

    void operator()(FFTDataDST &/*dst*/, FFTDataSRC &/*src*/, TransformDirection t_direction)
    {
        bool inverse;
        switch(t_direction)
        {
        case FORWARD:
            inverse = false;
            break;
        case BACKWARD:
            inverse = true;
            break;
        }
        if(!initialized || (old_direction != t_direction))
        {
            reinit(inverse);
            old_direction = t_direction;
        }

        fftwf_execute(this->plan_cpu_f_forward);
        fftwf_execute(this->plan_cpu_f_backward);
    }
};

template <int dim, cufftType_t cufft_type>
struct CUDAFFT<double, dim, cufft_type, cpu>
        : SciPAL::CUDAFFTbase<double, dim, cufft_type, cpu>
{
    typedef typename CUDAFFTbase<double, dim, cufft_type, cpu>::FFTDataDST FFTDataDST;
    typedef typename CUDAFFTbase<double, dim, cufft_type, cpu>::FFTDataSRC FFTDataSRC;

    //creates fftw handle object CPU
    fftw_plan plan_cpu_d_forward;
    fftw_plan plan_cpu_d_backward;
    int width, height;
    FFTDataDST &dst;
    FFTDataSRC &src;
    bool initialized;
    TransformDirection old_direction;
    const QString precision;
    int num_threads;

    CUDAFFT(int _height, int _width, FFTDataDST &_dst, FFTDataSRC &_src): width(_width),
        height(_height), dst(_dst), src(_src), initialized(false), precision("double")
    {
        //initializiation of multithreaded fftw
        dealii::ParameterHandler prm;
        prm.declare_entry("FFTW-Threads", "4",
                          dealii::Patterns::Integer(),
                          "sets the number of Threads, the FFTW uses");
        prm.read_input("fftw-param.prm");
        num_threads = prm.get_integer("FFTW-Threads");
        int error = fftw_init_threads();
#ifdef DEBUG
        std::cout<<"fftw says(0 bad, everything else ok):"<<error<<std::endl;
#endif
        fftw_plan_with_nthreads(num_threads);
        this->check_folder();

        //check if wisdom files exist if not create...
        QString wisdomname = QString("wisdom_%1_%2_%3_%4_forward")
                .arg(width).arg(height).arg(num_threads).arg(precision);
        //gives undefined reference ?!
        //int success = fftw_import_wisdom_from_filename(wisdomname.toStdString().c_str());
        int success = 1;
        if (success != 0)
        { //we already have a wisdom
            this->plan_init();
        }
        else
        { //no wisdom found, create new one
            std::cout<<"creating wisdom"<<std::endl;
          typedef SciPAL::VTraits<double, cpu>::cplxVector cplxVector;
            cplxVector initVec(width * height);
            initVec.reinit(width * height);
            cplxVector initVec2(width * height);//, std::complex<double>(1.0, 1.0));
            initVec2.reinit(width * height);

            fftw_plan plan_forward = fftw_plan_dft_2d(width, height,
                                                    reinterpret_cast<fftw_complex*>(initVec.data_ptr),
                                                    reinterpret_cast<fftw_complex*>(initVec.data_ptr),
                                                    FFTW_FORWARD,
                                                    FFTW_ESTIMATE/*FFTW_PATIENT*/);

//            fftw_export_wisdom_to_filename(QString("wisdom_%1_%2_%3_%4_forward")
//                                           .arg(width).arg(height).arg(num_threads).arg(precision)
//                                           .toStdString().c_str());

            fftw_plan plan_backward = fftw_plan_dft_2d(width, height,
                                                       reinterpret_cast<fftw_complex*>(initVec.data_ptr),
                                                       reinterpret_cast<fftw_complex*>(initVec.data_ptr),
                                                       FFTW_BACKWARD,
                                                       FFTW_PATIENT);

//            fftw_export_wisdom_to_filename(QString("wisdom_%1_%2_%3_%4_backward")
//                                           .arg(width).arg(height).arg(num_threads).arg(precision)
//                                           .toStdString().c_str());
            this->plan_init();
        }

    }

    ~CUDAFFT()
    {
        fftw_destroy_plan(this->plan_cpu_d_forward);
        fftw_destroy_plan(this->plan_cpu_d_backward);
    }

    void plan_init()
    {
        QString wisdomname_for = QString("fftw-wisdom/wisdom_%1_%2_%3_%4_forward")
                .arg(width).arg(height).arg(num_threads).arg(precision);
        QString wisdomname_back = QString("fftw-wisdom/wisdom_%1_%2_%3_%4_backward")
                .arg(width).arg(height).arg(num_threads).arg(precision);

        switch(dim)
        {
        case 1:
//            this->plan_cpu_d = fftw_plan_dft_1d(width,
//                                                reinterpret_cast<fftw_complex*>(src.data_ptr),
//                                                reinterpret_cast<fftw_complex*>(dst.data_ptr),
//                                                (inverse ? FFTW_BACKWARD : FFTW_FORWARD) , FFTW_ESTIMATE);

            break;
        case 2:

//            fftw_import_wisdom_from_filename(wisdomname_for.toStdString().c_str());
            this->plan_cpu_d_forward = fftw_plan_dft_2d(width, height,
                                                reinterpret_cast<fftw_complex*>(src.data_ptr),
                                                reinterpret_cast<fftw_complex*>(dst.data_ptr),
                                                FFTW_FORWARD,
                                                FFTW_PATIENT);

//            fftw_import_wisdom_from_filename(wisdomname_back.toStdString().c_str());
            this->plan_cpu_d_backward = fftw_plan_dft_2d(width, height,
                                                reinterpret_cast<fftw_complex*>(src.data_ptr),
                                                reinterpret_cast<fftw_complex*>(dst.data_ptr),
                                                FFTW_BACKWARD,
                                                FFTW_PATIENT);

            break;
        case 3:
            /*     plan_cpu = fftw_plan_dft_3d(int n0, int n1, int n2,
                    fftw_complex *in, fftw_complex *out,
                    int sign, unsigned flags);*/
            // cufftPlan3d(&plan, width ,depth, height , CUFFT_);
        default:
            break;
        }

    }

    void operator()(FFTDataDST &/*dst*/, FFTDataSRC &/*src*/,
                    TransformDirection t_direction)
    {
        switch(t_direction)
        {
        case FORWARD:
            fftw_execute((this->plan_cpu_d_forward));
            break;
        case BACKWARD:
            fftw_execute((this->plan_cpu_d_backward));
            break;
        }
    }
};

//Specialization for CUDA.
template <typename T, int dim, cufftType_t cufft_type>
struct CUDAFFT<T, dim, cufft_type, gpu_cuda>
        : SciPAL::CUDAFFTbase<T, dim, cufft_type, gpu_cuda>
{
    typedef typename CUDAFFTbase<T, dim, cufft_type, gpu_cuda>::FFTDataDST FFTDataDST;
    typedef typename CUDAFFTbase<T, dim, cufft_type, gpu_cuda>::FFTDataSRC FFTDataSRC;

    //create cufft handle object CUDA
    cufftHandle plan_cuda;
    int height, width;
    int batch2;

    PropagationKernels<typename PrecisionTraits<T, gpu_cuda>::ComplexType, gpu_cuda>
            pkernels;

    CUDAFFT(int _height, int _width, FFTDataDST &/*dst*/, FFTDataSRC &/*src*/, int batch = 1)
        :height(_height), width(_width), batch2(1), pkernels(4)
   {
        cufftResult result;
        switch(dim)
        {
        case 1:
            result = cufftPlan1d(&this->plan_cuda, width, cufft_type, batch);
            break;
        case 2:

            if(width*height > 4096 * 4096)
            {
                uint used_mem = 67108864; //64MB
                batch2 = (used_mem) / (width * 2 * sizeof(T));
                result = cufftPlan1d(&this->plan_cuda, width, cufft_type, batch2);
            }
            else
                result = cufftPlan2d(&this->plan_cuda, width, height, cufft_type);

            break;
        case 3:
            // result = cufftPlan3d(&plan_cuda, width ,depth, height , CUFFT_);
        default:
            break;
        }

        if(result!=CUFFT_SUCCESS)
        {
            switch(result)
            {
            case CUFFT_ALLOC_FAILED:
                std::cout<<"The allocation of GPU resources for the plan failed."<<std::endl;
                break;
            case CUFFT_INVALID_VALUE:
                std::cout<<"One or more invalid parameters were passed to the API."<<std::endl;
                break;
            case CUFFT_INTERNAL_ERROR:
                std::cout<<"An internal driver error was detected."<<std::endl;
                break;
            case CUFFT_INVALID_SIZE:
                std::cout<<"Either or both of the nx or ny parameters is not a supported size."<<std::endl;
                break;
            case CUFFT_SETUP_FAILED:
                std::cout<<"The CUFFT library failed to initialize."<<std::endl;
                break;
            case CUFFT_SUCCESS:
                break;
            default :
                std::cout<<"Unspecified error"<<std::endl;
            }
        }
    }

//    CUDAFFT(int _height, int _width, FFTDataDST &/*dst*/, FFTDataSRC &/*src*/):
//        height(_height), width(_width)
//    {
//        std::vector<int> size{height, width};

//        cufftPlanMany(&this->plan_cuda, 2, &size[0],
//                NULL, 1, width,
//                NULL, 1, width,
//                cufft_type, height);

//    }

    ~CUDAFFT()
    {
        cufftResult result = cufftDestroy(this->plan_cuda);

        if(result!=CUFFT_SUCCESS)
        {
            switch(result)
            {
            case CUFFT_INVALID_PLAN:
                std::cout<<"CUFFT is passed an invalid plan handle."<<std::endl;
                break;
            case CUFFT_SETUP_FAILED:
                std::cout<<"The CUFFT library failed to initialize."<<std::endl;
                break;
            case CUFFT_SUCCESS:
                break;
            }

        }
    }

    inline
    void operator()(FFTDataDST &dst, FFTDataSRC &src, TransformDirection t_direction)
    {
        if((src.leading_dim != src.n_rows) || (dst.leading_dim != dst.n_rows)) {
            std::cerr << "CUDAFFT::operator(): Memory must be tight" << std::endl;
        }

        bool inverse;
        switch(t_direction)
        {
        case FORWARD:
            inverse = false;
            break;
        case BACKWARD:
            inverse = true;
            break;
        }

        cufftResult result;

        if(width*height <= 4096 * 4096)
        {
            result =
            FFTTraits<cufft_type, gpu_cuda>::exec(this->plan_cuda,
                                                  dst.data_ptr,
                                                  src.data_ptr,
                                (inverse ? CUFFT_INVERSE : CUFFT_FORWARD) );
        }
        else
        {
            for(uint i = 0; i < dst.size(); i += width * batch2)
            {
                result =
                FFTTraits<cufft_type, gpu_cuda>::exec(this->plan_cuda,
                                                              &dst.data_ptr[i],
                                                              &src.data_ptr[i],
                                             (inverse ? CUFFT_INVERSE : CUFFT_FORWARD) );
            }



            pkernels.transpose2(dst.data_ptr, width, height, width*height);

            for(uint i = 0; i < dst.size(); i += width * batch2)
            {
                result =
                FFTTraits<cufft_type, gpu_cuda>::exec(this->plan_cuda,
                                                      &dst.data_ptr[i],
                                                      &src.data_ptr[i],
                                             (inverse ? CUFFT_INVERSE : CUFFT_FORWARD) );
            }
        }
//        typedef typename FFTTraits<cufft_type, gpu_cuda>::T_dst T_dst;
//        SciPAL::One<T> one;
//        std::complex<T> OnE(one());
//        FFTDataDST tmp; tmp.reinit(dst.size());
//        cublas::geam(CUBLAS_OP_T, CUBLAS_OP_N, height, width,
//             &OnE, dst.data_ptr, width,
//             &OnE, NULL, width,
//             tmp.data_ptr, width);

//        result = SciPAL::FFTTraits<cufft_type, gpu_cuda>::exec(this->plan_cuda,
//                                                                           dst.data_ptr,
//                                                                           src.data_ptr,
//                                                                           (inverse ? CUFFT_INVERSE : CUFFT_FORWARD) );

//        cudaDeviceSynchronize();

        //Output of possible error messages.
        if(result!=CUFFT_SUCCESS)
        {
            switch(result)
            {
            case CUFFT_INVALID_PLAN:
                std::cout<<"CUFFT is passed an invalid plan handle."<<std::endl;
                break;
            case CUFFT_ALLOC_FAILED:
                std::cout<<"CUFFT failed to allocate GPU memory."<<std::endl;
                break;
            case CUFFT_INVALID_TYPE:
                std::cout<<"The user requests an unsupported type."<<std::endl;
                break;
            case CUFFT_INVALID_VALUE:
                std::cout<<"The user specifies a bad memory pointer."<<std::endl;
                break;
            case CUFFT_INTERNAL_ERROR:
                std::cout<<"Internal driver errors."<<std::endl;
                break;
            case CUFFT_EXEC_FAILED:
                std::cout<<"CUFFT failed to execute an FFT on the GPU."<<std::endl;
                break;
            case CUFFT_SETUP_FAILED:
                std::cout<<"The CUFFT library failed to initialize."<<std::endl;
                break;
            case CUFFT_INVALID_SIZE:
                std::cout<<"The user specifies an unsupported FFT size."<<std::endl;
                break;
            case CUFFT_UNALIGNED_DATA:
                std::cout<<"Input or output does not satisfy texture alignment requirements."<<std::endl;
                break;
            case CUFFT_SUCCESS:
                break;
            }

        }
    }
};//end struct

} //namespace SciPAL end

//These explicit templates are enough to force the compiler to create the remaining
// instantiations by itself.
/*template class SciPAL::CUDAFFT<double, 2,CUFFT_Z2Z, cpu>;
template class SciPAL::CUDAFFTbase<double, 2,CUFFT_Z2Z, cpu>;

template class SciPAL::CUDAFFT<float, 2,CUFFT_C2C, cpu>;
template class SciPAL::CUDAFFTbase<float, 2,CUFFT_C2C, cpu>;*/

template class SciPAL::CUDAFFT<double, 2,CUFFT_Z2Z, gpu_cuda>;
template class SciPAL::CUDAFFTbase<double, 2,CUFFT_Z2Z, gpu_cuda>;

template class SciPAL::CUDAFFT<float, 2,CUFFT_C2C, gpu_cuda>;
template class SciPAL::CUDAFFTbase<float, 2,CUFFT_C2C, gpu_cuda>;

//template class SciPAL::CUDAFFT<float, 2,CUFFT_R2C, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<float, 2,CUFFT_R2C, gpu_cuda>;

//template class SciPAL::CUDAFFT<float, 2,CUFFT_C2R, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<float, 2,CUFFT_C2R, gpu_cuda>;

//inverse transforms
//template class SciPAL::CUDAFFT<double, 2,CUFFT_Z2Z, cpu>;
//template class SciPAL::CUDAFFTbase<double, 2,CUFFT_Z2Z, cpu>;

//template class SciPAL::CUDAFFT<float, 2,CUFFT_C2C, cpu>;
//template class SciPAL::CUDAFFTbase<float, 2,CUFFT_C2C, cpu>;


//template class SciPAL::CUDAFFT<double, 1, CUFFT_Z2D, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<double, 1,CUFFT_Z2D, gpu_cuda>;

//template class SciPAL::CUDAFFT<float, 1, CUFFT_C2R, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<float, 1, CUFFT_C2R, gpu_cuda>;

//template class SciPAL::CUDAFFT<double, 1, CUFFT_D2Z, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<double, 1,CUFFT_D2Z, gpu_cuda>;

//template class SciPAL::CUDAFFT<float, 1, CUFFT_R2C, gpu_cuda>;
//template class SciPAL::CUDAFFTbase<float, 1, CUFFT_R2C, gpu_cuda>;

#endif // CUDA_FFT_H
