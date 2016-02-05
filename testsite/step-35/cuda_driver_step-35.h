//@sect3{File: cuda_driver_step-35.h}
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

Copyright  Lutz Künneke, Jan Lebert, Stephan Kramer, Johannes Hagemann 2014-2015
*/

#ifndef CUDADriver_STEP_35_H
#define CUDADriver_STEP_35_H



#define USE_DEBLUR

#define DISABLE_PSF
#undef DISABLE_PSF


#define USE_FFT_CONV
#undef USE_FFT_CONV

#ifdef USE_FFT_CONV
#define USE_SHIFT_ARITH
#define SHIFT_FACTOR 1
#else
#define SHIFT_FACTOR 0
#endif

//std
#include<iostream>
#include <cmath>

//CUDA
#include <cuda_runtime.h>

//CUDA cufft
#include <cufft.h>

//shift for fft
#include <cufftShiftInterface.h>

//SciPAL
#include <base/PrecisionTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/cublas_Vector.h>

//Our stuff
#include "cuda_driver_step-35.hh"
#include <step-35/cuda_kernel_wrapper_step-35.cu.h>
#include "patch.h"
#include "cuda_helper.h"
#include "preprocessor_directives.h"
#include "ADMMParams.h"
#include <fftw3.h>
#include <smre_problem.hh>

//#include <csignal>

namespace step35 {
// We encapsulate each project into a dedicated namespace
// in order to be able to re-use parts of a test program in others.

//Forward declarations
template<typename Mpatch,typename T,ParallelArch c> class CUDADriver;
template<typename Mpatch,typename T,ParallelArch c> struct constructStruct;
template<typename Mpatch,typename T,ParallelArch c> struct runStruct;

//@sect4{Struct: constructStruct}
//template specialized struct for the constructor and destructor of
//CUDADriver class based on which Mpatch we're dealing with


// @sect4{Class: CUDADriver}
//
// Main class, sets up an instance of queue and info class. A number of threads is started which
// each handles a cuda stream and processes items obtained by the queue. The main thread adds items
// to the queue until everything is processed.
template<typename Mpatch,typename Mdouble, ParallelArch c>
class CUDADriver {
    friend  struct runStruct<Mpatch,Mdouble,c>;
public:
    //complex type: double2 for double, float2 for float
    //    typedef typename PrecisionTraits<Mdouble, c>::ComplexType complex;
    typedef CudaComplex<Mdouble> complex;

    typedef blas BW;

    typedef  SciPAL::Vector<Mdouble, BW> Vector;


private:
    // The data of the SMRE problem stored on the device.
    Vector __im_d, __e_d;

public:
    Vector x_d, z_d, x_old;

    dealii::Vector<Mdouble> x_h, e_h, z_h, im_h, generated_noise;

    // Access to noise contribution.
    const SciPAL::Vector<Mdouble,BW> & e_d() const { return __e_d; }

    SciPAL::Vector<Mdouble,BW> & writeable_e_d() { return __e_d; }

    // The input image should not be touched anywhere in the program. Therefore we only grant read access.
    const SciPAL::Vector<Mdouble,BW> & im_d() const { return __im_d; }

    // FIXME: remove:
    SciPAL::Vector<Mdouble,BW> & writeable_im_d() { return __im_d; }

    // FIXME: For the device side arrays use this:
    SciPAL::Vector<Mdouble, BW> tmp_d, tmp2_d, lag1, lag2, tmp_haar2, tmp_lagr, tmp_haar, tmp_lag1, tmp_lag2;



    // We put the convolution into a structure such tht we can use as a matrix in iteratuve solvers.
    struct Convolution {

        int depth, width, height;
        int cframesize; // image size for r2c images
        bool init; // for fftw

        cufftHandle *plan_fft,*iplan_fft;
        //        CUDAFFT<T, dim, cufft_type, gpu_cuda> FFT; //FIXME use that

        fftwf_plan plan_cpu_d_forward;
        fftwf_plan plan_cpu_d_backward;

        SciPAL::Vector<complex, BW> fm1;
        //Point spread function on device
        SciPAL::Vector<complex, BW> psf_fourier_transform_d;

        dealii::Vector<Mdouble> psf_h;

        // In cases where we only want to denoise
        // we have to bypass the convolution operation as then the blurrign operator is the identity.
        bool is_delta_peak;


        ~Convolution() {
#ifdef USE_FFT_CONV
            cufftDestroy(*plan_fft);
            cufftDestroy(*iplan_fft);
            delete plan_fft;
            delete iplan_fft;
            fftwf_destroy_plan(this->plan_cpu_f_forward);
            fftwf_destroy_plan(this->plan_cpu_f_backward);
#endif
        }

        Convolution(int d, int w, int h, double fwhm)
            :
              depth(d), width(w), height(h),
              cframesize(/*ext_*/depth*/*ext_*/width*(/*ext_*/height/2+1)/**sizeof(Mcomplex)*/),
              fm1(cframesize), init(false),
              psf_fourier_transform_d(fm1.size()),
              psf_h(width*height),
              is_delta_peak(false),
      #ifndef USE_FFT_CONV
              dof_handler(dealii::Point<3>(width, height, depth))
    #endif
        {
          {
                // TO DO: compute psf_h here

                SciPAL::Vector<Mdouble,BW> fpsf_tmp_d(psf_h.size());

                Mdouble psf_norm=0;

                for (int z=0; z<depth; z++ ) {
                    for (int x=0; x<width; x++) {
                        for (int y=0; y<height; y++) {
                            psf_h(z*height*width+x*height+y) = _psf(Mdouble(x), Mdouble(y), Mdouble(z), fwhm );
                            psf_norm += psf_h[z*height*width+x*height+y];
                        }
                    }
                }

                //norm
                psf_h *= 1. / psf_norm;
                //                                psf_h *= 65535;
                fpsf_tmp_d = psf_h;
                //                cufftShift_2D_impl(fpsf_tmp_d.data(), width, height);
                //Cufft setup
                //Create cufft plans
                plan_fft=new cufftHandle();
                iplan_fft=new cufftHandle();

#ifdef DOUBLE_PRECISION
                cufftPlan3d(plan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_D2Z);
                cufftPlan3d(iplan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_Z2D);
#else
                cufftPlan3d(plan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_R2C);
                cufftPlan3d(iplan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_C2R);
#endif


            }

#ifdef USE_FFT_CONV
            // TO DO: compute psf_h here
            SciPAL::Vector<Mdouble,BW> fpsf_tmp_d;
            fpsf_tmp_d = (psf_h);



            //Cufft setup
            //Create cufft plans
            plan_fft=new cufftHandle();
            iplan_fft=new cufftHandle();
#ifdef DOUBLE_PRECISION
            cufftPlan3d(plan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_D2Z);
            cufftPlan3d(iplan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_Z2D);
#else
            cufftPlan3d(plan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_R2C);
            cufftPlan3d(iplan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_C2C);
#endif


            if(typeid(BW) = typeid(cublas))
            {
                //Fourier transformation of psf using cuFFT
                cufftHandle plan;
#ifdef DOUBLE_PRECISION
                cufftPlan3d(&plan, ext_depth, ext_width, ext_height, CUFFT_D2Z);
                cufftExecD2Z(plan, fpsf_tmp_d, fpsf_d);
#else
                cufftPlan3d(&plan, /*ext_*/depth, /*ext_*/width, /*ext_*/height,  CUFFT_R2C);
                cufftExecR2C(plan, fpsf_tmp_d.data(), psf_fourier_transform_d.data());
#endif
                //Cleanup
                cufftDestroy(plan);
                getLastCudaError("cufft error!\n");
            } else
            {
                fftwf_plan plan =
                        fftwf_plan_dft_r2c_2d(width, height,
                                              fpsf_tmp_d.data(),
                                              reinterpret_cast<fftwf_complex*>(psf_fourier_transform_d.data()),
                                              FFTW_PATIENT);
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
            }
#else

            // Compute 1D PSF
            if (fwhm >0)
            {
                cut_off = std::ceil(1.2*fwhm); // corresponds approximately to a cut off radius of 3 std devs

                // In case of a PSF to narrow to be reasonably sampled on the lattice of pixels
                // we simply assume that it is a delta peak and consider the convolution operator as identity.

                psf_1d.reinit(2*cut_off + 1);

                // The magic number for converting FWHM to a standard deviation is taken
                // from <a href="http://mathworld.wolfram.com/GaussianFunction.html">wolfram.com</a>.
                double sigma_psf = std::max(1e-12, fwhm/2.3548);
                int x_0 = 0;
                double psf_norm = 0;

                for (int x = -cut_off; x <= cut_off; x++)
                {
                    psf_1d(x+cut_off) = exp(-(boost::math::pow<2>(x - x_0)  )/
                                            (2*boost::math::pow<2>(sigma_psf)));
                    psf_norm += psf_1d(x+cut_off);
                }
                psf_1d /= psf_norm;

                std::cout << psf_1d << std::endl;
            }
            else
                is_delta_peak = true;
#endif

        }
#ifdef USE_FFT_CONV
#else
        dealii::Vector<Mdouble> psf_1d;
        int cut_off;
        ImageDoFHandler dof_handler;
#endif

        void vmult(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            if (is_delta_peak)
                this->vmult_id(dst, src);
            else
                this->vmult_FFT(dst, src);
        }

        //    private:
        void vmult_id(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            dst = src;
            return;
        }

        void vmult_conv(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            // Assert that dst is large enough.
            // We assume that the image space has
            // the same dimension as the range space.
            if (dst.size() == 0)
                dst.reinit(src.size());


#ifdef USE_FFT_CONV
#ifdef DOUBLE_PRECISION
            cufftExecD2Z(*plan_fft, in, fm1.data());
#else
            cufftExecR2C(*plan_fft, const_cast<Mdouble*>(src.data()), fm1.data());
#endif
            checkCudaErrors(cudaDeviceSynchronize());

            //Convolve, multiply in Fourier space
            step35::Kernels<Mdouble> kernel;
            kernel.element_norm_product(fm1.data(), psf_fourier_transform_d.data(),
                                        /*inf->ext_*/width, /*inf->ext_*/height, /*inf->ext_*/depth);

            //Transform back
            // FIXME: replace by SciPAL's FFT wrappers
#ifdef DOUBLE_PRECISION
            cufftExecZ2D(*iplan_fft, fm1.data(), out);
#else
            cufftExecC2R(*iplan_fft, fm1.data(), dst.data());
#endif
            checkCudaErrors(cudaDeviceSynchronize());
            getLastCudaError("cufft error!\n");
#else

            dealii::Vector<Mdouble> tmp_dst(src.size()), tmp_src(dst.size());
            src.push_to(tmp_dst);

            // Number of rows before and after the current row
            // which have to be taken into account.
            int N = cut_off;
            for (int k=0; k< depth; k++)
                for (int i=0; i< width; i++)
                    for (int j=0; j< height; j++)
                    {
                        int xyz_pos = dof_handler.global_index(i, j, k);

                        tmp_src(xyz_pos) = 0.;

                        for (int p = -N; p <= N; p++)
                        {
                            if (j+p< 0 || j+p>=height)
                                continue;
                            int xyz_pos_p = dof_handler.global_index(i, j+p, k);
                            tmp_src(xyz_pos) += psf_1d[p+N]*tmp_dst(xyz_pos_p);
                        }
                    }

            // We abuse the host-side source array as a temporary one.
            // To avoid unphysical oscillations at the boundaries of the image we have to exclude a boundary layer from the convolution.
            for (int k=0; k< depth; k++)
                for (int i=N; i< width-N; i++)
                    for (int j=N; j< height-N; j++)
                    {
                        int xyz_pos = dof_handler.global_index(i, j, k);
                        tmp_dst(xyz_pos) = 0.;
                        for (int p = -N; p <= N; p++)
                        {
                            if (i+p< 0 || i+p>=width)
                                continue;

                            int xyz_pos_p = dof_handler.global_index(i+p, j, k);
                            tmp_dst(xyz_pos) += psf_1d[p+N]*tmp_src(xyz_pos_p);
                        }
                    }

            dst = tmp_dst;
#endif
        }

        void vmult_FFT(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            // Assert that dst is large enough.
            // We assume that the image space has
            // the same dimension as the range space.

            // setup fftw on first run
            //fftw plans
            if(init == false)
            {
                init = true;
//                fftwf_plan_with_nthreads();
                this->plan_cpu_d_forward
                        = fftwf_plan_dft_r2c_2d(width, height,
                                                (src.data()),
                                                reinterpret_cast<fftwf_complex*>(fm1.data()),
                                                FFTW_PATIENT);

                this->plan_cpu_d_backward
                        = fftwf_plan_dft_c2r_2d(width, height,
                                                reinterpret_cast<fftwf_complex*>(fm1.data()),
                                                (dst.data()),
                                                FFTW_PATIENT);
            }
            if (dst.size() == 0)
                dst.reinit(src.size());

            //            cufftShift_2D_impl(src.data(), width, height);
            if(typeid(BW) == typeid(cublas)) {
#ifdef DOUBLE_PRECISION
                cufftExecD2Z(*plan_fft, src.data(), fm1.data());
#else
                cufftExecR2C(*plan_fft, src.data(), fm1.data());
#endif
                checkCudaErrors(cudaDeviceSynchronize());
            } // cublas case
            else
            {
                fftwf_execute(plan_cpu_d_forward);
            }

            // multiplication in fourier space
            fm1 = fm1 && psf_fourier_transform_d;

            //Transform back
            // FIXME: replace by SciPAL's FFT wrappers
            if(typeid(BW) == typeid(cublas)) {
#ifdef DOUBLE_PRECISION
                cufftExecZ2D(*iplan_fft, fm1.data(), fm1.data(), CUFFT_INVERSE);
#else
                cufftExecC2R(*iplan_fft, fm1.data(), dst.data());
#endif
                checkCudaErrors(cudaDeviceSynchronize());
                getLastCudaError("cufft error!\n");
            }
            else
            {
                fftwf_execute(plan_cpu_d_backward);
            }

            dst *= 1./dst.size();
            //ugly...
            kernelConf* conf = (kernelConf*) malloc(sizeof(kernelConf));
            int threadsPerBlock_X = 16;
            int threadsPerBlock_Y = 16;
            conf->block = dim3(threadsPerBlock_X, threadsPerBlock_Y, 1);
            conf->grid = dim3(((width/2) / threadsPerBlock_X), ((width/2) / threadsPerBlock_Y), 1);
            cufftShift_2D_config_impl(dst.data(), width, height, conf);

            //            dst = SciPAL::abs(fm1);

        }

        //@sect5{Function: _psf}
        //@brief Returns a 2D gaussian convolution kernel in direct space
        Mdouble _psf(Mdouble x,Mdouble y,Mdouble z, Mdouble fwhm) {

            double sigma = std::max(1e-12, fwhm/2.3548);

            //Fixme
            //            if ( params.dim == 2 )
            //            {
            if ( z < 0.5 )//float conditional to distinguish z=0 vs. z=1
            {
                Mdouble tmp = exp(-(boost::math::pow<2>(x - width/2.) + boost::math::pow<2>(y - height/2.))/
                                  (2.0*boost::math::pow<2>(sigma/2.0)));
                //                if (tmp > 0.01)
                //                    std::cout<<"psf "<<tmp<<std::endl;

                return tmp;
            }
            else
                return 0;
            //            }
            //            else
            {

                //                    return exp(-(boost::math::pow<2>(x - sigma) + boost::math::pow<2>(y - sigma) +
                //                                 boost::math::pow<2>(z - sigma))/(TWO*boost::math::pow<2>(HALF*sigma)));
            }
        }

    }; // end of class definition of convolution


    Convolution convolution;

    // Compute the residual $I - \epsilon - A * x$. To be put into an expression.
    void residual(SciPAL::Vector<Mdouble, BW> & res,
                  const SciPAL::Vector<Mdouble, BW> & im,
                  const SciPAL::Vector<Mdouble, BW> & noise,
                  const SciPAL::Vector<Mdouble, BW> & est_im)
    {

        // A * x
        convolution.vmult(res, est_im);

#ifndef nFMM_VERSION
        // A * x + eps
        res += noise;
        // A * x + eps -I
        res -= im;

        // I - eps - A * x
        res *= -1;
#else
        res -= this->lag2; // noise;


#endif

    }


    ImageDoFHandler dof_handler;

    const int n_max_dykstra_steps;

    const Mdouble dykstra_Tol;

    const Mdouble sigma_noise;

    const int dim;

    Mdouble inv_gamma;

    regularisationType regType;


#ifdef TIME_KERNELS
    // Total GPU time for dykstra kernels in ms
    float kernel_time;
    // Output file where kernel execution time is written at shutdown
    std::ofstream out_time;
#endif

    // @sect5{Constructor: CUDADriver}
    //
    // @param cs_h array of $c_s$
    // @param dofs DoFHandler, provides information about image size.
    // @param params the list of runtime parameters
    // FIXME: cs_h must some a reference to some decent type not a raw pointer
    CUDADriver(std::vector<Mdouble> & cs_h,
               // FIXME: use ImageStack class developed for SOFI and ask that for nx, ny, nz, if needed
               const ImageDoFHandler& dofs,
               const ADMMParams& params)
        :
          x_h(dofs.n_dofs()),
          e_h(dofs.n_dofs()),
          z_h(dofs.n_dofs()),
          im_h(dofs.n_dofs()),
          __e_d(dofs.n_dofs()),
          x_d(dofs.n_dofs()),
          z_d(dofs.n_dofs()),
          __im_d(dofs.n_dofs()),
          //
          dof_handler(dofs),
          tmp_d(dofs.n_dofs()),
          tmp2_d(  dofs.n_dofs()),
          lag1(  dofs.n_dofs()),
          lag2(  dofs.n_dofs()),
          tmp_haar2(  dofs.n_dofs()),
          tmp_haar(  dofs.n_dofs()),
          tmp_lagr(  dofs.n_dofs()),
          convolution(dofs.pdepth(),  dofs.pwidth(),  dofs.pheight(), params.psf_fwhm),
          n_max_dykstra_steps(params.n_max_dykstra_steps),
          dykstra_Tol(params.dykstra_Tol),
          sigma_noise(params.gnoise),
          dim(params.dim),
          inv_gamma(params.inv_gamma),
          regType(params.regType)
    {

#ifdef TIME_KERNELS
        kernel_time = 0;
        //Write total dykstra kernel execution time at shutdown to file "kernel_time"
        out_time.open("kernel_time");
#endif


    }

    //@sect5{Destructor: CUDADriver}
    //Signal all threads to end, calls queue destructor
    ~CUDADriver () {



#ifdef TIME_KERNELS
        //Write total dykstra kernel execution time at shutdown to file "kernel_time"
        if (kernel_time)
            out_time << kernel_time << std::endl;
        out_time.close();
#endif

    }

    //@sect5{Function: nostream_handler}
    //@brief Wrapper if we use the dyadyc_dykstra Kernel.
    void iterate() {
#ifdef TIME_KERNELS
        cudaEvent_t start, stop;
        float time;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        cudaEventSynchronize(stop);
#endif

        step35::Kernels<Mdouble> kernels;

        if (this->dim == 2) {
            //Do the approximate method shifted by 2^so. This choice of shifts enables us to omit some sets in later instances of the kernel
            for (int so=0; so< 1 // 5
                 ; so++)
            {
                for (int z=0; z< dof_handler.pdepth(); z++) {
                    kernels.dyadic_dykstra(this->writeable_e_d().data(),
                                           this->im_d().data(),
                                           this->sigma_noise,
                                           dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth(),
                                           n_max_dykstra_steps, dykstra_Tol);
                }
                //  if (so == 1) AssertThrow(false, dealii::ExcInternalError());
            }
        }
        else { //TODO 3d
            for (int so=0; so<5; so++) {
                for (int z=0; z< dof_handler.pdepth(); z++) {
                    kernels.dyadic_dykstra(this->writeable_e_d().data(),
                                           this->im_d().data(),
                                           this->sigma_noise,
                                           dof_handler.pwidth(),  dof_handler.pheight(),  dof_handler.pdepth(), n_max_dykstra_steps, dykstra_Tol);
                }
            }
        }

#ifdef TIME_KERNELS
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time, start, stop);
        kernel_time += time;
#endif


    }

    //@sect5{Function: projection_gauss}
    //@brief CPU implementation of dykstras algorithm with Gauss projection, performs projection for a single set
    //@param d_q pointer to vector q
    //@param im_h pointer to image
    //@param cs_h array of weights
    //@param width width of the image
    //@param height height of the image
    //@param n edge length of the set
    //@param i0 starting row of the set
    //@param j0 starting col of the set
    void projection_gauss(Mdouble* d_q, Mdouble* im_h, Mdouble* cs_h, const int width, const int height,
                          const int n, const int i0, const int j0) {
        //Do the difference and square
        int size=n*n;
        std::vector<Mdouble> f(size);
        Mdouble square_sum;
        int is,js;
        //Square the vector
        for (int i=0; i<size; i++) {
            is = i/n + i0;
            if ( is >= width )
                std::cout << "Error in gauss Projection, dimension too large " << std::endl;
            js = i%n + j0;
            f[i]=im_h[is*height+js]-d_q[i];
            d_q[i]=-f[i];
            square_sum=square_sum+f[i]*f[i];
        }
        square_sum*=cs_h[n];
        //Check the condition
        if ( square_sum > 1.0 ) {
            square_sum=1.0/sqrt(square_sum);
            for (int i=0; i<size; i++) {
                //Do the projection
                is = i/n + i0;
                js = i%n + j0;
                im_h[is*height+js]=f[i]*square_sum;
                d_q[i]+=f[i];
            }
        } else {
            for (int i=0; i<size; i++) {
                is = i/n + i0;
                js = i%n + j0;
                im_h[is*height+js]=f[i];
                d_q[i]+=f[i];
            }
        }
    }



    //@sect5{Function: push_data}
    //@brief Push all data from host to device
    void push_data() {
        this->writeable_im_d() = this->im_h;
        this->x_d = this->x_h;
        this->z_d = this->z_h;
        this->writeable_e_d() = this->e_h;

    }

    //@sect5{Function: get_data}
    //@brief Pull all data from device to host
    void get_data() {
        this->im_d().push_to(this->im_h);
        this->x_d.push_to(this->x_h);
        this->z_d.push_to(this->z_h);
        this->e_d().push_to(this->e_h);

    }

    // Evaluation of the Euler-Lagrange equation for the minimization w.r.t.
    // to the primal variable, that is the image we want to reconstruct.
    // The mathematical expression is
    // \f{eqnarray}{ \text{dst}
    //                 & = &
    //                        A^T*\left\lbrack     \rho_1 \left(I-e-A*x\right) +\Upsilon_1      \right\rbrack
    //                        -
    //                        \rho_2 \frac{\delta J}{\delta x}(x_{old}) \,.
    //             \f}
    void argmin_x_rhs(SciPAL::Vector<Mdouble, BW> & dst,
                      const SciPAL::Vector<Mdouble, BW> & x_old,
                      const SciPAL::Vector<Mdouble, BW> & lag_1,
                      ADMMParams& params)
    {
        const Mdouble rho1 = params.rho1;

        const Mdouble reg_strength = params.reg_strength;


        // the residual
        // $\text{tmp}_d=I-e-A*x$
        this->residual(this->tmp2_d, this->im_d(), this->e_d(), x_old);

        // $\text{tmp}_d=A*\left(\left(I-e-A*x\right) \rho_1+\Upsilon_1\right)$
#ifndef nFMM_VERSION
        dst = rho1 * this->tmp2_d + lag_1;
#else
        dst =  Mdouble(1./params.inv_gamma) * lag_1 + this->tmp2_d;
        dst /= 0.99;
#endif
        convolution.vmult(dst, dst);

        // Add regularization
        step35::Kernels<Mdouble> kernel;
        kernel.tv_derivative(this->tmp2_d.data(), x_old.data(),
                             this->im_d().data(),
                     #ifndef nFMM_VERSION
                             reg_strength
                     #else
                             reg_strength / 0.99
                     #endif
                             ,
                             dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth());

#ifndef nFMM_VERSION
        this->tmp2_d *= params.rho2;
#else
        this->tmp2_d *= Mdouble((1./params.inv_gamma) / 0.99);
#endif
        dst -= this->tmp2_d;

        // For LJs's version:
        // 2nd constraint
        //            {
        //                  tmp2_d = this->z_d - this->x_old;
        //                  tmp2_d *= Mdouble(params.rho2);
        //                  tmp2_d += tmp_d;
        //                  tmp_d  = tmp2_d - lag2;
        //            }

        // test, whether sign reversal stems from here:
#ifdef nFMM_VERSION
        // It seems that we have to reverse the sign of the rhs for the FMM version to get convergence.
        dst *= -1;
#endif
    }


    void x_step_adaptive(ADMMParams & params)
    {

        const Mdouble rho1 = params.rho1;

        const Mdouble reg_strength = params.reg_strength;

        const Mdouble Tol = 1e-2;
        const int max_dt_iter = 10;

        SciPAL::Vector<Mdouble, BW> & x_tmp = this->tmp_lag2;

        SciPAL::Vector<Mdouble, BW> & k_1 = this->tmp_haar;
        SciPAL::Vector<Mdouble, BW> & k_2 = this->tmp_haar2;

        SciPAL::Vector<Mdouble, BW> & sol_diff = this->tmp_lagr;

        Mdouble err = 2*Tol;
        int n_max_admm_iter = 1000;
        int n_err_iter = 0;
        // For a given Lagrange parameter @p lag1 and a given estimate of the noise we compute
        // a new estimate for how the reconstructed image should look like.
        while (err > Tol)
        {
            // Save old solution.
            this->x_old = this->x_d;

            Mdouble Delta = 2*Tol;

            // Evaluate rhs of Euler-Lagrange for explicit Euler step.
            // This does not change over the course of the timestep adaption.
            // Hence, we have to compute it only once.
            argmin_x_rhs(k_1, this->x_old, this->lag1, params);


            int dt_iter = 1;
            Mdouble delta_t;
            // Th einne rloop does the time-step adaption.
            for (; dt_iter <= max_dt_iter; dt_iter++)
            {
                delta_t = 1./params.inv_gamma;

                // Euler step for intermediate solution.
                x_tmp = delta_t * k_1 + x_old;          // for reconstructing LJ: replace delta_t by Mdouble(1./delta_t)

                // For LJ:
                //this->x_d = x_tmp;
                //return;


                // Evaluate rhs of Euler-Lagrange
                argmin_x_rhs(k_2, x_tmp, this->lag1, params);

                // The difference of the two solutions is given by the differnce of the derivatives
                // and serves the timestep adaption.
                sol_diff = k_2;
                sol_diff -= k_1;
                sol_diff *= 0.5*delta_t;

                // Heun step
                k_2 += k_1;
                this->x_d = Mdouble(0.5*delta_t) * k_2 + x_old;

                // this->z_d = this->x_d;

                Delta = sol_diff.l2_norm()/std::sqrt(sol_diff.size());

                Mdouble new_dt = 0.9 * Tol / Delta;

                // To avoid excessive oscillations we limit the increase of the timestep.
                if (Delta < Tol) {
                    if (new_dt > 2*delta_t)
                        new_dt = 2*delta_t;
                }
                else
                    new_dt = 0.707 * delta_t;


                delta_t = params.inv_gamma = 1./(new_dt);

                if (Delta <= Tol)
                    break;
            }
            //  std::cout << "             n dt adaptions : " << dt_iter << ", new timestep : " << 1./params.inv_gamma << ", Delta : " << Delta << std::endl;

            // Convergence control for the outer iteration.
            tmp2_d = this->x_d;
            tmp2_d -= this->x_old;
            err = tmp2_d.l2_norm()/std::sqrt(tmp2_d.size());
            n_err_iter++;
            if (n_err_iter > n_max_admm_iter)
                break;
        }
        std::cout << "    n ADMM steps : " << n_err_iter << std::endl;


        this->z_d = this->x_d;
    }


    // FIXME: use SciPAL vectors and ETs!!!
    //@sect5{Function: x_step}
    //@brief Performs approximative argmin with respect to $x$
    // FIXME: component wise
    void x_step(const Mdouble rho1,const Mdouble rho2)
    {
        step35::Kernels<Mdouble> kernel;
        //$\text{tmp}_d=I$tmp2
        // // // kernel.reset(tmp_d.data(), inf->ext_num_pix);
        // // // kernel.sum(tmp_d.data(), this->im_d.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp}_d=I-e$
        // kernel.diff(tmp_d.data(), this->e_d().data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);

        // conv2(this->x_d.data(), tmp2_d.data()); // args: in, out

        // the residual
        //$\text{tmp}_d=I-e-A*x$
#ifdef USE_SHIFT_ARITH
        kernel.diff(tmp_d.data(), tmp2_d.data(),
                    SHIFT_FACTOR * inf->sigma,
                    inf->ext_width(),  inf->ext_height(),  inf->ext_depth());
#else

#endif
        // Test x -y
        {
            //            tmp2_d = Mdouble(-1.) * this->x_d + this->z_d;
            //            tmp_d =  this->z_d - this->x_d;
            //            tmp2_d -= tmp_d;
            //            std::cout << "tmp2_d.l2_norm : 0 =?= " << tmp2_d.l2_norm() << std::endl;
        }

        this->residual(tmp2_d, this->im_d(), this->e_d(), this->x_d);

        //$\text{tmp}_d=\left(I-e-A*x\right)*\rho_1$
        // kernel.mult(tmp_d.data(), rho1, inf->ext_num_pix);

        //$\text{tmp}_d=\left((im-e-A*x\right)*\rho_1+\Upsilon_1$
        // kernel.sum(tmp_d.data(), lag1.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp}_d=A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
        tmp_d = rho1 * tmp2_d + lag1;

        // conv2(tmp_d.data(), tmp_d.data());
        // FIXME: inplace works?
        convolution.vmult(tmp_d, tmp_d);

        //$\text{tmp2}_d=z$
        // kernel.reset(tmp2_d.data(), inf->ext_num_pix);
        // kernel.sum(tmp2_d.data(), this->z_d.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);



        //$\text{tmp2}_d=z-x$

        // kernel.diff(tmp2_d.data(), this->x_d.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        // tmp2_d -= this->x_d;

        //$\text{tmp2}_d=z-x$
        tmp2_d = this->z_d - this->x_d;



        //$\text{tmp2}_d=\left((z-x\right)*\rho_2$
        // kernel.mult(tmp2_d.data(), rho2, inf->ext_num_pix);
        tmp2_d *= rho2;

        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
#ifdef USE_SHIFT_ARITH
        kernel.sum(tmp2_d.data(), tmp_d.data(),
                   SHIFT_FACTOR * inf->sigma,
                   inf->ext_width(),  inf->ext_height(),  inf->ext_depth());
#else
        tmp2_d += tmp_d;
#endif

        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2$
        // kernel.diff(tmp2_d.data(), lag2.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        tmp2_d -= lag2;

        //$\text{tmp2}_d=\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        // kernel.mult(tmp2_d.data(), inf->gamma, inf->ext_num_pix);

        // tmp2_d *= inf->gamma;


        // THE FOLLOWING DOC EQ APPLIES TO HAVING DONE BEFORE : tmp2_d *= inf->gamma;
        //$x=x+\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        //kernel.sum(this->x_d.data(), tmp2_d.data(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        //  this->x_d += tmp2_d;
        this->x_d = this->inv_gamma * tmp2_d + this->x_d;

    }

    //@sect5{Function: update_lagrangian}
    //@brief Update the lagrangian estimates
    //$ \Upsilon_1^{r+1} = \Upsilon_1^r + \alpha_1 \left( I - A*x - \epsilon \right) $ \n
    //$\Upsilon_2^{r+1} = \Upsilon_2^r + \alpha_2 \left( x - z \right) $
    void update_lagrangian(Mdouble alpha1, Mdouble alpha2) {
        step35::Kernels<Mdouble> kernel;

        // kernel.reset(tmp_d.data(), inf->ext_num_pix);
        // kernel.sum(tmp_d.data(), this->x_d.data(),0, inf->ext_width, inf->ext_height, inf->ext_depth);

        //  conv2(tmp_d.data(), tmp_d.data());
        //Update the lagrangian estimates
#ifdef USE_SHIFT_ARITH
        kernel.update_lagrangian(lag1.data(), lag2.data(),
                                 SHIFT_FACTOR *  inf->sigma,
                                 inf->ext_width(),  inf->ext_height(),  inf->ext_depth(),
                                 alpha1, alpha2, this->e_d().data(), inf->writeable_im_d().data(),tmp_d.data(), this->x_d.data(), this->z_d.data());
#else


        // lag1 += alpha1 * (this->im_d() - this->e_d() - tmp_d);
        residual(tmp_lag1, this->im_d(), this->e_d(), this->x_d);
        lag1 =  alpha1 * tmp_lag1 + lag1;

#ifndef nFMM_VERSION

        // lag2 += alpha2 * (this->x_d - this->z_d);

        tmp_lag2 = this->x_d - this->z_d;
        lag2  = alpha2 * tmp_lag2 + lag2;
#endif
#endif
    }

    //@sect5{Function: dykstra_gauss}
    //@brief Wrapper around the whole projection procedure
    //@param rho parameter for the first constraint
    void dykstra_gauss(const Mdouble rho)
    {

#ifdef USE_SHIFT_ARITH
        convolution.vmult(tmp_d, this->x_d); // conv2(this->x_d.data(),tmp_d.data());


        step35::Kernels<Mdouble> kernel;
        //the value to be projected is copied into e_d

        kernel.prepare_e(inf->writeable_e_d().data(), this->im_d().data(), tmp_d.data(), lag1.data(), rho,
                         SHIFT_FACTOR * inf->sigma,
                         inf->ext_width(),  inf->ext_height(),  inf->ext_depth());
#else

#ifndef nFMM_VERSION
        dykstra_gauss_global_mem(rho); // dykstra_gauss_LNCS(rho);
#else
        dykstra_gauss_FMM(rho);
#endif
        // Will the generated noise be transformed?
        // this->writeable_e_d() = this->generated_noise;

#endif




    }




    void dykstra_gauss_LNCS(const Mdouble rho)
    {
        // $A * x_{k}$
        convolution.vmult(tmp_d, this->x_d);

        // $\epsilon = I -  A * x_{k} + \Upsilon_1/\rho$
        this->writeable_e_d() = Mdouble(1./rho) * lag1 + this->im_d();
        this->writeable_e_d() -= tmp_d;

        //Perform the Dykstra Algorithm
        iterate();

    }


    void dykstra_gauss_FMM(const Mdouble rho)
    {
        // $A * x_{k-1}$
        convolution.vmult(tmp_d, this->x_d);

        // $\epsilon = dt \cdot \Upsilon_1 + A * x_{k-1} $
        this->writeable_e_d() = Mdouble(1./this->inv_gamma) * lag1 + tmp_d;

        //Perform the Dykstra Algorithm
        this->iterate();

        // After the Dykstra iteration @p e_d contains the slack variable $v$ for the constraint $A * x - v = 0$
        this->lag2 = this->e_d();

        // compute noise contribution
        this->writeable_e_d() -= this->im_d();
    }


    void dykstra_gauss_global_mem(const Mdouble rho)
    {
        convolution.vmult(tmp_d, this->x_d);

        // $\epsilon = I -  A * x_{k} + \Upsilon_1/\rho$
        this->writeable_e_d() = Mdouble(1./rho) * lag1 + this->im_d();
        this->writeable_e_d() -= tmp_d;

        static const int n_scales = 10;
        static const int min_scale = 0;

        typedef SciPAL::Vector<Mdouble, BW> Vc;

        std::vector<Mdouble> zero_init(this->e_d().size());

        Vc h_iter (this->e_d().size()); h_iter = zero_init;

        Vc h_old (this->e_d().size()); h_old = zero_init;

        Vc h_init (this->e_d().size()); h_init = zero_init;

        Vc Q_full (this->e_d().size()); Q_full = zero_init;
        Vc Q_M (this->e_d().size()); Q_M = zero_init;

        // $h_0 = h$;
        h_init = this->e_d();
        h_old = h_init;

        Mdouble cs_weight =  //9.46624e-07; // 1024x1024, 0.9, c_S
                9.51511e-07; // 1024x1024, 0.9, c_Ss
        //9.53509e-07; // 1024x1024, 0.1, c_Ss

        int iter = 0;
        bool iterate = true;

        step35::Kernels<Mdouble> kernel;

        while (iterate)
        {

            // Start with estimator for whole image
            h_old -= Q_full;
            h_init = h_old;
            h_iter = h_old;
            // assume single subset for a moment
            Mdouble L2_norm = h_old.l2_norm(); // i.e. sqrt(sum_of_squares)

            // project
            h_iter *= this->sigma_noise/(L2_norm * std::sqrt(cs_weight) );

            Q_full = h_iter - h_old;
            h_old = h_iter;

            // Q_0 <- Q_full
            // initialize : d.h_init with global h_init

            // loop over subsets
            if(typeid(BW) == typeid(cublas))
            {
                kernel.dyadic_dykstra_fine_scale_part(h_iter.data(), h_old.data(),
                                                      Q_full.data(),
                                                      this->writeable_e_d().data(),
                                                      //    this->im_d().data(),
                                                      this->sigma_noise,
                                                      dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth(),
                                                      n_max_dykstra_steps, dykstra_Tol
                                                      );
            }
            else
            {
                kernel.dyadic_dykstra_fine_scale_part_cpu(h_iter.data(), h_old.data(),
                                                          Q_full.data(),
                                                          this->writeable_e_d().data(),
                                                          //    this->im_d().data(),
                                                          this->sigma_noise,
                                                          dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth(),
                                                          n_max_dykstra_steps, dykstra_Tol
                                                          );
            }

            // Convergence control.
            {
                h_init -= h_iter;
                Mdouble norm = h_init.l2_norm();

                iterate = ( norm > dykstra_Tol &&
                            iter < n_max_dykstra_steps);

                iterate ? h_init = h_iter : this->writeable_e_d() = h_iter;
                iter++;
            }
        }
    }


    //@sect5{Function: m_smoothing}
    //@brief Smoothing step
    //@param gamma Strength of the regularization, small choices render the algorithm unstable
    //TODO rho2
    void m_smoothing(Mdouble reg_strength, Mdouble rho2) {
        step35::Kernels<Mdouble> kernel;
        //Regularization by Haar Wavelet Coefficient sparsity
        if ( this->regType == haar ) {
            if ( dof_handler.pdepth() > 1) {//TODO 3d
                std::cerr << "Haar regularisation with 3d images is not yet implemented!" << std::endl;
                std::abort();
            }

            int nx2 = dof_handler.n_dofs_x_padded();
            int ny2 = dof_handler.n_dofs_y_padded();
            //kernel.reset(tmp_haar,inf->nx2*inf->ny2);
            tmp_haar = SciPAL::Vector<Mdouble,BW>( nx2 * ny2);
            tmp_lagr = SciPAL::Vector<Mdouble,BW>( tmp_haar.size());

            //Copy $x$ into the bigger temp variable while conserving its shape
            for (int i=0; i< nx2; i++) {
                if ( i <  dof_handler.pheight()) {
                    checkCudaErrors(cudaMemcpyAsync(&(tmp_haar.data()[i* ny2]), &(this->x_d.data()[i* dof_handler.pheight()]),
                            dof_handler.pwidth()*sizeof(Mdouble), cudaMemcpyDeviceToDevice));

                    checkCudaErrors(cudaMemcpyAsync(&(tmp_lagr.data()[i* ny2]), &(lag2.data()[i* dof_handler.pheight()]),
                            dof_handler.pwidth()*sizeof(Mdouble), cudaMemcpyDeviceToDevice));
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());

            //Forward 2D Haar Wavelet transform
            kernel.haar(tmp_haar.data(), tmp_haar2.data(), ny2);
            kernel.haar(tmp_lagr.data(), tmp_haar2.data(), ny2);
            kernel.soft_threshold(tmp_haar.data(), lag2.data(), tmp_haar.data(), rho2, reg_strength, nx2* ny2);
            //Backward 2D Haar Wavelet transform
            kernel.inverse_haar(tmp_haar.data(), tmp_haar2.data(), ny2);

            //Copy back, pay attention not to mess up the shape
            for (int i=0; i< nx2; i++) {
                if ( i <  dof_handler.pwidth() )
                    checkCudaErrors(cudaMemcpyAsync(&(this->z_d.data()[i* dof_handler.pheight()]), &(tmp_haar.data()[i* ny2]),
                            dof_handler.pwidth()*sizeof(Mdouble), cudaMemcpyDeviceToDevice));
            }
            checkCudaErrors(cudaDeviceSynchronize());
        }

        //Regularization by direct space sparsity
        if ( this->regType == sparse ) {
            this->z_d = this->x_d;
            kernel.soft_threshold(this->z_d.data(),lag2.data(),this->x_d.data(), rho2, reg_strength,   dof_handler.n_dofs());
            //kernel.tv_regularization(this->x_d,this->z_d,lag2,gamma,rho2,inf->ext_width,inf->ext_height,inf->ext_depth);
        }
        //Regularization by Fourier Space L_2 Norm
        if ( this->regType == quadratic ) {
            this->z_d = this->x_d;
            //the solution with smallest L_2 Norm is obtained, this corresponds to the pseudoinverse
            kernel.pseudo_inverse(this->z_d.data(), lag2.data(), rho2, reg_strength, dof_handler.n_dofs());
        }


        //Regularization by Fourier Space L_2 Norm
        if ( this->regType == TV ) {
            this->z_d = this->x_d;
            //checkCudaErrors(cudaMemcpy(this->z_d.data(), this->x_d.data(), inf->n_bytes_per_frame, cudaMemcpyDeviceToDevice));
            if (true)
                kernel.tv_regularization(this->x_d.data(), this->z_d.data(),
                                         lag2.data(), reg_strength, rho2,
                                         dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth());
            else
            {
                Mdouble dt = 1e-1;

                Mdouble Tol = 1e-4 * std::sqrt(tmp2_d.size()) / dt;
                int n_max_iter = 100;
                int n_iter = 0;
                Mdouble err = 2*Tol;

                while (err > Tol)
                {
                    kernel.tv_derivative(tmp2_d.data(), z_d.data(),
                                         this->im_d().data(),
                                         reg_strength,
                                         dof_handler.pwidth(), dof_handler.pheight(), dof_handler.pdepth());

                    tmp_d = rho2 * x_d + lag2;
                    tmp_d = (-rho2) * z_d + tmp_d;
                    tmp_d -= tmp2_d;

                    z_d = dt * tmp_d + z_d;

                    err = tmp_d.l2_norm();
                    n_iter++;
                    if (n_iter > n_max_iter)
                        break;

                }

                std::cout << "n iter in Reg term : " << n_iter << ", final error : " << (err / ( std::sqrt(tmp2_d.size()) / dt )) << std::endl;
            }
        }

    }
};
} // namespace step35 END

#endif // CUDADriver_STEP_35_H
