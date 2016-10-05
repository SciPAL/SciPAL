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

Copyright  Lutz Künneke and Jan Lebert 2014-2015, Stephan Kramer, Johannes Hagemann 2014-2016
*/

#ifndef CUDADriver_STEP_35_H
#define CUDADriver_STEP_35_H


//std
#include<iostream>
#include <cmath>

//CUDA
#include <cuda_runtime.h>

//CUDA cufft
#include <cufft.h>
#include <fftw3.h>

//shift for fft
#include <cufftShiftInterface.h>

//SciPAL
#include <base/PrecisionTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/cublas_Vector.h>

//Our stuff
#include <step-35/cuda_kernel_wrapper_step-35.cu.h>
#include <step-35/cuda_helper.h>
#include <step-35/ADMMParams.h>
#include <step-35/smre_problem.hh>

// deal.II
#include <deal.II/base/timer.h>


namespace step35 {
// We encapsulate each project into a dedicated namespace
// in order to be able to re-use parts of a test program in others.

// @sect4{Class: CUDADriver}
//
// Main class, sets up an instance of queue and info class. A number of threads is started which
// each handles a cuda stream and processes items obtained by the queue. The main thread adds items
// to the queue until everything is processed.
template<typename NumberType, typename BW >
class CUDADriver : protected step35::Kernels<NumberType, BW::arch> {

public:
    //@sect5{Function: dykstra_gauss}
    //@brief Wrapper around the whole projection procedure
    //@param rho parameter for the first constraint
    void dykstra_gauss(const NumberType rho)
    {
        // old dykstra_gauss and dykstra_gauss_global_mem merged.
        // The former was anyway only a wrapper of the latter.
        convolution.vmult(tmp_d, this->x_d);

        NumberType mean = 0;
        // $\epsilon = I -  A * x_{k} + \Upsilon_1/\rho$
        this->writeable_e_d() = Mdouble(1./rho) * lag1 + this->im_d();
        this->writeable_e_d() -= tmp_d;

//        //substract mean of noise
//            this->writeable_e_d().push_to(e_h);

//            mean = std::accumulate(e_h.begin(),e_h.end(), 0.)/e_h.size();
//            this->writeable_e_d() += -mean;


        static const int n_scales = 10;
        static const int min_scale = 0;

        typedef SciPAL::Vector<NumberType, BW> Vc;

        std::vector<NumberType> zero_init(this->e_d().size());

        Vc h_iter (this->e_d().size()); h_iter = zero_init;

        Vc h_old (this->e_d().size()); h_old = zero_init;

        Vc h_init (this->e_d().size()); h_init = zero_init;

        Vc Q_full (this->e_d().size()); Q_full = zero_init;
        Vc Q_M (this->e_d().size()); Q_M = zero_init;

        // $h_0 = h$;
        h_init = this->e_d();
        h_old = h_init;

        NumberType cs_weight =  //9.46624e-07; // 1024x1024, 0.9, c_S
                9.51511e-07; // 1024x1024, 0.9, c_Ss
        //9.53509e-07; // 1024x1024, 0.1, c_Ss

        int iter = 0;
        bool iterate = true;

        // FIXME: arch flag as template parameter in Kernels structure
        step35::Kernels<NumberType, BW::arch> kernel;
        NumberType norm = 2*dykstra_Tol;  //norm of residuals

        NumberType time = 0;

        while (iterate)
        {

            // Start with estimator for whole image
            h_old -= Q_full;
            h_init = h_old;
            h_iter = h_old;
            // assume single subset for a moment
            NumberType L2_norm = h_old.l2_norm(); // i.e. sqrt(sum_of_squares)

            // project
            h_iter *= this->sigma_noise/(L2_norm * std::sqrt(cs_weight) );

            Q_full = h_iter - h_old;
            h_old = h_iter;

            // Q_0 <- Q_full
            // initialize : d.h_init with global h_init


            dealii::Timer timer;
            timer.restart();
            // loop over subsets
            if(typeid(BW) == typeid(cublas))
            {
                this->dyadic_dykstra_fine_scale_part(h_iter.data(), h_old.data(),
                                                      Q_full.data(),
                                                      this->sigma_noise,
                                                      dof_handler.pwidth(),
                                                      dof_handler.pheight(),
                                                      dof_handler.pdepth()//,
//                                                      n_max_dykstra_steps, dykstra_Tol
                                                      );
            }
            else
            {
                this->dyadic_dykstra_fine_scale_part_cpu(h_iter.data(), h_old.data(),
                                                          Q_full.data(),
                                                          this->sigma_noise,
                                                          dof_handler.pwidth(),
                                                          dof_handler.pheight(),
                                                          dof_handler.pdepth()//,
                                                          );
            }
            time += timer.wall_time();



            //Convergence control.            {
                h_init -= h_iter;
                norm = h_init.l2_norm();

                iterate = ( norm > dykstra_Tol &&
                            iter < n_max_dykstra_steps);

                iterate ? h_init = h_iter : this->writeable_e_d() = h_iter;
                iter++;
            }
        }// iterate end

        std::cout<<"    n Dykstra steps used : " << iter <<
                   ", norm of noise increment : "<< norm <<" mean value of noise : "<<mean <<std::endl;
        std::cout<<"cumulative dykstra time for fine scale sweeps : "<< time<<std::endl;
    }


};
template<typename Mdouble, typename BW >
class CUDADriverMurks {

    static const ParallelArch c = BW::arch;
public:
    //complex type: double2 for double, float2 for float
    //    typedef typename PrecisionTraits<Mdouble, c>::ComplexType complex;
    typedef CudaComplex<Mdouble> complex;

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

        dealii::Vector<Mdouble> psf_1d;
        int cut_off;
        ImageDoFHandler dof_handler;

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
              dof_handler(dealii::Point<3>(width, height, depth))
        {

                // TO DO: compute psf_h here

                SciPAL::Vector<Mdouble,BW> fpsf_tmp_d(psf_h.size());

                Mdouble psf_norm=0;

                for (int z=0; z<depth; z++ ) {
                    for (int x=0; x<width; x++) {
                        for (int y=0; y<height; y++) {
                            psf_h(z*height*width+x*height+y) =
                                    _psf(Mdouble(x), Mdouble(y), Mdouble(z), fwhm );
                            psf_norm += psf_h[z*height*width+x*height+y];
                        }
                    }
                }
                std::cout<< "psf_norm: "<<psf_norm << std::endl;
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

                if(BW::arch == gpu_cuda)
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
                }
                else
                {
                    //Fourier transformation of psf using fftw
                    fftwf_plan plan =
                            fftwf_plan_dft_r2c_2d(width, height,
                                                  fpsf_tmp_d.data(),
                                                  reinterpret_cast<fftwf_complex*>(psf_fourier_transform_d.data()),
                                                  FFTW_PATIENT);
                    fftwf_execute(plan);
                    fftwf_destroy_plan(plan);
                }



//#ifdef USE_FFT_CONV

            // Compute 1D PSF
            if(fwhm > 1e-10)
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
//#endif


        }



        void vmult(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            if (is_delta_peak)
                this->vmult_id(dst, src);
            else
            {
                if(BW::arch == gpu_cuda )
                    this->vmult_FFT(dst, src);
                else
                    this->vmult_conv(dst, src);
            }
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


            dealii::Vector<Mdouble> tmp_dst(src.size()), tmp_src(dst.size());
            src.push_to(tmp_dst);

            // Number of rows before and after the current row
            // which have to be taken into account.
            int N = cut_off;
            for (int k=0; k< depth; k++)
#pragma omp parallel for
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
#pragma omp parallel for
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

        }

        void vmult_FFT(SciPAL::Vector<Mdouble, BW> &dst, const SciPAL::Vector<Mdouble, BW>& src)
        {
            // Assert that dst is large enough.
            // We assume that the image space has
            // the same dimension as the range space.

            // setup fftw on first run
            //fftw plans
            if(BW::arch == cpu)
            {
            if(init == false)
            {
                init = true;
//                fftwf_plan_with_nthreads();
                this->plan_cpu_d_forward
                        = fftwf_plan_dft_r2c_2d(width, height,
                                                (src.data()),
                                                reinterpret_cast<fftwf_complex*>(fm1.data()),
                                                FFTW_ESTIMATE);

                this->plan_cpu_d_backward
                        = fftwf_plan_dft_c2r_2d(width, height,
                                                reinterpret_cast<fftwf_complex*>(fm1.data()),
                                                (dst.data()),
                                                FFTW_ESTIMATE);
            }
            }
            if (dst.size() == 0)
                dst.reinit(src.size());

            //            cufftShift_2D_impl(src.data(), width, height);
            if(BW::arch == gpu_cuda) {
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
                dst *= 1./dst.size();
                //ugly...
                kernelConf conf;
                int threadsPerBlock_X = 16;
                int threadsPerBlock_Y = 16;
                conf.block = dim3(threadsPerBlock_X, threadsPerBlock_Y, 1);
                conf.grid = dim3(((width/2) / threadsPerBlock_X), ((width/2) / threadsPerBlock_Y), 1);
                cufftShift_2D_config_impl(dst.data(), width, height, &conf);
            }
            else
            {
                fftwf_execute(plan_cpu_d_backward);
            }

//            dst *= 1./dst.size();
            //ugly...
//            kernelConf conf;
//            int threadsPerBlock_X = 16;
//            int threadsPerBlock_Y = 16;
//            conf.block = dim3(threadsPerBlock_X, threadsPerBlock_Y, 1);
//            conf.grid = dim3(((width/2) / threadsPerBlock_X), ((width/2) / threadsPerBlock_Y), 1);
//            cufftShift_2D_config_impl(dst.data(), width, height, &conf);



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
                                  (2.0*boost::math::pow<2>(sigma)));
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

        // A * x + eps
        res += noise;
        // A * x + eps -I
        res -= im;

        // I - eps - A * x
        res *= -1;
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
    CUDADriverMurks(std::vector<Mdouble> & cs_h,
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
    ~CUDADriverMurks () {



#ifdef TIME_KERNELS
        //Write total dykstra kernel execution time at shutdown to file "kernel_time"
        if (kernel_time)
            out_time << kernel_time << std::endl;
        out_time.close();
#endif

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
    //       & = &
    //        A^T*\left\lbrack     \rho_1 \left(I-e-A*x\right) +\Upsilon_1      \right\rbrack
    //         -
    //        \rho_2 \frac{\delta J}{\delta x}(x_{old}) \,.
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
        dst = rho1 * this->tmp2_d + lag_1;

        convolution.vmult(dst, dst);

        // Add regularization
        step35::Kernels<Mdouble, BW::arch> kernel;
//        kernel.L1_derivative(this->tmp2_d.data(),
//                             x_old.data(),
//                             reg_strength,
//                             dof_handler.pwidth(),
//                             dof_handler.pheight(), dof_handler.pdepth());
        kernel.tv_derivative(this->tmp2_d.data(), x_old.data(),
                             this->im_d().data(),
                             reg_strength,
                             dof_handler.pwidth(),
                             dof_handler.pheight(), dof_handler.pdepth());

        this->tmp2_d *= params.rho2;

        dst -= this->tmp2_d;
    }


    void x_step_adaptive(ADMMParams & params)
    {

        const Mdouble rho1 = params.rho1;

        const Mdouble reg_strength = params.reg_strength;


        SciPAL::Vector<Mdouble, BW> & x_tmp = this->tmp_lag2;

        SciPAL::Vector<Mdouble, BW> & k_1 = this->tmp_haar;
        SciPAL::Vector<Mdouble, BW> & k_2 = this->tmp_haar2;

        SciPAL::Vector<Mdouble, BW> & sol_diff = this->tmp_lagr;

        Mdouble err = 2*params.heun_Tol;

        int n_err_iter = 0;
        // For a given Lagrange parameter @p lag1 and a given estimate of the noise we compute
        // a new estimate for how the reconstructed image should look like.
        while (err > params.heun_Tol)
        {
            // Save old solution.
            this->x_old = this->x_d;

            Mdouble Delta = 2*params.dt_adaption_Tol;

            // Evaluate rhs of Euler-Lagrange for explicit Euler step.
            // This does not change over the course of the timestep adaption.
            // Hence, we have to compute it only once.
            argmin_x_rhs(k_1, this->x_old, this->lag1, params);


            int dt_iter = 1;
            Mdouble delta_t;
            // The inner loop does the time-step adaption.
            for (/*dt_iter*/; dt_iter <= params.n_max_heun_dt_adaption_steps; dt_iter++)
            {
                delta_t = 1./params.inv_gamma;

                // Euler step for intermediate solution.
                x_tmp = delta_t * k_1 + x_old;

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

                Mdouble new_dt = 0.9 * params.dt_adaption_Tol / Delta;

                // To avoid excessive oscillations we limit the increase of the timestep.
                if (Delta < params.dt_adaption_Tol) {
                    if (new_dt > 2*delta_t)
                        new_dt = 2*delta_t;
                }
                else
                    new_dt = 0.707 * delta_t;


                delta_t = params.inv_gamma = 1./(new_dt);

                if (Delta <= params.dt_adaption_Tol)
                    break;
            }
            //  std::cout << "             n dt adaptions : " << dt_iter << ", new timestep : " << 1./params.inv_gamma << ", Delta : " << Delta << std::endl;

            // Convergence control for the outer iteration.
            tmp2_d = this->x_d;
            tmp2_d -= this->x_old;
            err = tmp2_d.l2_norm()/std::sqrt(tmp2_d.size());
            n_err_iter++;
            if (n_err_iter >= params.n_max_heun_steps)
                break;
        }
        std::cout << "    n Heun steps used : " << n_err_iter <<
                  ", norm of solution increment : "<< err<< std::endl;


        this->z_d = this->x_d;
    }




    //@sect5{Function: dykstra_gauss}
    //@brief Wrapper around the whole projection procedure
    //@param rho parameter for the first constraint
    void dykstra_gauss(const Mdouble rho)
    {
        // old dykstra_gauss and dykstra_gauss_global_mem merged.
        // The former was anyway only a wrapper of the latter.
        convolution.vmult(tmp_d, this->x_d);

        Mdouble mean = 0;
        // $\epsilon = I -  A * x_{k} + \Upsilon_1/\rho$
        this->writeable_e_d() = Mdouble(1./rho) * lag1 + this->im_d();
        this->writeable_e_d() -= tmp_d;

//        //substract mean of noise
//            this->writeable_e_d().push_to(e_h);

//            mean = std::accumulate(e_h.begin(),e_h.end(), 0.)/e_h.size();
//            this->writeable_e_d() += -mean;


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

        // FIXME: arch flag as template parameter in Kernels structure
        step35::Kernels<Mdouble, BW::arch> kernel;
        Mdouble norm = 2*dykstra_Tol;  //norm of residuals

        Mdouble time = 0;

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

            dealii::Timer timer;
            timer.restart();
            // loop over subsets
            if(typeid(BW) == typeid(cublas))
            {
                kernel.dyadic_dykstra_fine_scale_part(h_iter.data(), h_old.data(),
                                                      Q_full.data(),
                                                      this->sigma_noise,
                                                      dof_handler.pwidth(),
                                                      dof_handler.pheight(),
                                                      dof_handler.pdepth()//,
//                                                      n_max_dykstra_steps, dykstra_Tol
                                                      );
            }
            else
            {
                kernel.dyadic_dykstra_fine_scale_part_cpu(h_iter.data(), h_old.data(),
                                                          Q_full.data(),
                                                          this->sigma_noise,
                                                          dof_handler.pwidth(),
                                                          dof_handler.pheight(),
                                                          dof_handler.pdepth()//,
                                                          );
            }
            time += timer.wall_time();

            //Convergence control.
            {
                h_init -= h_iter;
                norm = h_init.l2_norm();

                iterate = ( norm > dykstra_Tol &&
                            iter < n_max_dykstra_steps);

                iterate ? h_init = h_iter : this->writeable_e_d() = h_iter;
                iter++;
            }
        }// iterate end
        std::cout<<"    n Dykstra steps used : " << iter <<
                   ", norm of noise increment : "<< norm <<" mean value of noise : "<<mean <<std::endl;
        std::cout<<"cumulative dykstra time for fine scale sweeps : "<< time<<std::endl;
    }



};
} // namespace step35 END

#endif // CUDADriver_STEP_35_H
