#ifndef SMRE_PROBLEM_HH
#define SMRE_PROBLEM_HH

#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>

#include <step-35/ADMMParams.h>
#include <step-35/cuda_kernel_wrapper_step-35.cu.h>

//SciPAL
#include <base/PrecisionTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/cublas_Vector.h>


//Boost used for Noise simulation
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>

namespace step35 {




struct ConvolutionBase {

    virtual void vmult(/*dst, src*/) = 0;

    virtual void Tvmult (/*dst, src*/) = 0;

};

struct DiscreteConvolution : public ConvolutionBase {};

struct FFTConvolution : public ConvolutionBase {};



// We put the convolution into a structure such tht we can use as a matrix in iteratuve solvers.
/*

  Params:
   - psf fwhm (depends on operator)
   - image w/ discretized psf
   - dimension

*/
template <typename NumberType, typename BW>
struct Convolution {

    typedef CudaComplex<NumberType> complex;

    typedef  SciPAL::Vector<NumberType, BW> Vector;


    int depth, width, height;
    int cframesize; // image size for r2c images
    bool init; // for fftw

  // FIXME: reactivate  cufftHandle *plan_fft,*iplan_fft;
    //        CUDAFFT<T, dim, cufft_type, gpu_cuda> FFT; //FIXME use that

    // FIXME: reactivate  fftwf_plan plan_cpu_d_forward;
     // FIXME: reactivate fftwf_plan plan_cpu_d_backward;

    SciPAL::Vector<complex, BW> fm1;
    //Point spread function on device
    SciPAL::Vector<complex, BW> psf_fourier_transform_d;

    dealii::Vector<NumberType> psf_h;

    // In cases where we only want to denoise
    // we have to bypass the convolution operation as then the blurring operator is the identity.
    bool is_delta_peak;

    dealii::Vector<NumberType> psf_1d;
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


    Convolution(const ImageDoFHandler & dh, double fwhm)
        :
          depth(dh.depth_p()), width(dh.width_p()), height(dh.height_p()),
          cframesize(/*ext_*/depth*/*ext_*/width*(/*ext_*/height/2+1)/**sizeof(Mcomplex)*/),
          fm1(cframesize), init(false),
          psf_fourier_transform_d(fm1.size()),
          psf_h(width*height),
          is_delta_peak(false),
          dof_handler(dh)
    {

            SciPAL::Vector<NumberType,BW> fpsf_tmp_d(psf_h.size());

            NumberType psf_norm=0;

            for (int z=0; z<depth; z++ ) {
                for (int x=0; x<width; x++) {
                    for (int y=0; y<height; y++) {
                        psf_h(z*height*width+x*height+y) =
                                _psf(NumberType(x), NumberType(y), NumberType(z), fwhm );
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
      // FIXME: reactivate        plan_fft=new cufftHandle();
        // FIXME: reactivate      iplan_fft=new cufftHandle();
#ifdef DOUBLE_PRECISION
            cufftPlan3d(plan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_D2Z);
            cufftPlan3d(iplan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_Z2D);
#else
   // FIXME: reactivate           cufftPlan3d(plan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_R2C);
   // FIXME: reactivate           cufftPlan3d(iplan_fft, /*inf->ext_*/depth, /*inf->ext_*/width, /*inf->ext_*/height,  CUFFT_C2R);
#endif

            if(BW::arch == gpu_cuda)
            {
                //Fourier transformation of psf using cuFFT
        // FIXME: reactivate          cufftHandle plan;
#ifdef DOUBLE_PRECISION
                cufftPlan3d(&plan, ext_depth, ext_width, ext_height, CUFFT_D2Z);
                cufftExecD2Z(plan, fpsf_tmp_d, fpsf_d);
#else
    // FIXME: reactivate              cufftPlan3d(&plan, /*ext_*/depth, /*ext_*/width, /*ext_*/height,  CUFFT_R2C);
      // FIXME: reactivate            cufftExecR2C(plan, fpsf_tmp_d.data(), psf_fourier_transform_d.data());
#endif
                //Cleanup
     // FIXME: reactivate             cufftDestroy(plan);
    // FIXME: reactivate              getLastCudaError("cufft error!\n");
            }
            else
            {
                //Fourier transformation of psf using fftw
    // FIXME: reactivate              fftwf_plan plan =
//                        fftwf_plan_dft_r2c_2d(width, height,
//                                              fpsf_tmp_d.data(),
//                                              reinterpret_cast<fftwf_complex*>(psf_fourier_transform_d.data()),
//                                              FFTW_PATIENT);
//                fftwf_execute(plan);
//                fftwf_destroy_plan(plan);
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
                psf_1d(x+cut_off) = exp(-std::pow( x - x_0, 2 )/
                                        (2*sigma_psf*sigma_psf) );
                psf_norm += psf_1d(x+cut_off);
            }
            psf_1d /= psf_norm;

            std::cout << psf_1d << std::endl;
        }
        else
            is_delta_peak = true;
//#endif


    }



    void vmult(SciPAL::Vector<NumberType, BW> &dst, const SciPAL::Vector<NumberType, BW>& src)
    {
        if (is_delta_peak)
            this->vmult_id(dst, src);
        else
        {
//            if(BW::arch == gpu_cuda )
//                this->vmult_FFT(dst, src);
//            else
                this->vmult_conv(dst, src);
        }
    }

    //    private:
    void vmult_id(SciPAL::Vector<NumberType, BW> &dst, const SciPAL::Vector<NumberType, BW>& src)
    {
        dst = src;
        return;
    }

    void vmult_conv(SciPAL::Vector<NumberType, BW> &dst, const SciPAL::Vector<NumberType, BW>& src)
    {
        // Assert that dst is large enough.
        // We assume that the image space has
        // the same dimension as the range space.
        if (dst.size() == 0)
            dst.reinit(src.size());


        dealii::Vector<NumberType> tmp_dst(src.size()), tmp_src(dst.size());
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

    void vmult_FFT(SciPAL::Vector<NumberType, BW> &dst, const SciPAL::Vector<NumberType, BW>& src)
    {
        // Assert that dst is large enough.
        // We assume that the image space has
        // the same dimension as the range space.
#ifdef USE_FFT_AGAIN
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

#endif

    }

    //@sect5{Function: _psf}
    //@brief Returns a 2D gaussian convolution kernel in direct space
    NumberType _psf(NumberType x,NumberType y,NumberType z, NumberType fwhm) {

        double sigma = std::max(1e-12, fwhm/2.3548);

        //Fixme
        //            if ( params.dim == 2 )
        //            {
        if ( z < 0.5 )//float conditional to distinguish z=0 vs. z=1
        {
            NumberType tmp = exp( -(std::pow(x - width/2., 2) + std::pow(y - height/2., 2) )
                                  /
                              (2.0*std::pow(sigma, 2)));
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





// @sect3{Class: StatInverseProblem}
//
// This class models a linear algebraic, statistical inverse problem $Y = A*x + \epsilon$,
// where $Y$ contains the measured data, $A$ is the ill-posed operator and $\epsilon$ is the noise.
// This class should contain the host-side representation of $x$, $Y$ and $\epsilon$.
// The purpose of this class is to store the problem-specific data, i.e. solution, right-hand side, operator and
// noise field, and to provide a function which evaluates the residual. In this way the details of the
// inverse problem are completely encapsulated and the augmented Lagrangian can use the inverse problem in an abstract manner.
// Where it is actually stored internally is a non-public thing.
template <typename NumberType, typename BW>
/*

  params:
  - conv op
  - measured data
  - std dev noise
  - dimension

*/
class StatInverseProblem : public dealii::Subscriptor {
public:

    static const ParallelArch c = BW::arch;
public:
    //complex type: double2 for double, float2 for float
    //    typedef typename PrecisionTraits<NumberType, c>::ComplexType complex;
    typedef CudaComplex<NumberType> complex;

    typedef  SciPAL::Vector<NumberType, BW> Vc;

    StatInverseProblem(const AugmentedLagrangianParams & prm)
        :
          convolution(prm.image_as_read.dofs, prm.psf_fwhm)
    {
        rhs = prm.image_as_read.data;

        solution.reinit(prm.image_as_read.data.size());
        eps.reinit(prm.image_as_read.data.size());

        // Consider input image as ideal and use it to generate a synthetic testcase.
        if (prm.simulate)
        {
            // Simulate imaging process
            // by convolution with PSF
            this->convolution.vmult(rhs, rhs );

            // Add Gaussian noise.
            //Constant seed to get reproducable tests
            boost::mt19937 rng;
            // FIXME : reactivate
            //Seed rng
//            if ( params.time_seed )
//                rng.seed(time(NULL));
            //Prepare a normal distribution with mean 0 and standard deviation @p prm.std_dev_noise
            boost::normal_distribution<> nd(0.0, prm.std_dev_noise);
            boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

            std::vector<float> generated_noise(rhs.size(), 0.);

            size_t n_elements = rhs.size();
            for (size_t i = 0; i < n_elements; i++)
                generated_noise[i] = var_nor();

            solution = generated_noise;
            rhs += solution;
            solution = rhs;
        }
    }

protected:
    Vc solution, eps, rhs;

    Convolution<NumberType, BW> convolution;

public:
    // Compute the residual $I - \epsilon - A * x$. To be put into an expression.
    void residual(Vc & res,
                  const Vc & im,
                  const Vc & noise,
                  const Vc & est_im)
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

    const  Vc & measured_data () const { return rhs; }

    const Vc & signal_estimator () const { return solution; }

    const Vc & noise () const { return eps; }

    Vc & signal_estimator () { return solution; }

    Vc & noise () { return eps; }

};

// This class mainly adds the functionality for defining the functional the optimum of which yields the solution to the SMRE problem.
// Although the cost functional only HAS an instance of the @p StatInverseProblem we nevertheless use the latter as public subclass.
// As long as we do not plan to extend @p CostFunctional such that it can work on different types of inverse problems this design does not
// imply any restrictions. We could still make the StatInverseProblem a template parameter and still use it then as base class.
// As the Kernels structure contains details of the implementation we make it a non-public base class.
//
// cf. https://en.wikipedia.org/wiki/Augmented_Lagrangian_method
/*

  Params
  - reg strength
  - reg type
  - penalty parameter rho (for || Y - A*x + e ||^2)

*/
template <typename NumberType, typename BW>
class AugmentedLagrangian : public StatInverseProblem<NumberType, BW>, protected step35::Kernels<NumberType, BW::arch> {

public:

    typedef StatInverseProblem<NumberType, BW> Base1;

     typedef typename Base1::Vc Vc;

    AugmentedLagrangian(const AugmentedLagrangianParams & alp)
        :
          Base1(alp),
          params(&alp)
    {
        this->m_lag_mult.reinit(this->solution);
        this->tmp_res.reinit(this->solution);
    }


    // implements all sorts of evaluating $\mathcal{L}(x,\epsilon)$

    // Evaluation of the Euler-Lagrange equation for the minimization w.r.t.
    // to the primal variable, that is the image we want to reconstruct.
    // The mathematical expression is
    // \f{eqnarray}{ \text{dst}
    //                 & = &
    //                        A^T*\left\lbrack \rho \left(Y - \epsilon -A*x\right) + \Upsilon \right\rbrack
    //                        -
    //                        \frac{\delta R}{\delta x}(x_{old}) \,.
    //             \f}
    void argmin_x_rhs(Vc & dst,
                      const Vc & x_old);


    // The derivative with respect to the noise component yields a subdifferential and effectively is a projection onto the
    // intersection of convex sets.
    void argmin_eps_G();


    // @sect5{Function: update_lagrangian}
    // @brief Update the lagrange multiplier
    // $\rho$ is the penalty parameter
    // $ \Upsilon_{r+1} = \Upsilon_r + \rho \left( Y - A*x - \epsilon \right) $
    void update_lag_multiplier();

    const Vc & lag_mult() const { return m_lag_mult; }


    // The second template parameter's sole purpose is to generate a
    // unique id for the smart pointer.
    dealii::SmartPointer<const AugmentedLagrangianParams, const AugmentedLagrangianParams> params;

protected:
    // In addition to the data of the inverse problem the ADMM requires two additional vectors.
    Vc tmp_res, m_lag_mult;

    void dykstra_gauss();

};


template <typename NumberType, typename BW>
void AugmentedLagrangian<NumberType, BW>::argmin_x_rhs(Vc &dst,
                                                   const Vc &x_old)
{
    // the residual
    // $\text{tmp2}_d = Y - \epsilon - A*x_{old}$
    this->residual(this->tmp_res, this->measured_data(), this->noise(), x_old);

    // $\text{tmp}_d = A * \left(
    //                            \rho \left(Y - \epsilon -A*x \right) + \Upsilon_1
    //                    \right)$
    //
    dst = params->penalty * this->tmp_res + m_lag_mult;

    // This is the only place where we actually have to know something about the internal structure of the inverse problem.
    this->convolution.vmult(dst, dst);

    // Add regularization
 const ImageDoFHandler & dof_handler = this->params->image_as_read.dofs;
    // TODO: switch over regularizations
    this->tv_derivative(this->tmp_res.data(), x_old.data(),
                         this->measured_data().data(),
                         params->regularization_strength,
                         dof_handler.width_p(),
                         dof_handler.height_p(), dof_handler.depth_p());

  // WHat's this?  this->tmp2_d *= params.rho2;

    dst -= this->tmp_res;

}



template <typename NumberType, typename BW>
void AugmentedLagrangian<NumberType, BW>::argmin_eps_G()
{
    this->dykstra_gauss();
}




//@sect5{Function: dykstra_gauss}
//@brief Wrapper around the whole projection procedure
//@param rho parameter for the first constraint

template <typename NumberType, typename BW>
void AugmentedLagrangian<NumberType, BW>::dykstra_gauss()
{
    const ImageDoFHandler & dof_handler = this->params->image_as_read.dofs;

    // old dykstra_gauss and dykstra_gauss_global_mem merged.
    // The former was anyway only a wrapper of the latter.
    this->convolution.vmult(tmp_res, this->signal_estimator());

    NumberType mean = 0;

    // $\epsilon = I -  A * x_{k} + \Upsilon_1/\rho$
    this->noise() = NumberType(1./params->penalty) * m_lag_mult + this->measured_data();

    this->noise() -= tmp_res;


    //        //substract mean of noise
//            this->writeable_e_d().push_to(e_h);

//            mean = std::accumulate(e_h.begin(),e_h.end(), 0.)/e_h.size();
//            this->writeable_e_d() += -mean;


    static const int n_scales = 10;
    static const int min_scale = 0;

    typedef SciPAL::Vector<NumberType, BW> Vc;

    std::vector<NumberType> zero_init(this->noise().size());

    // TODO: make attributes
    Vc h_iter (this->noise().size()); h_iter = zero_init;

    Vc h_old (this->noise().size()); h_old = zero_init;

    Vc h_init (this->noise().size()); h_init = zero_init;

    Vc Q_full (this->noise().size()); Q_full = zero_init;
    Vc Q_M (this->noise().size()); Q_M = zero_init;

    // $h_0 = h$;
    h_init = this->noise();
    h_old = h_init;

    NumberType cs_weight =  //9.46624e-07; // 1024x1024, 0.9, c_S
            9.51511e-07; // 1024x1024, 0.9, c_Ss
    //9.53509e-07; // 1024x1024, 0.1, c_Ss

    int iter = 0;
    bool iterate = true;

    // FIXME: arch flag as template parameter in Kernels structure

    NumberType norm = 2*this->params->dykstra_Tol;  //norm of residuals

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
        h_iter *= this->params->std_dev_noise/(L2_norm * std::sqrt(cs_weight) );
#ifdef ET_AVAILABLE

        Q_full = h_iter - h_old;
#else
        Q_full = h_iter;
        Q_full -= h_old;
#endif

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
                                                  this->params->std_dev_noise,
                                                  dof_handler.width_p(),
                                                  dof_handler.height_p(),
                                                  dof_handler.depth_p()//,
//                                                      n_max_dykstra_steps, dykstra_Tol
                                                  );
        }
        else
        {
            this->dyadic_dykstra_fine_scale_part_cpu(h_iter.data(), h_old.data(),
                                                      Q_full.data(),
                                                     this->params->std_dev_noise,
                                                      dof_handler.width_p(),
                                                      dof_handler.height_p(),
                                                      dof_handler.depth_p()//,
                                                      );
        }
        time += timer.wall_time();


        //Convergence control.
        {
            h_init -= h_iter;
            norm = h_init.l2_norm();

            iterate = ( norm > this->params->dykstra_Tol &&
                        iter < this->params->n_max_dykstra_steps);

            iterate ? h_init = h_iter : this->noise() = h_iter;
            iter++;
        }

    }// iterate end

    std::cout<<"    n Dykstra steps used : " << //iter <<
               ", norm of noise increment : " //<< norm
            <<" mean value of noise : "<<mean <<std::endl;
    std::cout<<"cumulative dykstra time for fine scale sweeps : "<< time<<std::endl;
}






template <typename NumberType, typename BW>
void AugmentedLagrangian<NumberType, BW>::update_lag_multiplier()
{
    // TODO: access to problem data
    this->residual(tmp_res, this->measured_data(), this->noise(), this->signal_estimator() );
    this->m_lag_mult =  params->penalty * tmp_res +  this->m_lag_mult;
}


//  Vc& Vc::operator = <BX<BX<L*Vc >, +, Vc > >(SciPAL::Expr<BX<BX<L*Vc >, +, Vc > > const&)



// @sect3{Class: ADMMStepper}
//
// The overall problem of reconstructing an image using the SMRE approach is put into a class
// which administrates memory management, I/O and running the solver.
// From a structural point of view it is very similar to the examples of stationary PDEs in deal.II.
// To some extent this is no spurprise as mathematically we try to solve an inverse problem
// which basically follows the same pattern of solving a PDE. The only difference is that applying the
// operator representing the inverse problem usually means to compute a convolution or evaluating an
// integral equation rather than just computing derivatives.
// Much of the things in this class can be generalized to fit other use cases as well adn could be factored
// out into a base class InverseProblem. However, for the time being we do not do this and reserve this for a later step.
template <typename NumberType, typename BW>
class ADMMStepper {

public:
    ADMMStepper(const ADMMParams & prm)
        :
          params (&prm)
    {}

    typedef SciPAL::Vector<NumberType, BW> Vc;


    // one step the ADMM
    void step(AugmentedLagrangian<NumberType, BW> & sip);

protected:


    // This member function runs Heun's method to find a new image reconstruction
    // for a given noise.
    void x_step_adaptive(AugmentedLagrangian<NumberType, BW> & sip);


    // dump intermediate images and convergence data
    void output();


    // The second template parameter's sole purpose is to generate a
    // unique id for the smart pointer.
    dealii::SmartPointer<const ADMMParams, const ADMMParams> params;




    NumberType dt;


    Vc  x_tmp, x_old;

    Vc tmp2_d;

    Vc  k_1;
    Vc  k_2;

    Vc  sol_diff;

};

}



// @sect3{Function: setup_system}
//
// We allocate the memory for the various solution vectors (i.e. images) and set the initial condition.
template <typename NumberType, typename BW>
void step35::ADMMStepper<NumberType, BW>::x_step_adaptive(AugmentedLagrangian<NumberType, BW> & sip)
{

    NumberType err = 2*params->heun_Tol;

    int n_err_iter = 0;
    // For a given Lagrange parameter @p lag1 and a given estimate of the noise we compute
    // a new estimate for how the reconstructed image should look like.
    while (err > params->heun_Tol)
    {
        // Save old solution.
        this->x_old = sip.signal_estimator();

        NumberType Delta = 2*params->dt_adaption_Tol;

        // Evaluate rhs of Euler-Lagrange for explicit Euler step.
        // This does not change over the course of the timestep adaption.
        // Hence, we have to compute it only once.
        sip.argmin_x_rhs(k_1, this->x_old);


        int dt_iter = 1;
        NumberType delta_t;
        // The inner loop does the time-step adaption.
        for (/*dt_iter*/; dt_iter <= params->n_max_heun_dt_adaption_steps; dt_iter++)
        {
            delta_t = params->pseudo_timestep;

            // Euler step for intermediate solution.
            x_tmp = delta_t * k_1 + x_old;

            // Evaluate rhs of Euler-Lagrange
            sip.argmin_x_rhs(k_2, x_tmp);

            // The difference of the two solutions is given by the differnce of the derivatives
            // and serves the timestep adaption.
            sol_diff = k_2;
            sol_diff -= k_1;
            sol_diff *= 0.5*delta_t;

            // Heun step
            k_2 += k_1;
            sip.signal_estimator() = NumberType(0.5*delta_t) * k_2 + x_old;


            Delta = sol_diff.l2_norm()/std::sqrt(sol_diff.size());

            NumberType new_dt = 0.9 * params->dt_adaption_Tol / Delta;

            // To avoid excessive oscillations we limit the increase of the timestep.
            if (Delta < params->dt_adaption_Tol) {
                if (new_dt > 2*delta_t)
                    new_dt = 2*delta_t;
            }
            else
                new_dt = 0.707 * delta_t;

             delta_t = const_cast<ADMMParams&>(*params).pseudo_timestep = new_dt;

            if (Delta <= params->dt_adaption_Tol)
                break;
        }

#ifdef DEBUG
        std::cout << "             n dt adaptions : " << dt_iter << ", new timestep : " << params->pseudo_timestep << ", Delta : " << Delta << std::endl;
#endif

        // Convergence control for the outer iteration.
        tmp2_d = sip.signal_estimator();
        tmp2_d -= this->x_old;
        err = tmp2_d.l2_norm()/std::sqrt(tmp2_d.size());
        n_err_iter++;
        if (n_err_iter >= params->n_max_heun_steps)
            break;
    }
    std::cout << "    n Heun steps used : " << n_err_iter <<
              ", norm of solution increment : "<< err<< std::endl;
}


template <typename NumberType, typename BW>
void step35::ADMMStepper<NumberType, BW>::step(step35::AugmentedLagrangian<NumberType, BW> &sip)
{


    //Argmin w.r.t. x
 // /* still in CUDADriver : */ // driver.
    this->x_step_adaptive(sip);


    //Argmin w.r.t. e, is equivalent to a projection
    // if ( params.gnoise > 0 )
        /* still in CUDADriver*/ //dykstra_gauss(params.rho1);
    sip.argmin_eps_G();

    //Update the Lagrange Multipliers
    /* still in CUDADriver*/ //update_lagrangian(1./params.inv_gamma, params.alpha2);

    sip.update_lag_multiplier();


}







#endif // SMRE_PROBLEM_HH

