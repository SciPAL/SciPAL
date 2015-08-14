//@sect3{File: step-35.hh}
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

Copyright  Lutz Künneke, Jan Lebert 2014
*/

#ifndef STEP35_HH
#define STEP35_HH

//std
#include<string>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include <time.h>


//Boost used for Noise simulation
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/pow.hpp>

//TIFF used for IO
#include<tiff.h>
#include<tiffio.h>

//Our stuff
#include "extremeValueStatisticsGenerator.h"
#include "ADMMParams.h"
#include <cuda_driver_step-35.h>
#include <cuda_driver_step-35.hh>

// SciPal includes
//
// This is an auxiliary structure which collects data related to the GPUs
// available in the computer on which this program gets executed.
// For details about the cuda<something> functions have a look at
// the CUDA reference manual.
#include <base/GPUInfo.h>

//DEBUG
#include <chrono>

namespace step35 {

//@sect4{Class: ADMM}
//
// To make this test facility extendible, we implement
// a class for a simple user interface. Its primary tasks are
// - management of run-time parameters by a simple text-based parameter file
// - setting device parameters according to the user's parameters
// - preprocessing and output of results
template<typename T>
class ADMM {

public:

    ADMM(int argc, char *argv[], SciPAL::GPUInfo &g);
    void run();

    T _psf(T x, T y, T z, T sigma);
    void create_psf();
    void read_image(std::string path, double gnoise);
    template<typename Q>
    void write_image(std::string path, std::vector<Q> & in, int mheight, int mwidth, int mdepth,double gnoise);
    template<typename Mpatch>
    void add_gaussian_noise(step35::CUDADriver<Mpatch, T, gpu_cuda> &driver);

private:
    SciPAL::GPUInfo & gpuinfo;
    //TIFF object used for IO
    TIFF *tif;
    //TIFF options, will be saved to make the output identical to the input
    int sampleperpixel, width, height, photom, bps, depth;
    //Input buffer for uint16 tiff image
    uint16* buf_in;
    //pwidth: width after convolution with point spread function
    //pwidth0: width after padding on next multiple of 32
    int pwidth, pwidth0;
    //pheight: height after convolution with point spread function
    //pheight0: height after padding on next multiple of 32
    int pheight, pheight0;
    //pdepth: depth after convolution with psf
    //pdepth0: depth after padding on next multiple of 32
    int pdepth, pdepth0;
    //Array holding the signal
    // FIXME: store images in QImage. Applies to other attributes as well.
    std::vector<T> input_image;
    //Array holding the signal as read from the tif
    std::vector<T> image_as_read;
    //Array holding the point spread function
    std::vector<T> psf;

    //DEBUG
    std::chrono::high_resolution_clock::time_point clock1;
    std::chrono::high_resolution_clock::time_point clock2;

    template<typename field_patch, typename driver_patch>
    void __run(extremeValueStatisticsGenerator<field_patch, T, gpu_cuda> &field, step35::CUDADriver<driver_patch, T, gpu_cuda> &driver);

protected:
    ADMMParams<T> params;
};
}

// @sect5{Constructor: ADMM}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
template<typename T>
step35::ADMM<T>::ADMM(int argc, char *argv[], SciPAL::GPUInfo &g)
    : gpuinfo(g) {
    //DEBUG
    clock1 = std::chrono::high_resolution_clock::now();

    // Declare and read parameters from a file. Basically, the parameter
    // file must have the same name as the binary. The extension has
    // to be ".prm". What has been read will be dumped into a log file.
    dealii::ParameterHandler prm_handler;

    ADMMParams<T>::declare(prm_handler);

    // Get the current working directory ...
    QDir cwd = QDir::current();

    // i.e. where the program has been started.
    QDir launch_dir = cwd;

    // By default, the parameter file has the same name as the binary
    std::string prm_filename;
    if (argc == 1) {
        prm_filename  = argv[0];
        prm_filename += ".prm";

        QFileInfo tmp(prm_filename.c_str());
        if (!tmp.exists()) {
            std::cerr << "No parameter file found! Creating default one as " << prm_filename
                      << " and abort!" << std::endl;
            std::ofstream create_default(prm_filename);
            prm_handler.print_parameters (create_default,
                                          dealii::ParameterHandler::Text);
            create_default.close();
            std::abort();
        }
    } else {
        // Whatever gets passed as first command line argument is considered as path
        // to a parameter file.

        // We convert the sequence of characters into something more meaningful.
        QFileInfo tmp(argv[1]);

        std::cout << "Given parameter file : " << argv[1] << std::endl;


        // Before we proceed, let us figure out whether the given parameter file exists.
        // Note: If the file is a symlink that points to a non existing file,
        // false is returned as well.
        if(!tmp.exists()) {
            std::cerr << "The following parameter file does not exist:\n"
                      << argv[1] << std::endl;

            qFatal("Cannot proceed without proper path to parameter file");
        }

        // Next, we subdivide the given filename into its path and filename
        // so that the corresponding subdirectories can be created.
        QString prm_path = tmp.absolutePath();
        cwd.setPath(prm_path);
        cwd.makeAbsolute();
        prm_filename = tmp.fileName().toStdString();

        std::cout << "Parameter file path : "
                  << tmp.absolutePath().toStdString().c_str()
                  << std::endl;
    }

    std::cout << "Parameter file : " << prm_filename  << std::endl;

    prm_handler.read_input(argv[1]); // (prm_path+ QDir::separator() ).toStdString() + prm_filename);

    QDir::setCurrent(launch_dir.absolutePath());

    this->params.get(prm_handler);

    //Create log file

    prm_filename += ".log";
    std::ofstream log_out_text(("./" + QString(prm_filename.c_str()).split("/").last()).toStdString().c_str());
    prm_handler.print_parameters (log_out_text,
                                  dealii::ParameterHandler::Text);
}

//@sect5{Function: create_psf}
//@brief Creates point spread function
template<typename T>
void step35::ADMM<T>::create_psf() {
    // psf=new T[pwidth*pheight*pdepth]; // there is no corresponding call to delete in this function. Hence this may be the origin of a memory leak.
    this->psf.resize(pwidth*pheight*pdepth, T(0));
    T psf_norm=0;
     // FIXME: use scipal::Vector and write psf_norm = psf.l1_norm or psf.l2_norm, std::accumulate is also possible
    for (int z=0; z<pdepth; z++ ) {
        for (int x=0; x<pwidth; x++) {
            for (int y=0; y<pheight; y++) {
                psf[z*pheight*pwidth+x*pheight+y] = _psf((T)x,(T)y,(T)z,params.sigmaf);
                psf_norm += psf[z*pheight*pwidth+x*pheight+y];
            }
        }
    }
    // FIXME: use scipal::Vector and write psf *= 1./psf_norm
    for (int x=0; x<pwidth; x++) {
        for (int y=0; y<pheight; y++) {
            for (int z=0; z<pdepth; z++ ) {
                psf[z*pheight*pwidth+x*pheight+y] = psf[z*pheight*pwidth+x*pheight+y]/psf_norm;
            }
        }
    }
}

//@sect5{Function: _psf}
//@brief Returns a 2D gaussian convolution kernel in direct space
template<typename T>
T step35::ADMM<T>::_psf(T x,T y,T z,T sigma) {

    if ( params.dim == 2 )
    {
        if ( z < 0.5 )//float conditional to distinguish z=0 vs. z=1
        {
            if ( abs(x-sigma) <= std::ceil(sigma )&& abs(y-sigma) <= std::ceil(sigma))
            {
                return exp(-(boost::math::pow<2>(x - sigma) + boost::math::pow<2>(y - sigma))/
                           (TWO*boost::math::pow<2>(HALF*sigma)));
            }
            else
                return 0;
        }
        else
            return 0;
    }
    else
    {
        if ( abs(x-sigma) <= sigma && abs(y-sigma) <= sigma && abs(z-sigma) <= sigma )
        {
            return exp(-(boost::math::pow<2>(x - sigma) + boost::math::pow<2>(y - sigma) +
                         boost::math::pow<2>(z - sigma))/(TWO*boost::math::pow<2>(HALF*sigma)));
        }
        else
            return 0;
    }
}

//@sect5{Function: read_image}
//@brief Read tiff directory as 2D or 3D image
//@param path path to the tiff directory
//@param gnoise standard deviation of Gaussian noise
template<typename T>
void step35::ADMM<T>::read_image(std::string path, double gnoise) {
    std::cout << "Reading in image : " << path.c_str() << std::endl;
    tif=TIFFOpen(path.c_str(), "r");
    TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&sampleperpixel);
    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);
    TIFFGetField(tif,TIFFTAG_PHOTOMETRIC,&photom);
    TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bps);
    //get the depth
    depth=0;
    do
    {
        depth++;
    } while ( TIFFReadDirectory(tif) );
    pdepth=depth;
    pdepth0=depth;
    pwidth=width;
    pheight=height;
    pwidth0=pwidth;
    pheight0=pheight;
    //We use a Haar Wavelet implementation which works on quadratic images
    pwidth=std::max(pwidth,pheight);
    pheight=std::max(pwidth,pheight);
    //reinit the tif
    tif=TIFFOpen(path.c_str(), "r");
    buf_in=(uint16*)_TIFFmalloc(sizeof(uint16)*width);

    // FIXME: use math-aware Vectors not containers
    input_image.resize(pwidth*pheight*pdepth, T(0));
    image_as_read.resize(input_image.size(), T(0));

    int dd=0;
    do
    {
        for (int y=0; y<height; y++) {
            if(TIFFReadScanline(tif,buf_in,y,0) < 0) {
                std::cerr << "reading error\n";
                return;
            }
            for (int x=0; x<width; x++) {
                if (params.anscombe) {
                    //Do the anscombe transformation
                    input_image[dd*pwidth*pheight+x*pheight+y]=sqrt(buf_in[x]+3.0/8.0)*gnoise;
                }
                else {
                    input_image[dd*pwidth*pheight+x*pheight+y]=buf_in[x];
                }
#ifdef READ_IMG_DEBUG
                if (x%465==0)
                    printf("I(%d, %d, %d) : %f, buf: %d", x, y, dd, image[dd*pwidth*pheight+x*pheight+y], buf_in[x]);
#endif
            }
        }
        std::cout << "Reading tiff layer " << dd << std::endl;
        dd++;
    } while ( TIFFReadDirectory(tif) );
}


//@sect5{Function: write_image}
//@brief write image to a .tif file
//@param path path to image
//@param in what to read with dimensions pheight, pwidth
//@param mheight output height
//@param mwidth output width
//@param mdepth number of layers in image stack
//@param gnoise standard deviation of Gaussian noise
template<typename T> template<typename Q>
void step35::ADMM<T>::write_image(std::string path, std::vector<Q> &in, int mheight, int mwidth, int mdepth, double gnoise) {
    uint16* pbuf=(uint16*)_TIFFmalloc(sizeof(uint16)*mwidth);
    TIFF *imout=TIFFOpen(path.c_str(),"w");
    bps=16;
    TIFFSetField(imout, TIFFTAG_ROWSPERSTRIP,TIFFDefaultStripSize(tif, mwidth*sampleperpixel));
    TIFFSetField(imout, TIFFTAG_IMAGEWIDTH, mwidth);  // set the width of the image
    TIFFSetField(imout, TIFFTAG_IMAGELENGTH, mheight);    // set the height of the image
    TIFFSetField(imout, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
    TIFFSetField(imout, TIFFTAG_BITSPERSAMPLE, bps);    // set the size of the channels
    TIFFSetField(imout, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image->
    TIFFSetField(imout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, photom);
    TIFFSetField(imout, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    for (int z=0; z<mdepth; z++) {
        for (int y=0; y<mheight; y++) {
            for (int x=0; x<mwidth; x++) {
                if ( in[z*pheight*pwidth+x*pheight+y] > 0 )
                    if (params.anscombe) {
                        // Inverse anscombe transformation
                        pbuf[x]=(uint16)round((in[z*pheight*pwidth+x*pheight+y]
                                               *in[z*pheight*pwidth+x*pheight+y])/gnoise);
                    }
                    else
                        pbuf[x]=(uint16)round(in[z*pheight*pwidth+x*pheight+y]);
                else
                    pbuf[x]=0;
            }
            if (TIFFWriteScanline(imout, pbuf, y, 1) < 0) {
                std::cerr << "writing error" << std::endl;
                break;
            }
        }
        if ( z < mdepth-1 ) {
            // Prepare next layer
            TIFFWriteDirectory(imout);
            TIFFSetField(imout, TIFFTAG_ROWSPERSTRIP,TIFFDefaultStripSize(tif, mwidth*sampleperpixel));
            TIFFSetField(imout, TIFFTAG_IMAGEWIDTH, mwidth);  // set the width of the image
            TIFFSetField(imout, TIFFTAG_IMAGELENGTH, mheight);    // set the height of the image
            TIFFSetField(imout, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
            TIFFSetField(imout, TIFFTAG_BITSPERSAMPLE, bps);    // set the size of the channels
            TIFFSetField(imout, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image->
            TIFFSetField(imout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, photom);
            TIFFSetField(imout, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        }
    }
    TIFFClose(imout);
    _TIFFfree(pbuf);
}

//@sect5{Function: add_gaussian_noise}
//@brief Simulates dataset by adding gaussian noise, the whole driver is given to used on device convolution
template<typename T> template<typename Mpatch>
void step35::ADMM<T>::add_gaussian_noise (step35::CUDADriver<Mpatch, T, gpu_cuda> &driver) {
    //Constant seed to get reproducable tests
    boost::mt19937 rng;
    //Seed rng
    if ( params.time_seed )
        rng.seed(time(NULL));
    //Prepare a normal distribution with mean 0 and standard deviation params.gnoise
    boost::normal_distribution<> nd(0.0, params.gnoise);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

    driver.conv2(driver.inf->im_d.array().val(), driver.inf->im_d.array().val());
    //Push data from device to host
    driver.get_data();
    //Now add gaussian noise
    for (int z=0; z<pdepth; z++) {
        for (int x=0; x<pwidth; x++) {
            for (int y=0; y<pheight; y++) {
                driver.inf->im_h[z*pwidth*pheight+x*pheight+y] += var_nor();
                //The simulated data is meant to represent a noisy image, which is naturally non-negative everywhere
                if (driver.inf->im_h[z*pwidth*pheight + x*pheight + y] < 0)
                    driver.inf->im_h[z*pwidth*pheight + x*pheight + y] = 0;
            }
        }
    }
    //Push data from host to device
    driver.push_data();

    //Write the noisy image for control reasons
    write_image("noise.tif", driver.inf->im_h, pheight, pwidth, pdepth, params.gnoise);
    std::cout << "noise written" << std::endl;
}

// @sect5{Function: run}
//
// Read in settings and data and perform the algorithm
template<typename T>
void step35::ADMM<T>::run() {
    std::cout << "Starting run" << std::endl;
    //Quantile of the extreme value statistics in (0,1). The larger the smoother the result
    // FIXME: why is this not a runtime parameter accessible from the parameter file?
    const T alpha_quant=0.9;
    //Inverse stabilization factor of linearized ADMM.
    //The smaller the slower and more stable the ADMM
     // FIXME: why is this not a runtime parameter accessible from the parameter file?
    T gamma_fac=0.1;

    //Read in the image
    read_image(params.imagename, params.gnoise);

    //Create the psf
    create_psf();

    T *cs,*qalpha_ret;
    //A reusable part of the statics generation will be stored here and written to disk
    qalpha_ret = new T[1];

    //At this point we have to decide if we do an approximated dykstra method.
    //The approximated dykstra method will be much faster, but use much less frames
    //for Dykstra's projection algorithm. \n
    //* The exact algorithm uses our classes cset_small and cluster<cset_small> \n
    //* The approximation uses cset_dyadic \n
    //Our extremeValueStatisticsGenerator and CUDADriver classes will use specialized strategies based on the template instantiation

    //Approximated algorithm
    if ( params.do_approx ) {
        typedef cset_dyadic field_patch;
        typedef cset_dyadic driver_patch;

        //Get an instance of the extremeValueStatisticsGenerator class \n
        //* Generates the patches and puts them in clusters if needed \n
        //* Calculates the weights $c_s$
        std::cout << "Creating field\n";
        extremeValueStatisticsGenerator<field_patch, T, gpu_cuda> field(pwidth, pheight, pdepth,  input_image, params.sigma, params.step, 1024);//only kept for statitics generation atm

        std::cout << "Getting quantile\n";
        cs = field.get_quantile_gauss(params.gnoise, alpha_quant, params.step, pwidth, pheight, pdepth, qalpha_ret);

        //Prepare the driver
        std::cout << "Setting up driver\n";
        step35::CUDADriver<driver_patch, T, gpu_cuda> driver(field.croot, input_image, psf, cs, pwidth, pheight, pdepth,
                                                             params.step, gamma_fac, params.sigma, params.regType,
                                                             params.dim);
        //This function templates to whatever Dykstra Flavour we have chosen
        std::cout << "Starting run\n";
        __run<field_patch, driver_patch>(field, driver);
    }
    //Exact algorithm
    else {
#define USE_EXACT
#ifdef USE_EXACT
        typedef cset_small field_patch;
        typedef cluster<cset_small, T> driver_patch;

        //Get an instance of the extremeValueStatisticsGenerator class \n
        //* Generates the patches and puts them in clusters if needed \n
        //* Calculates the weights $c_s$
        std::cout << "Creating field\n";
        extremeValueStatisticsGenerator<field_patch, T, gpu_cuda> field(pwidth, pheight, pdepth, input_image, params.sigma, params.step, 1024);//only kept for statitics generation atm

        std::cout << "Getting quantile\n";
        cs = field.get_quantile_gauss(params.gnoise, alpha_quant, params.step, pwidth, pheight, pdepth, qalpha_ret);

        //Prepare the driver
        std::cout << "Setting up driver\n";
        step35::CUDADriver<driver_patch, T, gpu_cuda> driver(field.cluster_root, input_image, psf, cs, pwidth, pheight, pdepth,
                                                             params.step, gamma_fac, params.sigma, params.regType,
                                                             params.dim);
        //This function templates to whatever Dykstra Flavour we have chosen
        __run<field_patch, driver_patch>(field, driver);
#else
       AssertThrow(false, dealii::ExcNotImplemented());
#endif
    }
}

//@sect5{Function: __run}
//@brief Second part of the run function, templatized
//this enables a unified interface of the driver for exact and approximative method
//but one has to split the run method in two functions
template<typename T> template<typename field_patch, typename driver_patch>
void step35::ADMM<T>::__run (extremeValueStatisticsGenerator<field_patch, T, gpu_cuda> &field,
                             step35::CUDADriver<driver_patch, T, gpu_cuda> &driver) {
    //Used for shifted indexing
    int ti,tj,tk;
    //Used to count up to the next report iteration
    int next_report=0;

    //It is convenient to do the image blurring now
   // FIXME: check dimensions
   std::copy(input_image.begin(), input_image.end(), image_as_read.begin());
    //std::copy(&(image[0]), &(image[pwidth*pheight*pdepth]), orig);

    //Update pwidth and pheight as the driver might add extra padding
    //@param pwidth width as padded by the driver
    //@param pheight height as padded by the driver
    //@param pdepth depth as padded by the driver
    pwidth  = driver.inf->ext_width;
    pheight = driver.inf->ext_height;
    pdepth  = driver.inf->ext_depth;

    driver.get_data();
    //Simulates dataset by adding gaussian noise
    if( params.simulate ) {
        add_gaussian_noise<driver_patch>(driver);
    }

    //Copy to determine relative change in each iteration
    std::vector<T> prev_image(input_image.size(), T(0));
   // T *prev=new T[pwidth*pheight*pdepth];
//    for (int k=0;k<depth;k++)   {
//        for (int i=0; i<pwidth; i++) {
//            for (int j=0; j<pheight; j++) {
//                prev[k*pwidth*pheight+i*pheight+j]=0;
//            }
//        }
//    }
    //Residual of the current step, initialize as residual>tolerance
    T res=TWO*params.tol;

    //Log file
    std::ofstream gain_out("gain.txt");

    //DEBUG
    clock2 = std::chrono::high_resolution_clock::now();
    typedef std::chrono::duration<float> secs;
    std::ofstream dt_out("startup_time");
    dt_out    << std::chrono::duration_cast<secs>(clock2 - clock1).count();
    dt_out.close();

    //Main loop of the ADMM Algorithm
    //Iteration counter
    int iter=1;
    //Will be used to store the constraint violations
    T c1,c2=0;
    while (  ( /*res > params.tol &&*/ iter < params.max_it) || iter < 5 ) {
        //Argmin w.r.t. x
        driver.x_step(params.rho1, params.rho2);
        //driver.x_step_ET(params.rho1, params.rho2);

        //Argmin w.r.t. z
        driver.m_smoothing(params.regInt, params.rho2);

        //Argmin w.r.t. e, is equivalent to a projection
        if ( params.gnoise > 0 )
            driver.dykstra_gauss(params.rho1);

        //Update the lagrange Multipliers
        driver.update_lagrangian(params.alpha1, params.alpha2);

        //Report progress
        if ( iter >= next_report ) {
            //Calculate the change in this step and check the constraints
            res = 0;
            c1 = 0;
            c2 = 0;
            driver.conv2(driver.inf->x_d.array().val(), driver.tmp_d.array().val());
            driver.get_data();
            //Calculate the change w.r.t. the last reported result
            //Calculate the violation of constraints
            for (int i=0; i<pwidth; i++) {
                for (int j=0; j<pheight; j++) {
                    for (int k=0; k<pdepth; k++) {
                        //The image shifts by sigma during the convolution
                        //therefore access it shifted
                        ti = i+params.sigma;
                        tj = j+params.sigma;
                        tk = k+params.sigma;
                        while ( ti >= pwidth )
                            ti=ti-pwidth;
                        while ( tj >= pheight )
                            tj=tj-pheight;
                        while ( tk >= pdepth )
                            tk = tk-pdepth;
                        //Violation of data constraint
                        //if (i < 50)
                          //  if (j < 20)
                        // std::cout << driver.inf->im_h[i*pheight+j] << " " << driver.inf->e_h[i*pheight+j] << " " << driver.tmp_h[ti*pheight+tj] << " " << i << " " << j << std::endl;
                        c1 += boost::math::pow<2>(driver.inf->im_h[k*pwidth*pheight+i*pheight+j]  - driver.inf->e_h[k*pwidth*pheight+i*pheight+j]
                                                  - driver.tmp_h[tk*pwidth*pheight+ti*pheight+tj]);
                        //Violation of smoothness constraint
                        //std::cout << driver.inf->z_h[i*pheight+j] << " " << driver.inf->x_h[i*pheight+j] << " " << i << " " << j << endl;
                        c2 += boost::math::pow<2>(driver.inf->z_h[k*pwidth*pheight+i*pheight+j] - driver.inf->x_h[k*pwidth*pheight+i*pheight+j]);
                        //Residual towards last reported image
                        res+= boost::math::pow<2>(driver.inf->x_h[k*pwidth*pheight+i*pheight+j] - prev_image[k*pwidth*pheight+i*pheight+j]);
                    }
                }
            }
            c1 = c1/((T)(pwidth*pheight*pdepth));
            c2 = c2/((T)(pwidth*pheight*pdepth));
            res=res/((T)(pwidth*pheight*pdepth))+c1;
            // FIXME: vectors in drif.inf
            std::copy(&(driver.inf->x_h[0]), &(driver.inf->x_h[pwidth*pheight*pdepth]), &prev_image[0]);
            gain_out << iter << " " << c1 << " " << c2 << " " <<  std::endl;
            std::cout << "Iteration: " << iter << " | Constraint1: " << c1 << " | Constraint2: "
                      << c2 << " MR-level " << field.split << std::endl;
            next_report = iter + params.report_interval;

            //Print the current estimate to a tiff file
            if (params.do_control)
            {
                QString img_out (QString(params.out_imagename.c_str()).replace(".", QString("-" + QString::number(iter))+"."));

                write_image<T>(img_out.toStdString(), driver.inf->z_h, pheight0, pwidth0, pdepth0, 1.0);//params.gnoise);

            }
        }
        iter=iter+1;
    } // End of main loop

    write_image<T>(params.out_imagename, driver.inf->z_h, pheight0, pwidth0, pdepth0, params.gnoise);
    std::cout << "Done." << std::endl;
}

#endif // STEP35_HH
