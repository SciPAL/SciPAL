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

Copyright  Lutz Künneke, Jan Lebert 2014-2015,
           Stephan Kramer and Johannes Hagemann 2014-2016
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

// deal.II
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

//Boost used for Noise simulation
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/pow.hpp>

//TIFF used for IO
#include<tiff.h>
#include<tiffio.h>

//Our stuff
#include <ADMMParams.h>
#include <cuda_driver_step-35.h>
#include <cuda_driver_step-35.hh>

#include <smre_problem.hh>

// SciPal includes
//
// This is an auxiliary structure which collects data related to the GPUs
// available in the computer on which this program gets executed.
// For details about the cuda<something> functions have a look at
// the CUDA reference manual.
#include <base/GPUInfo.h>

// for shifting of fft
#include <cufftShift/Src/cufftShiftInterface.h>
//DEBUG
#include <chrono>

namespace step35 {

// @sect4{Class: ImageIO}
//
// The input and output of the raw data and the results is managed by a separate class.
// Although this program is primarily used for image processing
// encapsulating the I/O in a class of its own allows to add further types of usages
// while leaving the actual solver untouched.
class ImageIO {

public:
    template<typename T>
    void read_image(std::vector<T> &input_image,
                    std::vector<T> &image_as_read,
                    ImageDoFHandler & dof_handler,
                    const std::string path,
                    const double gnoise,
                    bool use_anscombe=false);

    template<typename PixelDataType>
    void write_image(std::string path,
                     dealii::Vector<PixelDataType> & in,
                     const ImageDoFHandler & dof_handler,
                     double gnoise,
                     bool use_anscombe=false);

// protected:
    //TIFF options, will be saved to make the output identical to the input
    int sampleperpixel, width, height, photom, bps, depth;

    //Input buffer for uint16 tiff image
    uint16* buf_in;
    //pwidth: width after convolution with point spread function
    //pwidth0: width after padding on next multiple of 32
    // FIXME : pwidth0, pheight0 never got padded. This is in ImageInfo, cf. ext_*

    //int pwidth; // , pwidth0;

    //pheight: height after convolution with point spread function
    //pheight0: height after padding on next multiple of 32

    // int pheight; //, pheight0;

    //pdepth: depth after convolution with psf
    //pdepth0: depth after padding on next multiple of 32

    // int pdepth; // , pdepth0;

    int default_strip_size;

    void print() const
    {
        std::cout << "state of ImageIO:\n"
                  << "--------------------"
                  << "samples per pixel : " << sampleperpixel << " \n"
                  << "width : " << width << "\n"
                  // << "pwidth : " << pwidth << "\n"
                 //  << "pwidth0 : " << pwidth0 << "\n"
                  << "height : "<< height <<"\n"
                //  << "pheight : "<< pheight <<"\n"
                 //  << "pheight0 : "<< pheight0 <<"\n"
                  << "default_strip_size : " << default_strip_size <<"\n" << std::endl;
    }

};

//@sect5{Function: read_image}
//
// Since the image defines the computational domain we also pass in the dof_handler and initialize it there.
// @brief Read tiff directory as 2D or 3D image
// @param path: path to the tiff directory
// @param gnoise standard deviation of Gaussian noise
template<typename T>
void step35::ImageIO::read_image(std::vector<T> & input_image,
                                 std::vector<T> & image_as_read,
                                 ImageDoFHandler & dof_handler,
                                 const std::string input_file, const double gnoise,  bool use_anscombe)
{
    std::cout << "Reading in image : " << input_file.c_str() << std::endl;

    // The first we need for reading a tif file is a TIFF object.
    // Since the TIFF library is implemented in plain C we have to declare a
    // raw pointer to it ...
    TIFF *tif_in;

    // and instantiate the object by opening the path to the tif file.
    tif_in=TIFFOpen(input_file.c_str(), "r");

    // Then, we have to get the information about the amount of data to read.
    TIFFGetField(tif_in,TIFFTAG_SAMPLESPERPIXEL,&sampleperpixel);
    TIFFGetField(tif_in,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tif_in,TIFFTAG_IMAGELENGTH,&height);
    TIFFGetField(tif_in,TIFFTAG_PHOTOMETRIC,&photom);
    TIFFGetField(tif_in,TIFFTAG_BITSPERSAMPLE,&bps);

    // By looking up the number of entries in the tif directory we get the depth
    depth=0;
    do
    {
        depth++;
    }
    while ( TIFFReadDirectory(tif_in) );


    // pdepth=depth;
    // pdepth0=depth;
   // pwidth=width;
   // pheight=height;
    // pwidth0=pwidth;
    // pheight0=pheight;

    default_strip_size = TIFFDefaultStripSize(tif_in, /*p*/width /*mwidth*/ *sampleperpixel);

    // FIXME: is this stil true? Is Haar really only for squares?
    // Since our Haar wavelet implementation only works on quadratic images we need a computatonal domain which is quadratic.
    // pwidth=std::max(pwidth, pheight);
    // pheight=std::max(pwidth, pheight);

    // Immediately after reading the dimensions of the image stack we can initialize the @p dof_handler;
    dealii::Point<3> bounding_box( width, height, depth);
    int padding_multiplier = 32;
    dof_handler.reinit(bounding_box, padding_multiplier);




    //reinit the tif
    tif_in=TIFFOpen(input_file.c_str(), "r");
    buf_in=(uint16*)_TIFFmalloc(sizeof(uint16)*width);

    // FIXME: use math-aware Vectors not containers
#ifdef nUSE_DOF_HANDLER
    input_image.resize(pwidth * pheight * pdepth, T(0));
#else
    input_image.resize(dof_handler.n_dofs() /*reinit to padded size*/);
#endif

    image_as_read.resize(input_image.size(), T(0));

    std::cout << "Image sizes right after reading : \n";
    this->print();
    int dd=0;
    do
    {
        for (int y=0; y<height; y++) {
            if(TIFFReadScanline(tif_in,buf_in,y,0) < 0) {
                std::cerr << "reading error\n";
                return;
            }
            for (int x=0; x<width; x++) {
#ifdef nUSE_DOF_HANDLER
                int xyz_pos = dd*pwidth*pheight+x*pheight+y;
#else
                int global_index = dof_handler.global_index(x, y, dd);
                int xyz_pos = global_index;
#endif
                std::ostringstream message("index mismatch; ");
                        message << xyz_pos << " != " << global_index << std::ends;
//                AssertThrow(xyz_pos == global_index, dealii::ExcMessage(
//                                message.str().c_str()
//                                ) );

                // FIXME: Anscombe trafo is not an IO operation. Move somewher else.
                if (use_anscombe)
                {
                    //Do the anscombe transformation
                    input_image[xyz_pos]=sqrt(buf_in[x]+3.0/8.0)*gnoise;
                }
                else {
                    input_image[xyz_pos]=buf_in[x];
                }
#ifdef READ_IMG_DEBUG
                if (x%465==0)
                    printf("I(%d, %d, %d) : %f, buf: %d", x, y, dd, image[xyz_pos], buf_in[x]);
#endif
            }
        }
        std::cout << "Reading tiff layer " << dd << std::endl;
        dd++;
    } while ( TIFFReadDirectory(tif_in) );

      TIFFClose(tif_in);
}


//@sect5{Function: write_image}
//@brief write image to a .tif file
//@param path path to image includign filename
//@param in what to read with dimensions pheight, pwidth
//@param mheight output height
//@param mwidth output width
//@param mdepth number of layers in image stack
//@param gnoise standard deviation of Gaussian noise
template<typename PixelDataType>
void step35::ImageIO::write_image(std::string path,
                                  dealii::Vector<PixelDataType> &in,
                                  const ImageDoFHandler &dof_handler, double gnoise,  bool use_anscombe)
{
    int mwidth = dof_handler.n_dofs_x();

    int mheight = dof_handler.n_dofs_y();

    int mdepth = dof_handler.n_dofs_z();

    uint16* pbuf=(uint16*)_TIFFmalloc(sizeof(uint16)*mwidth);
    TIFF *imout=TIFFOpen(path.c_str(),"w");
    bps=16;

    TIFFSetField(imout, TIFFTAG_ROWSPERSTRIP, default_strip_size);
    TIFFSetField(imout, TIFFTAG_IMAGEWIDTH, mwidth);  // set the width of the image
    TIFFSetField(imout, TIFFTAG_IMAGELENGTH, mheight);    // set the height of the image
    TIFFSetField(imout, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
    TIFFSetField(imout, TIFFTAG_BITSPERSAMPLE, bps);    // set the size of the channels
    TIFFSetField(imout, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image->
    TIFFSetField(imout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, photom);
    TIFFSetField(imout, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(imout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    for (int z=0; z<mdepth; z++)
    {
        for (int y=0; y<mheight; y++)
        {
            for (int x=0; x<mwidth; x++)
            {
                #ifdef nUSE_DOF_HANDLER
                int xyz_pos = z*pheight*pwidth+x*pheight+y;
#else
                int xyz_pos = dof_handler.global_index(x, y, z);
#endif
                PixelDataType intensity_value_xyz = in[xyz_pos]; //FIXME think about negative values
                if ( intensity_value_xyz > 0 )
                {    // FIXME: reenable Anscombe trafo
                    if (use_anscombe)
                    {
                        // Inverse anscombe transformation
                        pbuf[x]= (uint16) round((intensity_value_xyz
                                               *intensity_value_xyz)/gnoise);
                    }
                    else
                        pbuf[x]= (uint16)round(intensity_value_xyz);
                }
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
            TIFFSetField(imout, TIFFTAG_ROWSPERSTRIP, default_strip_size); // TIFFDefaultStripSize(tif, mwidth*sampleperpixel));
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


//@sect4{Class: ADMM}
//
// To make this test facility extendible, we implement
// a class for a simple user interface. Its primary tasks are
// - management of run-time parameters by a simple text-based parameter file
// - setting device parameters according to the user's parameters
// - preprocessing and output of results
template<typename T, typename BW>
class ADMM {
public:

    ADMM(int argc, char *argv[], SciPAL::GPUInfo &g);
    void run();

    T _psf(T x, T y, T z, T sigma);
    void create_psf();

    void add_blur_and_gaussian_noise(step35::CUDADriver<T, BW> &driver);

    // The ADMM algoritmh is put into a separate class which inherits from deal.II's base class
    // for iterative solvers. The advantage of this is, that we do not have to explain the design in detail
    // as this is already done in the deal.II doc.
    class ADMMSolver : public dealii::Solver<SciPAL::Vector<T, BW> > {

        typedef SciPAL::Vector<T, BW> VECTOR;

        typedef dealii::Solver< > Base;

    public:
        ADMMSolver (dealii::SolverControl & solver_control) : Base(solver_control) {}

        template<class MATRIX , class PRECONDITIONER >
        void 	solve (const MATRIX &A, VECTOR &x, const VECTOR &b, const PRECONDITIONER &precondition)
        {

        }

    protected:
        std::vector<T> prev_image; // (input_image.size(), T(0));

    };


private:
    SciPAL::GPUInfo & gpuinfo;

    ImageDoFHandler dof_handler;

    //Array holding the signal
    // FIXME: store images in QImage. Applies to other attributes as well.
    std::vector<T> input_image;
    //Array holding the signal as read from the tif
    std::vector<T> image_as_read;
    //Array holding the point spread function
    dealii::Vector<T> psf;

    //DEBUG
    std::chrono::high_resolution_clock::time_point clock1;
    std::chrono::high_resolution_clock::time_point clock2;

    void __run(step35::CUDADriver<T, BW> &driver);

    ImageIO image_io;
protected:
    SimParams params;


    template<typename Driver>
    void set_initial_condition(Driver &driver);
};
}

// @sect5{Constructor: ADMM}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
template<typename T, typename BW>
step35::ADMM<T, BW>::ADMM(int argc, char *argv[], SciPAL::GPUInfo &g)
    : gpuinfo(g)
{
    //DEBUG
    clock1 = std::chrono::high_resolution_clock::now();

    // Declare and read parameters from a file. Basically, the parameter
    // file must have the same name as the binary. The extension has
    // to be ".prm". What has been read will be dumped into a log file.
    dealii::ParameterHandler prm_handler;

    SimParams::declare(prm_handler);

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

    // Create toplevel run directory

    this->params.run_dir.setPath( this->params.run_dir.absolutePath() +  QDir::separator() +
    QString(QDateTime::currentDateTime().toString("yyyy-MM-dd/hh_mm_ss")
            ).remove("."));

    cwd.setPath(this->params.run_dir.absolutePath());

    std::cout << "path to run directory : " << this->params.run_dir.absolutePath().toStdString().c_str() << std::endl;


    // The following lets a directory make its own path.
    if (!cwd.exists())
        cwd.mkpath( "." );


    // After the run directory we create the log directory.
    this->params.prm_log_dir = this->params.run_dir.absolutePath() + QDir::separator() + "log";
    if (!this->params.prm_log_dir.exists())
        this->params.prm_log_dir.mkpath(".");

    std::cout << "log path : " << this->params.prm_log_dir.absolutePath().toStdString().c_str() << std::endl;

    // Now, change to the run directory
    QDir::setCurrent(cwd.absolutePath());

    // ... and write what has been actually read
    // into log file. Basically, this is just another parameter file
    // and thus could be used again as input to another run after stripping the .log suffix.
    std::string master_log_file = (this->params.prm_log_dir.absolutePath() + QDir::separator()
                                   +
                                  (QFileInfo(prm_filename.c_str()).fileName()
                                   + ".log") ).toStdString();

    std::cout << "log file : " << master_log_file.c_str()
                                   << std::endl;

    std::ofstream log_master_prm( master_log_file.c_str() );
    prm_handler.print_parameters (log_master_prm,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    QDir::setCurrent(this->params.run_dir.absolutePath());

    //set number of threads
    omp_set_num_threads(params.n_omp_threads);
}


// @sect5{Function: set_initial_condition}
//
//
template<typename T, typename BW>
template<typename Driver>
void step35::ADMM<T, BW>::set_initial_condition(Driver &driver)
{
    int depth = dof_handler.n_dofs_z();

    int width = dof_handler.n_dofs_x();

    int height = dof_handler.n_dofs_y();


    int __ext_width = dof_handler.n_dofs_x_padded();
    int __ext_height = dof_handler.n_dofs_y_padded();


    // FIXME: I do not want this written in components
    for (int z=0; z<depth; z++) {
        for (int x=0; x<width; x++) {
            for (int y=0; y<height; y++) {
                if ( z < depth && x < width && y < height )
                    driver.im_h[z*__ext_width*__ext_height + x*__ext_height + y] =
                            input_image[z*width*height + x*height + y];
            }
        }
    }

    driver.x_h = driver.im_h;

//    if ( arch == gpu_cuda ) {

        driver.writeable_im_d() = driver.im_h;
        driver.x_d = driver.x_h;
        driver.z_d = driver.z_h;
        driver.writeable_e_d() = driver.e_h;

//        cs_d = cs_vec;
//        //Copy $c_s$ to constant CUDA memory
//        step35::Kernels<Mdouble> kernel;
//        kernel.set_cs(const_cast<Mdouble*>(&cs_vec[0]), cs_n);

//        std::cout << "weights c_s : \n";
//        cs_d.print();

     //   checkCudaErrors(cudaDeviceSynchronize());
//    }

}
//@sect5{Function: add_gaussian_noise}
//@brief Simulates dataset by adding gaussian noise, the whole driver is given to used on device convolution
template<typename T, typename BW>
void step35::ADMM<T, BW>::add_blur_and_gaussian_noise (step35::CUDADriver<T, BW> &driver) {
    //Constant seed to get reproducable tests
    boost::mt19937 rng;
    //Seed rng
    if ( params.time_seed )
        rng.seed(time(NULL));
    //Prepare a normal distribution with mean 0 and standard deviation params.gnoise
    boost::normal_distribution<> nd(0.0, params.gnoise);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);



   driver.convolution.vmult(driver.writeable_im_d(), driver.writeable_im_d() ); //

    //Push data from device to host
    driver.get_data();

    driver.generated_noise.reinit(driver.im_h);

    dealii::Vector<T> gn_offset(driver.generated_noise.size());
    gn_offset = 0.;

    driver.generated_noise = 0.;
    //Now add gaussian noise
#ifdef nUSE_DOF_HANDLER
    for (int z=0; z<image_io.pdepth; z++) {
        for (int x=0; x<image_io.pwidth; x++) {
            for (int y=0; y<image_io.pheight; y++) {
                int xyz_pos = z*image_io.pwidth*image_io.pheight + x*image_io.pheight + y;
#else
    for (int z=0; z<dof_handler.n_dofs_z(); z++ ) {
        for (int x=0; x<dof_handler.n_dofs_x_padded(); x++) {
            for (int y=0; y<dof_handler.n_dofs_y_padded(); y++)
            {
                int xyz_pos = dof_handler.global_index(x, y, z);
#endif


              driver.generated_noise[xyz_pos] = var_nor();
               driver.im_h[xyz_pos] += driver.generated_noise[xyz_pos];
               //shift output of noise distribution to avoid negative numbers
               gn_offset[xyz_pos] = 6*params.gnoise + driver.generated_noise[xyz_pos];
               // Add offset so that we something when plotting the noise component alone.

                //The simulated data is meant to represent a noisy image, which is naturally non-negative everywhere
                if (  driver.im_h[xyz_pos] < 0)
                   driver.im_h[xyz_pos] = 0;
            }
        }
    }

    //Push data from host to device
    driver.push_data();

    //Write the noisy image for control reasons
    image_io.write_image("input.tif", driver.im_h,
                         // image_io.pheight, image_io.pwidth, image_io.pdepth,
                         dof_handler,
                         params.gnoise, params.anscombe);

    image_io.write_image("noise_component.tif", gn_offset,
                         // image_io.pheight, image_io.pwidth, image_io.pdepth,
                         dof_handler,
                         params.gnoise, false/*params.anscombe*/);

    image_io.write_image("psf.tif", driver.convolution.psf_h,
                         // image_io.pheight, image_io.pwidth, image_io.pdepth,
                         dof_handler,
                         params.gnoise, params.anscombe);
//    dealii::Vector<float> bla(input_image.begin(), input_image.end());// bla= (input_image);
//    image_io.write_image("bla.tif", bla,
//                         // image_io.pheight, image_io.pwidth, image_io.pdepth,
//                         dof_handler,
//                         params.gnoise, params.anscombe);




    std::cout << "noise written" << std::endl;

    // Set the initial value of x to the simulated input.
    // driver.x_d = driver.im_d();
}

// @sect5{Function: run}
//
// Read in settings and data and perform the algorithm
template<typename T, typename BW>
void step35::ADMM<T, BW>::run() {

    std::cout << "Starting run" << std::endl;

    //Read in the image
    image_io.read_image(this->input_image, this->image_as_read, dof_handler,
                        params.input_image, params.gnoise, params.anscombe);

    // Debugging info about image sizes.
    image_io.print();


    std::vector<T> cs;


    //At this point we have to decide if we do an approximated dykstra method.
    //The approximated dykstra method will be much faster, but use much less frames
    //for Dykstra's projection algorithm. \n
    //* The exact algorithm uses our classes cset_small and cluster<cset_small> \n
    //* The approximation uses cset_dyadic \n
    //Our extremeValueStatisticsGenerator and CUDADriver classes will use specialized strategies based on the template instantiation

    //Approximated algorithm
    if ( params.do_approx ) {
        //Prepare the driver
        std::cout << "Setting up driver\n";
        step35::CUDADriver<T, BW> driver(cs,
                                                     #ifdef nUSE_DOF_HANDLER
                                                             image_io.pwidth, image_io.pheight, image_io.pdepth,
                                                     #else
                                                             dof_handler,
                                                     #endif
                                                             params);

        this->set_initial_condition(driver);


        //This function templates to whatever Dykstra Flavour we have chosen
        std::cout << "Starting run\n";
        __run(driver);
    }
    //Exact algorithm
    else {

        AssertThrow(false, dealii::ExcNotImplemented());

    }
}

//@sect5{Function: __run}
//@brief Second part of the run function, templatized
//this enables a unified interface of the driver for exact and approximative method
//but one has to split the run method in two functions
template<typename T, typename BW>
void step35::ADMM<T, BW>::__run (step35::CUDADriver<T, BW> &driver)
{
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
#ifdef nUSE_DOF_HANDLER
    image_io.pwidth  = driver.ext_width();
    image_io.pheight = driver.ext_height();
    image_io.pdepth  = driver.ext_depth();
#endif


    driver.get_data();
    //Simulates dataset by adding gaussian noise
    if( params.simulate ) {
        add_blur_and_gaussian_noise(driver);
    }

    //Copy to determine relative change in each iteration
    dealii::Vector<T> prev_image(input_image.size());
    prev_image = 0;

    SciPAL::Vector<T, BW> res_tmp1, res_tmp2;

    //Residual of the current step, initialize as residual>tolerance
    T res=2.0*params.solver_control.tolerance();

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
    int iter = 0;
    //Will be used to store the constraint violations
    T c1 = 0, c2 = 0;

    {
        std::vector<T> zero_init(driver.x_d.size(), 0);
        driver.x_d = zero_init; //FIXME: try init with driver.im_h
        driver.z_d = zero_init;
    }

    // set initial conditions
    driver.lag2 = driver.im_h;

    // FIXME: unsatisfying solution to disable the checking of the residual in the sc
    params.solver_control.clear_failure_criterion();
    while (  // ( /*res > params.tol &&*/ iter < params.max_it) || iter < 5
             iter < params.solver_control.max_steps()
    //       params.solver_control.check(iter, res) == dealii::SolverControl::iterate
             )
    {
        // std::cout << "entering step " << iter << std::endl;


        //Argmin w.r.t. x
   driver.x_step_adaptive(params);


        //Argmin w.r.t. e, is equivalent to a projection   
        if ( params.gnoise > 0 )
            driver.dykstra_gauss(params.rho1);

        //Update the Lagrange Multipliers


        driver.update_lagrangian(
                    1./
                                 params.inv_gamma
                                 ,
                                 params.alpha2);


        //Report progress
        // iter = params.solver_control.last_step();
        if ( iter <3 || (iter+1) % params.report_interval == 0 )
        {
            //Calculate the change in this step and check the constraints
            res = 0;
            c1 = 0;
            c2 = 0;
            driver.convolution.vmult(driver.tmp_d, driver.x_d); // conv2(driver.x_d.array().val(), driver.tmp_d.array().val());
            driver.get_data();
            //Calculate the change w.r.t. the last reported result
            //Calculate the violation of constraints

 // driver.im_d - driver.e_d - driver.tmp_d
            // c1  = L2_norm (driver.im_h - driver.e_h - driver.tmp_h);

            driver.residual(res_tmp1, driver.im_d(), driver.e_d(), driver.x_d);

            c1 = res_tmp1.l2_norm();

            // c2  = L2_norm (driver.z_h - driver.x_h);
            res_tmp2 = driver.z_d;
            res_tmp2 -= driver.x_d;
            c2 = res_tmp2.l2_norm();

            // res = L2_norm (driver.x_h - prev_image);
            res_tmp1 = driver.x_d;
            res_tmp2 = prev_image;
            res_tmp1 -= res_tmp2;
            res = res_tmp1.l2_norm();

            c1 = c1/std::sqrt( (T)(dof_handler.n_dofs()) );
            c2 = c2/std::sqrt( (T)(dof_handler.n_dofs()) );
            res=res/std::sqrt( (T)(dof_handler.n_dofs()) ) + c1;

            res_tmp1 = driver.generated_noise;
            res_tmp1 -= driver.e_d();
            T c3 = res_tmp1.l2_norm()/std::sqrt( (T)(dof_handler.n_dofs()) );

            // FIXME: vectors in drif.inf

            // std::copy(driver.x_h.begin(), driver.x_h.end(), prev_image.begin() );
            driver.x_d.push_to(prev_image);
            gain_out << iter << " " << c1 << " " << c2 << " " <<  std::endl;
            std::cout << "Iteration: " << iter << " | Constraint1 ||im - e - tmp||_2   : " << c1 << " | Constraint2 : ||x - z||_2  : "
                      << c2 << ", L2-distance to generated noise : " << c3 << std::endl;


            //Print the current estimate to a tiff file
            if (params.do_control)
            {
                // As name for the intermediate result we use what was given in the parameter file and append the zero-padded iteration counter.
                QString img_out (QString(params.out_imagename.c_str()).replace(".", QString("-" +
                                                                                            QString("%1").arg(iter, 6, 10, QChar('0'))
                                                                                                      )+"."));

                // FMM version gives the wrong sign, hence : driver.z_h *= -1;
                image_io.write_image(img_out.toStdString(), driver.z_h,
                                     dof_handler, // image_io.pheight0, image_io.pwidth0, image_io.pdepth0,
                                     1.0,//params.gnoise
                                     params.anscombe);



#ifdef X_STEP_DEBUG
        // print contents of x_d and x_h before step


        {
            QString x_d_fname = QString ("im_d-%1d.tif").arg(iter);

            dealii::Vector<T> x_d_tmp(driver.im_d().size());
            driver.im_d().push_to(x_d_tmp);
            image_io.write_image(x_d_fname.toStdString(), x_d_tmp, dof_handler, 1.0 /*dummy value for Gaussian noise*/, params.anscombe);
        }
#else

        {
            if (false)
            {
                QString x_d_fname = QString ("x_d-%1d.tif").arg(iter);

                dealii::Vector<T> x_d_tmp(driver.x_d.size());
                driver.x_d.push_to(x_d_tmp);
                image_io.write_image(x_d_fname.toStdString(), x_d_tmp, dof_handler, 1.0 /*dummy value for Gaussian noise*/, params.anscombe);
            }

            if (false)
            {
                QString x_d_fname = QString ("z_d-%1d.tif").arg(iter);

                dealii::Vector<T> x_d_tmp(driver.z_d.size());
                driver.z_d.push_to(x_d_tmp);
                image_io.write_image(x_d_fname.toStdString(), x_d_tmp, dof_handler, 1.0 /*dummy value for Gaussian noise*/, params.anscombe);
            }

              if (true)
            {
                QString x_d_fname = QString ("e_d-%1.tif").arg(iter, 6, 10, QChar('0'));

                dealii::Vector<T> x_d_tmp(driver.e_d().size());
                driver.e_d().push_to(x_d_tmp);
                x_d_tmp.add(6*params.gnoise);
                image_io.write_image(x_d_fname.toStdString(), x_d_tmp
                                     , dof_handler, 1.0 /*dummy value for Gaussian noise*/, params.anscombe);
            }
        }
#endif



        }
        }
        iter++;
    } // End of main loop

       image_io.write_image("final-" + params.out_imagename, driver.z_h,
                            dof_handler,
                            // image_io.pheight0, image_io.pwidth0, image_io.pdepth0,
                            params.gnoise, params.anscombe);
    std::cout << "Done." << std::endl;
}

#endif // STEP35_HH
