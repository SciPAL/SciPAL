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

Copyright  S. C. Kramer , J. Hagemann  2010 - 2016
*/


#ifndef SIMPARAMS_H
#define SIMPARAMS_H

#include <QDir>
#include <QTime>


// deal.II
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/base/point.h>

// SciPAL
#include <base/ParallelArch.h>


//TIFF used for IO
#include<tiff.h>
#include<tiffio.h>


namespace step35 {






// @sect3{Class: ImageDoFHandler}
//
// Images or stacks of them (aka movies) are bounded domains in 2D or 3D. As there numerical representation is
// 1-dimensional, we need a structure which encapsulates the mapping from the multidimensional coordinates to the
// linear memory address. In analogy to finite element methods this is called a degrees-of-freedom handler.
// template<int dim>
class ImageDoFHandler {

public:
    static const int dim = 3;
    typedef dealii::Point<dim> Point;

    // The default CTor initializes everything to zero and requires to call reinit() at some later time.
    ImageDoFHandler() {}

    // The CTor takes the true dimensions of the image stack
    // and additionally a flag, whether these should be extended for accomodating some padding
    // in order to meet particular rquirements concerning memory alignment.
    // If the padding multiplier is set to a positive value, then this value is used for
    // computing the leading dimension for the fastest running index, that is x.
    // On input the z component of the bounding box has to be at least 1 as a 2D image represents one plane in 3D.
    ImageDoFHandler(const Point& bounding_box, int padding_mult=-1)
        :
          bbox(bounding_box),
          padding_multiplier(padding_mult)
    {
        this->reinit(bounding_box, padding_multiplier);
    }

    // The copy constructor for copying the DoF information around easily.
    ImageDoFHandler(const ImageDoFHandler & other)
        :
          bbox(other.bbox),
          padding_multiplier(other.padding_multiplier)
    {
        this->reinit(other.bbox, other.padding_multiplier);
    }

    // The reinit function does the same as the constructor.
    // On input the z component of the bounding box has to be at least 1 as a 2D image represents one plane in 3D.
    void reinit(const Point& bounding_box, int padding_mult=-1)
    {
        this->bbox = bounding_box;
        this->padding_multiplier = padding_mult;

        _leading_dim_x = int(bbox[0]);
        _leading_dim_y = int(bbox[1]);

        if (padding_multiplier > 0)
        {
            _leading_dim_x = padding_multiplier * (int(bbox[0] + padding_multiplier - 1)/padding_multiplier);
            _leading_dim_y = padding_multiplier * (int(bbox[1] + padding_multiplier - 1)/padding_multiplier);
        }

        // Last but not least, check whether all values are positive.
        AssertThrow(bbox[0] > 0, dealii::ExcMessage("x-dim is zero or negative."));
        AssertThrow(bbox[1] > 0, dealii::ExcMessage("y-dim is zero or negative."));
        AssertThrow(_leading_dim_y > 0, dealii::ExcMessage("y-dim is zero or negative."));
        AssertThrow(bbox[2] > 0, dealii::ExcMessage("z-dim is zero or negative."));

        _n_dofs_per_plane = _leading_dim_x * _leading_dim_y;
        _n_dofs = _n_dofs_per_plane * n_dofs_z();
        _n_active_dofs = n_dofs_x() * n_dofs_y() * n_dofs_z();
    }


    // We still need a few helper functions for convencience.
    int n_dofs() const { return _n_dofs; }

    int n_active_dofs() const { return _n_active_dofs; }

    // The true width of an image in the 2D plane.
    int n_dofs_x() const { return bbox[0]; }

    // The width after padding.
    int width_p() const { return n_dofs_x_padded(); }

    // The true height of an image in the 2D plane and its padded counterpart.
    int n_dofs_y() const { return bbox[1]; }
    int height_p() const { return n_dofs_y_padded(); }

    // The height of a stack of images. Also designated as depth.
    int n_dofs_z() const { return bbox[2]; }
    int depth_p() const { return n_dofs_z(); }

    // Repeats of width_p() and height_p() for situations when talking about the number of DoFs is more intuitive.
    int n_dofs_x_padded() const { return _leading_dim_x; }

    int n_dofs_y_padded() const { return _leading_dim_y; }



    // For the lexicographic ordering we assume @p y to the fastest running index, followed by @p x and then @p z.
     // Within the 2D plane ordering is thus column-major.
    int global_index(const int x, const int y, const int z=0) const
    {
        // Exceptions are expensive, hence we test the correctness of the input only in DEBUG mode.
        Assert(x >= 0 && x < _leading_dim_x, dealii::ExcIndexRange(x, 0, _leading_dim_x));
        Assert(y >= 0 && y < _leading_dim_y, dealii::ExcIndexRange(y, 0, _leading_dim_y));
        Assert(z >= 0 && z < bbox[2], dealii::ExcIndexRange(z, 0, bbox[2]));

        return z * _n_dofs_per_plane
                +
                x * _leading_dim_y
                +
                y;
    }


    void print() const
    {
        std::cout << " w " << this->n_dofs_x()
                  << " h " << this->n_dofs_y()
                  << " d " << this->n_dofs_z()
                  << " leading dim x " << this->n_dofs_x_padded()
                  << " leading dim y " << this->n_dofs_y_padded()
                  << std::endl;

    }

protected:
    int _leading_dim_x, _leading_dim_y, _n_dofs, _n_dofs_per_plane, _n_active_dofs;
    Point bbox;
    int padding_multiplier;
};



// This structure unites data and dimensions. Similar to QImage, except that it is capable of holding 3D images.
struct TiffImage {

    TiffImage () {}

    // The ordinary copy CTor just copies the other image @p o.
    // If @p copy_data is set to false, only the dimensions are copied but not the
    // content of the pixels.
    TiffImage (const TiffImage & o, bool copy_data = true)
        :
          dofs(o.dofs),
          sample_per_pixel(o.sample_per_pixel),
          photo_m(o.photo_m),
          bps(o.bps),
          default_strip_size(o.default_strip_size)
    {
        if (copy_data)
            data = o.data;
    }


//    TiffImage(const std::string & path)
//    {
//        ImageIO().read_image(*this, path)
//    }

    ImageDoFHandler dofs;

    std::vector<float> data;



    // ==================

    // FIXME: the following attributes must be part of an image and not of the IODevice
    // TIFF options, will be saved to make the output identical to the input
    int sample_per_pixel, photo_m, bps;

    int default_strip_size;


    void print() const
    {
        std::cout << "state of TIFFImage:\n"
                  << "--------------------"
                  << "samples per pixel : " << sample_per_pixel << " \n"
                  << "bps : " << bps << "\n"
                  << "photo_m : "<< photo_m <<"\n"
                  << "default_strip_size : " << default_strip_size <<"\n" << std::endl;
    }


    // ==================

};





// @sect4{Class: ImageIO}
//
// The input and output of the raw data and the results is managed by a separate class.
// Although this program is primarily used for image processing
// encapsulating the I/O in a class of its own allows to add further types of usages
// while leaving the actual solver untouched.
class ImageIO {

public:
    // template<typename T>
    static void read_image(TiffImage &image_to_read,
                           const std::string path);

    // template<typename PixelDataType>
    static void write_image(const std::string path,
                            const TiffImage & image_to_write);

// protected:



};

//@sect5{Function: read_image}
//
// Since the image defines the computational domain we also pass in the dof_handler and initialize it there.
// @brief Read tiff directory as 2D or 3D image
// @param path: path to the tiff directory
// @param gnoise standard deviation of Gaussian noise
// template<typename T>
void step35::ImageIO::read_image(TiffImage & image_to_read,
                                 const std::string input_file)
{
    std::cout << "Reading in image : " << input_file.c_str() << std::endl;

    // Abbreviation
    ImageDoFHandler & dof_handler = image_to_read.dofs;


    // The first we need for reading a tif file is a TIFF object.
    // Since the TIFF library is implemented in plain C we have to declare a
    // raw pointer to it ...
    TIFF *tif_in;

    // and instantiate the object by opening the path to the tif file.
    tif_in = TIFFOpen(input_file.c_str(), "r");

    int width, height;

    // Then, we have to get the information about the amount of data to read.
    TIFFGetField(tif_in, TIFFTAG_SAMPLESPERPIXEL, &image_to_read.sample_per_pixel);
    TIFFGetField(tif_in, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif_in, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif_in, TIFFTAG_PHOTOMETRIC, &image_to_read.photo_m);
    TIFFGetField(tif_in, TIFFTAG_BITSPERSAMPLE, &image_to_read.bps);

    // By looking up the number of entries in the tif directory we get the depth
    int depth = 0;
    do
    {
        depth++;
    }
    while ( TIFFReadDirectory(tif_in) );

    // actual length of scan line
    image_to_read.default_strip_size = TIFFDefaultStripSize(tif_in, width *image_to_read.sample_per_pixel);

    // Immediately after reading the dimensions of the image stack we can initialize the @p dof_handler;
    dealii::Point<3> bounding_box(width, height, depth);

    // Adjust to CUDA mem constraints
    int padding_multiplier = 32;
    image_to_read.dofs.reinit(bounding_box, padding_multiplier);

    //reinit the tif
    tif_in=TIFFOpen(input_file.c_str(), "r");

    //Input buffer for uint16 tiff image
    uint16* buf_in;

    buf_in = (uint16*)_TIFFmalloc(sizeof(uint16)*width);

    // FIXME: use math-aware Vectors not containers
    image_to_read.data.resize(image_to_read.dofs.n_dofs() /*reinit to padded size*/, 0);


    std::cout << "Image sizes right after reading : \n";
    image_to_read.print();
    int dd=0;
    do
    {
        for (int y=0; y < height; y++)
        {
            if(TIFFReadScanline(tif_in, buf_in, y, 0) < 0)
            {
                std::cerr << "reading error\n";
                return;
            }
            // FIXME: use std::copy
            for (int x=0; x < width; x++) {

                int global_index = dof_handler.global_index(x, y, dd);
                int xyz_pos = global_index;

                std::ostringstream message("index mismatch; ");
                        message << xyz_pos << " != " << global_index << std::ends;
                Assert(xyz_pos == global_index, dealii::ExcMessage(
                                message.str().c_str()
                                ) );

                // FIXME: Anscombe trafo is not an IO operation. Move somewher else.
//                if (use_anscombe)
//                {
//                    //Do the anscombe transformation
//                    input_image[xyz_pos]=sqrt(buf_in[x]+3.0/8.0)*gnoise;
//                }
//                else {
//                    input_image[xyz_pos]=buf_in[x];
//                }

                         image_to_read.data[xyz_pos] = buf_in[x];
#ifdef READ_IMG_DEBUG
                if (x%465==0)
                    printf("I(%d, %d, %d) : %f, buf: %d", x, y, dd, image[xyz_pos], buf_in[x]);
#endif
            }

            // Padd by using periodic boundary conditions.
             for (int x = width; x < dof_handler.width_p(); x++)
             {
                 int global_index = dof_handler.global_index(x, y, dd);
                 int xyz_pos = global_index;

                 std::ostringstream message("index mismatch; ");
                         message << xyz_pos << " != " << global_index << std::ends;
                 Assert(xyz_pos == global_index, dealii::ExcMessage(
                                 message.str().c_str()
                                 ) );
                   image_to_read.data[xyz_pos] = buf_in[x - width];
             }
        }

//        // Padd bottom rows
        for (int y = height; y < dof_handler.height_p(); y++)
        {
            for (int x = 0; x < dof_handler.width_p(); x++)
            {
           int src_xyz = dof_handler.global_index(x, y- height, dd);
           int dst_xyz = dof_handler.global_index(x, y, dd);

             image_to_read.data[dst_xyz] = image_to_read.data[src_xyz];

           }
        }

        std::cout << "Reading tiff layer " << dd << std::endl;
        dd++;
    } while ( TIFFReadDirectory(tif_in) );

    _TIFFfree(buf_in);

    TIFFClose(tif_in);
}


//@sect5{Function: write_image}
//@brief write image to a .tif file
//@param path path to image includign filename
//@param in what to read with dimensions height_p, width_p
//@param mheight output height
//@param mwidth output width
//@param mdepth number of layers in image stack
//@param gnoise standard deviation of Gaussian noise
void step35::ImageIO::write_image(const std::string path,
                                  const TiffImage & image_to_write
                                  // dealii::Vector<PixelDataType> &out_data,
                                  // const ImageDoFHandler &dof_handler
                                  )
{

    typedef float PixelDataType;

    const ImageDoFHandler & dof_handler = image_to_write.dofs;

    const std::vector<float> & out_data = image_to_write.data;

    AssertThrow(out_data.size() > 0, dealii::ExcMessage("Bummer"));

    int out_width = dof_handler.n_dofs_x_padded(); // dof_handler.n_dofs_x();

    int out_height = dof_handler.n_dofs_y_padded(); // dof_handler.n_dofs_y();

    int out_depth = dof_handler.n_dofs_z();

    uint16* pbuf=(uint16*)_TIFFmalloc(sizeof(uint16)*out_width);
    TIFF *tiff_out=TIFFOpen(path.c_str(),"w");

    image_to_write.print();

    // FIXME
    TIFFSetField(tiff_out, TIFFTAG_ROWSPERSTRIP, image_to_write.default_strip_size);
    TIFFSetField(tiff_out, TIFFTAG_IMAGEWIDTH, out_width);  // set the width of the image
    TIFFSetField(tiff_out, TIFFTAG_IMAGELENGTH, out_height);    // set the height of the image
    TIFFSetField(tiff_out, TIFFTAG_SAMPLESPERPIXEL, image_to_write.sample_per_pixel);   // set number of channels per pixel
    TIFFSetField(tiff_out, TIFFTAG_BITSPERSAMPLE, image_to_write.bps);    // set the size of the channels
    TIFFSetField(tiff_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image->
    TIFFSetField(tiff_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff_out, TIFFTAG_PHOTOMETRIC, image_to_write.photo_m);
    TIFFSetField(tiff_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tiff_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    for (int z=0; z<out_depth; z++)
    {
        for (int y=0; y<out_height; y++)
        {
            for (int x=0; x<out_width; x++)
            {
                #ifdef nUSE_DOF_HANDLER
                int xyz_pos = z*height_p*width_p+x*height_p+y;
#else
                int xyz_pos = dof_handler.global_index(x, y, z);
#endif
                PixelDataType intensity_value_xyz = out_data[xyz_pos]; //FIXME think about negative values
                if ( intensity_value_xyz > 0 )
                {    // FIXME: reenable Anscombe trafo
//                    if (use_anscombe)
//                    {
//                        // Inverse anscombe transformation
//                        pbuf[x]= (uint16) round((intensity_value_xyz
//                                               *intensity_value_xyz)/gnoise);
//                    }
//                    else
                        pbuf[x]= (uint16)round(intensity_value_xyz);
                }
                else
                    pbuf[x]=0;
            }
            if (TIFFWriteScanline(tiff_out, pbuf, y, 1) < 0) {
                std::cerr << "writing error" << std::endl;
                break;
            }
        }
        if ( z < out_depth-1 ) {
            // Prepare next layer
            TIFFWriteDirectory(tiff_out);
            TIFFSetField(tiff_out, TIFFTAG_ROWSPERSTRIP, image_to_write.default_strip_size); // TIFFDefaultStripSize(tif, mwidth*sampleperpixel));
            TIFFSetField(tiff_out, TIFFTAG_IMAGEWIDTH, out_width);  // set the width of the image
            TIFFSetField(tiff_out, TIFFTAG_IMAGELENGTH, out_height);    // set the height of the image
            TIFFSetField(tiff_out, TIFFTAG_SAMPLESPERPIXEL, image_to_write.sample_per_pixel);   // set number of channels per pixel
            TIFFSetField(tiff_out, TIFFTAG_BITSPERSAMPLE, image_to_write.bps);    // set the size of the channels
            TIFFSetField(tiff_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image->
            TIFFSetField(tiff_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tiff_out, TIFFTAG_PHOTOMETRIC, image_to_write.photo_m);
            TIFFSetField(tiff_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(tiff_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        }
    }
    TIFFClose(tiff_out);
    _TIFFfree(pbuf);
}






// ========================================



// @sect4{enum: Regularization}
// @param haar regularization by haar wavelet sparsity
// @param sparse regularization by direct space sparsity
// @param quadratic Minimize 2-Norm
enum Regularization { sparse/*L1*/, quadratic /* L2*/, TV, Haar };



// We use deal.II's Subscriptor as base class in order to use deal.II's SmartPointer mechanism.
struct AugmentedLagrangianParams : public dealii::Subscriptor {


// Data and parameters for actual inverse problem
    std::string path_to_measured_image;

   TiffImage image_as_read;

    //Wether or not to simulate noise on a test image
    bool simulate; // TODO: add to declare/get
    //Dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    int dim;  // TODO: add to declare/get

    //Do an anscombe before and after
    bool anscombe; // TODO: add to declare/get


    double psf_fwhm;
    double std_dev_noise;

    // Wether or not to use time seeded random numbers
    // when generating synthetic test data
    bool time_seed; // TODO: add to declare/get


    // Parameters for abstract cost functional
    Regularization reg_type;

    double regularization_strength;

    double penalty;



    // First of all, we declare the most important parameter of the whole algorithm: How
    // probable it is that the residual is follows a chi-squared distribution.
    double alpha_quantile;

    // Threshold for L2-norm of solution increment in Dykstra.
    double dykstra_Tol;

    // Maximum number of Dykstra iterations
    int n_max_dykstra_steps;



    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);
};

void AugmentedLagrangianParams::declare(dealii::ParameterHandler &prm)
{
    prm.declare_entry ("FWHM of PSF",
                       "3",
                       dealii::Patterns::Double(),
                       "Full-width half-maximum of the PSF in units of pixeledge lnegths. In case of a PSF to narrow to be reasonably sampled "
                       "on the lattice of pixels we simply assume that it is a delta peak "
                       "and consider the convolution operator as identity.");


    prm.declare_entry ("Std. dev. of Gaussian noise",
                       "0",
                       dealii::Patterns::Double(),
                       "Estimate of the standard deviation of Gaussian noise. "
                       "If simulate = true, it will be used to add Gaussian noise with the given value as "
                       "standard deviation.");



    prm.declare_entry ("Regularization parameter",
                       "1.0",
                       dealii::Patterns::Double(),
                       "Weight given to the regularisation term.");

      prm.declare_entry ("Regularization method", "quadratic",
                       dealii::Patterns::Selection("quadratic|haar|sparse|TV"),
                       "Regularisation type");



      prm.declare_entry ("Penalty parameter",
                         "1.0",
                         dealii::Patterns::Double(),
                         "Weight given to the penalization of the quadratic fidelity term.");


      prm.declare_entry ("Image to reconstruct",
                         "",
                         dealii::Patterns::FileName(),
                         "path to the .tif image");

      // TODO
      prm.enter_subsection ("Input data");
      {

        // TODO
          prm.declare_entry ("Spatial dimension", "2",
                             dealii::Patterns::Integer(),
                             "Spatial Dimension of the tiff-stack. "
                             "Up to 3 the image is considered as something purely spatial.");

          // TODO : add flag for time-dependent data, e.g. 1+1, 2+1 or 3+1 dim


           // TODO
          prm.declare_entry("Anscombe", "false",
                            dealii::Patterns::Bool(),
                            "Do an Anscombe Transformation on input image");



      }
      prm.leave_subsection ();


      prm.enter_subsection ("Generation of unit test from real image");
      {

          prm.declare_entry ("Generate noisy image",
                             "true",
                             dealii::Patterns::Bool(),
                             "If set to false the input is treated as real data, if true input will be treated as test image and blurring and noise are added");

          prm.declare_entry ("Time seed",
                             "false",
                             dealii::Patterns::Bool(),
                             "If false simulated noise has a constant seed, if true the seed is taken from the clock");

      }
      prm.leave_subsection ();

      prm.enter_subsection("Dykstra");
      {
          prm.declare_entry ("N max Dykstra iterations",
                             "100",
                             dealii::Patterns::Integer(),"");
          prm.declare_entry ("Dykstra tolerance",
                             "1e-4",
                             dealii::Patterns::Double(),"Finish when |x_r - x_{r-1}| < tolerance");

          prm.declare_entry("Alpha quantile",
                            "0.9",
                            dealii::Patterns::Double(1e-8, 1.0-1e-8), "Principal parameter. Allowed range : (0,1). Default : 0.9. The higher the more accurate.");

      }
      prm.leave_subsection();


}

void AugmentedLagrangianParams::get(dealii::ParameterHandler &prm)
{

      this->psf_fwhm = prm.get_double("FWHM of PSF");

      this->std_dev_noise = prm.get_double("Std. dev. of Gaussian noise");

      this->regularization_strength = prm.get_double("Regularization parameter");

    this->penalty = prm.get_double("Penalty parameter");

     this->path_to_measured_image = prm.get ("Image to reconstruct");

    std::string temp_reg_type = prm.get("Regularization method");

    if (temp_reg_type == "quadratic") {
        this->reg_type = quadratic;
    }
    else if (temp_reg_type == "haar") {
        this->reg_type = Haar;
    }
    else if (temp_reg_type == "sparse") {
        this->reg_type = sparse;
    }
    else if (temp_reg_type == "TV") {
        this->reg_type = TV;
    }


    prm.enter_subsection("Generation of unit test from real image");
    {
        this->simulate = prm.get_bool("Generate noisy image");
        this->time_seed = prm.get_bool("Time seed");
    }
    prm.leave_subsection ();


    prm.enter_subsection("Dykstra");
    {
        this->n_max_dykstra_steps = prm.get_integer("N max Dykstra iterations");

        this->dykstra_Tol = prm.get_double("Dykstra tolerance");

        this->alpha_quantile = prm.get_double("Alpha quantile");

    }
    prm.leave_subsection();




     // TODO
    ImageIO img_io;


    img_io.read_image(image_as_read, path_to_measured_image);

    img_io.write_image("./output.tiff", image_as_read);



            // generate synthetic test data here

}



    // @sect3{Class: SimParams}
    //
    // This class contains parameters to steer the simulation.
    struct SimParams : public dealii::Subscriptor {

        SimParams() {}

        static void declare(dealii::ParameterHandler & prm);

        void get(dealii::ParameterHandler & prm);

        // A useful thing a default implementation can provide
        // is the location where potential results should be stored.
        // This becomes important when the parametric dependence of a
        // physical problem gets investigated.
        QDir run_dir;

        // The other thing is the directory where the parameters the
        // get logged which the program has actually read.
        QDir prm_log_dir;

        // Last but not least we provide a flag for the architecture.
        ParallelArch arch;

        // number of omp threads used for computation
        int n_omp_threads;

        // Since we want to avoid accidental copy operations
        // we make the copy CTor and the assignment operator private.
    private:
        SimParams (const SimParams & /*other*/) {}

        SimParams & operator = (const SimParams & /*other*/) { return *this; }
    };


    void SimParams::declare(dealii::ParameterHandler &prm)
    {
        prm.declare_entry("Run directory",
                          QString(
                              QString("./")
                              +
                              QString("results-"
                                      +
                                      QDateTime::currentDateTime().toString("ddd-yyyy-MM-dd/hh_mm_ss")
                                      ).remove(".")).toStdString(),
                          dealii::Patterns::Anything(),
                          "Specify a directory where the results of "
                          "the test are to be stored. This can be either "
                          "an absolute path or path relative to the directory "
                          "where the program has been started. The default is a "
                          "subdirectory called results-<date> where <date> will be replaced "
                          "by the date at which the program has been started. "
                          "this simplifies keeping the projects directory clean. "
                          "A prm file generated with the default values for the parameters gives an example.");

        prm.declare_entry("Use time tag",
                          "true",
                          dealii::Patterns::Bool(),
                          "If set to true, the results of the run will be stored in a subdirectory "
                          "of a given run directory which is named after the date and time at the "
                          "start of the program. This allows to group runs without "
                          "overwriting previous results.");

        prm.declare_entry("Architecture", "cuda",
                          dealii::Patterns::Selection("cpu|cuda"),
                          "Whether to use GPU acceleration for the SMRE. \n"
                          "cpu: do it on the CPU (implies OpenMP)\n"
                          "cuda: use your GPU (default).\n"
                          );

        prm.declare_entry ("Number of OpenMP threads",
                           "1",
                           dealii::Patterns::Integer(),
                           "");
    }


    void SimParams::get(dealii::ParameterHandler &prm)
    {
        std::string tmp_run_dir = prm.get("Run directory");

        bool use_time_tag = prm.get_bool("Use time tag");

        // Depending on the user's input the path is extended.
        // We can safely use a slash as QT will convert it into the
        // separator appropriate for the operating system.
        // Basically, this is a plain copy from step-27.
        if (use_time_tag)
            tmp_run_dir += (
                    QString("/"
                            +
                            QDateTime::currentDateTime().toString("ddd-yyyy-MM-dd/hh_mm_ss")
                            ).remove(".")).toStdString();

        // Get the run directory's absolute path and create it, if it does not exist.
        this->run_dir.setPath(tmp_run_dir.c_str() );
        this->run_dir.makeAbsolute();

        if(!this->run_dir.exists())
            this->run_dir.mkpath(".");


        std::string arch_string = prm.get("Architecture").c_str();

        std::cout << "chosen architecture : " << arch_string.c_str() << std::endl;

        if(arch_string == "cuda"){
            this->arch = gpu_cuda;
        } else if (arch_string == "cpu") {
            this->arch = cpu;
        } else {
            // Throw a fit if an invalid architecture has been provided
    #ifdef QT_NO_DEBUG
            AssertThrow(false,
                        dealii::StandardExceptions::ExcMessage("Unknown Architecture string. Must be \"cuda\", \"cpu\", or \"both\".") );
    #else
            Assert(false,
                   dealii::StandardExceptions::ExcMessage("Unknown Architecture string. Must be \"cuda\", \"cpu\", or \"both\".") );
    #endif
        }

        this->n_omp_threads = prm.get_integer("Number of OpenMP threads");

    }




#ifdef USE_OLDVERSION

    // @sect4{Class: SimParams}
    //
    // We collect the generic parameters of the simulation in a separate class.
    // For simplicity the method-specific ones are inherited from ADMMParams and their
    // declaration to and retrieval from the dealii::ParameterHandler is done from the
    // respective functions of the SimParams class.
    struct SimParams : public ADMMParams {

        //Input image name
        std::string input_image;

        //Output image name
        std::string out_imagename;

        //Wether or not to print control images during reports
        bool do_control;

        // One of the basic parameters of a simulation
        // is the location where potential results should be stored.
        // This becomes important when the parametric dependence of a
        // physical problem gets investigated.
        QDir run_dir;

        // We also need a directory where the parameter settings of a run should be logged.
        // This will always be @p run_dir + "/log".
        QDir prm_log_dir;

        static void declare(dealii::ParameterHandler &prm);

        void get(dealii::ParameterHandler &prm);
    };




    void SimParams::declare(dealii::ParameterHandler &prm)
    {
        ADMMParams::declare(prm);

        prm.enter_subsection ("Input");
        {

            prm.declare_entry ("Image to reconstruct",
                               "",
                               dealii::Patterns::FileName(),
                               "path to the .tif image");
        }
        prm.leave_subsection ();


        prm.enter_subsection ("Output");
        {
            prm.declare_entry ("control",
                               "false",
                               dealii::Patterns::Bool(),
                               "If true, preliminary results of the output image are saved to the run directory.");
            prm.declare_entry ("output image",
                               "control.tif",
                               dealii::Patterns::FileName(),
                               "Where should we put the ouput image? Will be a tiff image");

            prm.declare_entry("Run directory",
                              QString(
                                  QString("./")
                                  +
                                  QString("results-"
                                          +
                                          QDateTime::currentDateTime().toString("ddd-yyyy-MM-dd/hh_mm_ss")
                                          ).remove(".")).toStdString(),
                              dealii::Patterns::Anything(),
                              "Specify a directory where the results of "
                              "the test are to be stored. This can be either "
                              "an absolute path or path relative to the directory "
                              "where the program has been started. The default is "
                              "subdir called results-<date> where <date> will be replaced "
                              "by the date at which the program has been started. "
                              "this simplifies keeping the projects directory clean");
        }
        prm.leave_subsection ();
    }


    void SimParams::get(dealii::ParameterHandler &prm)
    {
        this->ADMMParams::get(prm);

        prm.enter_subsection ("Input");
        {
            this->input_image = prm.get ("Image to reconstruct");
        }
        prm.leave_subsection ();


        prm.enter_subsection("Output");
        {
            this->do_control = prm.get_bool("control");
            this->out_imagename = prm.get("output image");

            this->run_dir.setPath(prm.get("Run directory").c_str() );
            this->run_dir.makeAbsolute();

            if(!this->run_dir.exists())
                this->run_dir.mkpath(".");
        }
        prm.leave_subsection ();
    }

#endif


}


#endif // SIMPARAMS_H
