//@sect3{File: ADMMParams.h}
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

#ifndef SIMPARAMS_H
#define SIMPARAMS_H

//std
#include <iostream>

//Qt
#include <QDir>
#include <QTime>

// deal.II
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>

namespace step35 {
//@sect4{enum: regularisationType}
//@param haar regularization by haar wavelet sparsity
//@param sparse regularization by direct space sparsity
//@param quadratic Minimize 2-Norm
enum regularisationType {haar, sparse, quadratic};

// @sect4{Class: ADMMParams}
//
// This class contains parameters read from control file, default: step-35.prm
template<typename T>
class ADMMParams {

  public:
    //Width of point spread function
    int sigma;
    //Width of point spread function as double
    T sigmaf;
    //Stopping tolerance of ADMM
    T tol;
    //First Lagrangian update parameter
    T alpha1;
    //Second Lagrangian update parameter
    T alpha2;
    //First constraint quadratic weight
    T rho1;
    //Second constraint quadratic weight
    T rho2;
    //Regularization parameter, bigger is smoother
    T regInt;
    //Wether or not to use time seeded random numbers
    bool time_seed;
    //stabilization parameter of the ADMM
    T gamma_fac;
    //Maximum number of iterations
    int max_it;
    //Used for shifted indexing
    int ti,tj;
    //Report interval of ADMM
    int report_interval;
    //Wether or not to print control images during reports
    bool do_control;
    //Wether or not to use approximative projections
    bool do_approx;
    //Wether or not to simulate noise on a test image
    bool simulate;
    //Dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    int dim;
    //Used in statistics generation
    int step;
    //Standard deviation of gaussian noise
    T gnoise;
    //Input image name
    std::string imagename;
    //Output image name
    std::string out_imagename;
    //Do an anscombe before and after
    bool anscombe;
    //Type of regularisation used in SMRE
    regularisationType regType;
    //Stochastic or incomplete Dykstra
    bool stoch_dyk;
    //Width of rectangle
    int power_x;
    //Resolution
    int resol;
    //Decoposition rule for higher resolution
    bool both_halving;

    ADMMParams() {}


    static void declare(dealii::ParameterHandler &prm);

    void get(dealii::ParameterHandler &prm);

    // The only useful thing a default implementation can provide
    // is the location where potential results should be stored.
    // This becomes important when the parametric dependence of a
    // physical problem gets investigated.
    QDir run_dir;

    // Secondly, we need a directory where the parameter settings of a run should be logged.
    // This will always be @p run_dir + "/log".
    QDir prm_log_dir;
  private:
    ADMMParams (const ADMMParams & /*other*/) {}

    ADMMParams & operator = (const ADMMParams & /*other*/) {
        return *this;
    }
};

//@sect5{Function: declare}
//@brief Declaration our parameters
template<typename T>
void ADMMParams<T>::declare(dealii::ParameterHandler &prm) {
    prm.enter_subsection ("input data");
    {
        prm.declare_entry ("image",
                           "",
                           dealii::Patterns::FileName(),
                           "path to the .tif image");
        prm.declare_entry ("sigma",
                           "3",
                           dealii::Patterns::Double(),
                           "PSF spread");
        prm.declare_entry ("regularization",
                           "1.0",
                           dealii::Patterns::Double(),
                           "intensity of regularisation");
        prm.declare_entry ("rho1",
                           "6.192",
                           dealii::Patterns::Double(),
                           "stabilises first constraint");
        prm.declare_entry ("rho2",
                           "1.8",
                           dealii::Patterns::Double(),
                           "stabilises second constraint");
        prm.declare_entry ("alpha1",
                           "1.2",
                           dealii::Patterns::Double(),
                           "enforces first constraint");
        prm.declare_entry ("alpha2",
                           "0.12",
                           dealii::Patterns::Double(),
                           "enforces second constraint");
        prm.declare_entry ("gaussian noise",
                           "0",
                           dealii::Patterns::Double(),
                           "estimate of gaussian noise standard deviation. If simulate = true, will be used to add gaussian noise with standard deviation of ...");
        prm.declare_entry ("regType", "quadratic",
                           dealii::Patterns::Selection("quadratic|haar|quadratic"),
                           "Regularisation type -- TODO");
        prm.declare_entry ("dim", "2",
                           dealii::Patterns::Integer(),
                           "Dimension of the images (if TIFF-stack process each slice individually if dim = 2)");
        prm.declare_entry("anscombe", "false",
                          dealii::Patterns::Bool(),
                          "Do an Anscombe Transformation on input image");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("stoch dykstra");
    {
        prm.declare_entry ("stoch dyk","true",
                            dealii::Patterns::Bool(),
                            "stochastic dykstra on, else old incomplete dykstra is executed");
        prm.declare_entry ("power x","-1",
                           dealii::Patterns::Integer(),
                           "decomposition of image into rectangles with width 2^power_x "
                           "power_x = -1 means stochastically chosen decompositions");
        prm.declare_entry ("resol","-1",
                           dealii::Patterns::Integer(),
                           "random resolution, if resol = -1, all resolutions if resol = 0, "
                           "else consider only certain resoution depth (between 1 and power_x)");
        prm.declare_entry ("both halving","false",
                           dealii::Patterns::Bool(),
                           "for rectangles choose way of decomposing into chunks for higher resolution "
                           "both_halving = true means for next resolution to consider in algorithm the chunks are "
                           "created by halving both sides of the rectangle until one side is 1 pixel, then only the other is "
                           "divided by 2, else first the longest side is divided by 2 until a square is achieved, then both sides are halved "
                           "until their length is 1");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("simulate dataset from real image");
    {
        prm.declare_entry ("simulate",
                           "true",
                           dealii::Patterns::Bool(),
                           "If set to false the input is treated as real data, if true input will be treated as test image and blurring and noise are added");
        prm.declare_entry ("time seed",
                           "false",
                           dealii::Patterns::Bool(),
                           "If false simulated noise has a constant seed, if true the seed is taken from the clock");

    }
    prm.leave_subsection ();

    prm.enter_subsection ("program flow control");
    {
        prm.declare_entry ("approx",
                           "false",
                           dealii::Patterns::Bool(),
                           "do a small Dykstra in approximation");
        prm.declare_entry ("MRdepth",
                           "15",
                           dealii::Patterns::Integer(),
                           "largest patch edge length if not using small dykstra in approximation");
        prm.declare_entry ("maximum iterations",
                           "10000",
                           dealii::Patterns::Integer(),"Set maximum number of iterations");
        prm.declare_entry ("tolerance",
                           "1e-3",
                           dealii::Patterns::Double(),"Finish when |x_r - x_{r-1}| < tolerance");
        prm.declare_entry ("report interval",
                           "1",
                           dealii::Patterns::Integer(),
                           "Reporting progress in intervals of ... Iterations");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("output");
    {
        prm.declare_entry ("control",
                           "false",
                           dealii::Patterns::Bool(),
                           "save preliminary results of the output image");
        prm.declare_entry ("output image",
                           "control.tif",
                           dealii::Patterns::FileName(),
                           "where should we put the ouput image? Will be a tiff image");
        // The most generic parameter for a simulation is a string which indicates
        // where the simulation should store its results. By default, it contains
        // the current date and uses the time as subdirectory.
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

//@sect5{Function: get}
//@brief Read the parameters from file or get default value
template<typename T>
void ADMMParams<T>::get(dealii::ParameterHandler &prm) {

    prm.enter_subsection("input data");
    {
        this->alpha1 = prm.get_double("alpha1");
        this->alpha2 = prm.get_double("alpha2");
        this->rho1 = prm.get_double("rho1");
        this->rho2 = prm.get_double("rho2");
        this->dim = prm.get_integer("dim");

        this->gnoise = prm.get_double("gaussian noise");


        this->sigmaf = prm.get_double("sigma");
        this->sigmaf = 2.0*this->sigmaf; //the gaussian is calculated with the input sigma, its support has edge length 4 sigma
        this->sigma=std::ceil(this->sigmaf);

        this->imagename = prm.get("image");
        this->regInt = prm.get_double("regularization");

        this->anscombe = prm.get_bool("anscombe");

        this->dim = prm.get_integer("dim");

        std::string temp_regType;
        temp_regType = prm.get("regType");
        if (strcmp(temp_regType.c_str(),"quadratic") == 0 ) {
            this->regType = quadratic;
        }
        else if (strcmp(temp_regType.c_str(),"haar") == 0 ) {
            this->regType = haar;
        }
        else if (strcmp(temp_regType.c_str() ,"sparse") == 0 ) {
            this->regType = sparse;
        }
    }
    prm.leave_subsection ();

    prm.enter_subsection("stoch dykstra");
    {
        this->stoch_dyk = prm.get_bool("stoch dyk");
        this->power_x = prm.get_integer("power x");
        this->resol = prm.get_integer("resol");
        this->both_halving = prm.get_bool("both halving");
    }
prm.leave_subsection ();

    prm.enter_subsection("simulate dataset from real image");
    {
        this->simulate = prm.get_bool("simulate");
        this->time_seed = prm.get_bool("time seed");
    }
    prm.leave_subsection ();

    prm.enter_subsection("program flow control");
    {
        this->step = prm.get_integer("MRdepth");
        this->report_interval = prm.get_integer("report interval");
        this->max_it = prm.get_integer("maximum iterations");
        this->do_approx = prm.get_bool("approx");
        this->tol = prm.get_double("tolerance");
    }
    prm.leave_subsection ();

    prm.enter_subsection("output");
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
}

#endif // SIMPARAMS_H
