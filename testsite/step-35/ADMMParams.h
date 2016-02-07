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

#include <deal.II/lac/solver_control.h>

namespace step35 {
//@sect4{enum: regularisationType}
//@param haar regularization by haar wavelet sparsity
//@param sparse regularization by direct space sparsity
//@param quadratic Minimize 2-Norm
enum regularisationType {haar, sparse, quadratic, TV};

// @sect4{Class: ADMMParams}
//
// This class contains parameters read from control file, default: step-35.prm
class ADMMParams {

public:
    // First of all, we declare the most important parameter of the whole algorithm: How
    // probable it is that the residual is chi-squared distributed.
    double alpha_quantile;

    //Width of point spread function as double
    double psf_fwhm;

    //Stopping tolerance of ADMM
    // double tol;

    //First Lagrangian update parameter
    double alpha1;
    //Second Lagrangian update parameter
    double alpha2;
    //First constraint quadratic weight
    double rho1;
    //Second constraint quadratic weight
    double rho2;
    //Regularization parameter, bigger is smoother
    double reg_strength;
    //Wether or not to use time seeded random numbers
    bool time_seed;

    // Inverse stabilization parameter of the ADMM.
    double inv_gamma;

    // Maximum number of Dykstra iterations
    int n_max_dykstra_steps;

    //Maximum number of Heun steps, number of iteration steps to compute the reconstructed image for given noise.
    int n_max_heun_steps;

    // In each time step of the heun algorithm we adapt the time step this many times.
    int n_max_heun_dt_adaption_steps;

    // Threshold for L2-norm of solution increment in Dykstra.
    double dykstra_Tol;

    // Tolerance of the norm of the difference of two successive Heun steps.
    double heun_Tol;

    // Tolerance of the norm of the difference of two successive time step attempts within a Heun step.
    double dt_adaption_Tol;


    //Used for shifted indexing
    int ti,tj;
    //Report interval of ADMM
    int report_interval;

    // number of omp threads used for computation
    int n_omp_threads;

    //Wether or not to use approximative projections
    bool do_approx;
    //Wether or not to simulate noise on a test image
    bool simulate;
    //Dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    int dim;
    //Used in statistics generation
    int step;
    //Standard deviation of gaussian noise
    double gnoise;

    //Do an anscombe before and after
    bool anscombe;
    //Type of regularisation used in SMRE
    regularisationType regType;

    ADMMParams() {}


    static void declare(dealii::ParameterHandler &prm);

    void get(dealii::ParameterHandler &prm);

    dealii::SolverControl solver_control;

private:
    ADMMParams (const ADMMParams & /*other*/) {}

    ADMMParams & operator = (const ADMMParams & /*other*/) {
        return *this;
    }
};

//@sect5{Function: declare}
//@brief Declaration our parameters
void ADMMParams::declare(dealii::ParameterHandler &prm) {


    prm.enter_subsection ("Solver control");

    dealii::SolverControl::declare_parameters(prm);
    prm.enter_subsection("Heun");
    {
        prm.declare_entry ("N max time steps",
                           "1000",
                           dealii::Patterns::Integer(),"Maximum number of Heun steps, number of iteration steps to compute the reconstructed image for given noise.");

        prm.declare_entry ("N max time step adaptions",
                           "10",
                           dealii::Patterns::Integer(),"");

        prm.declare_entry ("Tolerance",
                           "1e-4",
                           dealii::Patterns::Double(),
                           "Tolerance of the norm of the difference of two successive Heun steps.");

        prm.declare_entry ("Time step adaption tolerance",
                           "1e-4",
                           dealii::Patterns::Double(),
                           "Tolerance of the norm of the difference of two successive time step attempts within a Heun step.");

    }
    prm.leave_subsection();

    prm.enter_subsection("Dykstra");
    {
        prm.declare_entry ("N max Dykstra iterations",
                           "100",
                           dealii::Patterns::Integer(),"");
        prm.declare_entry ("Dykstra tolerance",
                           "1e-4",
                           dealii::Patterns::Double(),"Finish when |x_r - x_{r-1}| < tolerance");
    }
    prm.leave_subsection();
    prm.leave_subsection ();


    prm.enter_subsection ("input data");
    {
        prm.declare_entry ("PSF FWHM",
                           "3",
                           dealii::Patterns::Double(),
                           "Full-width half-maximum of the PSF. In case of a PSF to narrow to be reasonably sampled "
                           "on the lattice of pixels we simply assume that it is a delta peak "
                           "and consider the convolution operator as identity.");
        prm.declare_entry ("regularization",
                           "1.0",
                           dealii::Patterns::Double(),
                           "intensity of regularisation");


        prm.declare_entry("Alpha quantile",
                          "0.9",
                          dealii::Patterns::Double(1e-8, 1.0-1e-8), "Principal parameter. Allowed range : (0,1). Default : 0.9. The higher the more accurate.");

        prm.declare_entry("Gamma",
                          ".1",
                          dealii::Patterns::Double(),
                          "Stabilization parameter. Similar to a pseudo-timestep. The smaller, the more stable but also the more slower the convergence." );



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
        prm.declare_entry ("Regularization method", "quadratic",
                           dealii::Patterns::Selection("quadratic|haar|sparse|TV"),
                           "Regularisation type");
        prm.declare_entry ("dim", "2",
                           dealii::Patterns::Integer(),
                           "Dimension of the images (if TIFF-stack process each slice individually if dim = 2)");
        prm.declare_entry("Anscombe", "false",
                          dealii::Patterns::Bool(),
                          "Do an Anscombe Transformation on input image");
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
        /*
        prm.declare_entry ("maximum iterations",
                           "10000",
                           dealii::Patterns::Integer(),"Set maximum number of iterations");
        prm.declare_entry ("tolerance",
                           "1e-3",
                           dealii::Patterns::Double(),"Finish when |x_r - x_{r-1}| < tolerance");
        */

        prm.declare_entry ("report interval",
                           "1",
                           dealii::Patterns::Integer(),
                           "Reporting progress in intervals of ... Iterations");

        prm.declare_entry ("Number of OpenMP threads",
                           "1",
                           dealii::Patterns::Integer(),
                           "");

    }
    prm.leave_subsection ();
}

//@sect5{Function: get}
//@brief Read the parameters from file or get default value
void ADMMParams::get(dealii::ParameterHandler &prm) {


    prm.enter_subsection ("Solver control");

    solver_control.parse_parameters(prm);
    prm.enter_subsection("Heun");
    {
        this->n_max_heun_steps = prm.get_integer("N max time steps");

        this->n_max_heun_dt_adaption_steps = prm.get_integer ("N max time step adaptions");

        this->heun_Tol = prm.get_double("Tolerance");

        this->dt_adaption_Tol = prm.get_double("Time step adaption tolerance");
    }
    prm.leave_subsection();

    prm.enter_subsection("Dykstra");
    {
        this->n_max_dykstra_steps = prm.get_integer("N max Dykstra iterations");

        this->dykstra_Tol = prm.get_double("Dykstra tolerance");
    }
    prm.leave_subsection();
    prm.leave_subsection ();

    prm.enter_subsection("input data");
    {
        this->alpha1 = prm.get_double("alpha1");
        this->alpha2 = prm.get_double("alpha2");

        this->alpha_quantile = prm.get_double("Alpha quantile");

        this->inv_gamma = 1./prm.get_double("Gamma");

        this->rho1 = prm.get_double("rho1");
        this->rho2 = prm.get_double("rho2");
        this->dim = prm.get_integer("dim");

        this->gnoise = prm.get_double("gaussian noise");


        this->psf_fwhm = prm.get_double("PSF FWHM");

        this->reg_strength = prm.get_double("regularization");

        this->anscombe = prm.get_bool("Anscombe");

        this->dim = prm.get_integer("dim");

        std::string temp_regType;
        temp_regType = prm.get("Regularization method");
        if (strcmp(temp_regType.c_str(),"quadratic") == 0 ) {
            this->regType = quadratic;
        }
        else if (strcmp(temp_regType.c_str(),"haar") == 0 ) {
            this->regType = haar;
        }
        else if (strcmp(temp_regType.c_str() ,"sparse") == 0 ) {
            this->regType = sparse;
        }
        else if (strcmp(temp_regType.c_str() ,"TV") == 0 ) {
            this->regType = TV;
        }
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
        // this->max_it = prm.get_integer("maximum iterations");

        this->do_approx = prm.get_bool("approx");
        this->n_omp_threads = prm.get_integer("Number of OpenMP threads");
        // this->tol = prm.get_double("tolerance");
    }
    prm.leave_subsection ();


}


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

} // namespace step35 END

#endif // SIMPARAMS_H
