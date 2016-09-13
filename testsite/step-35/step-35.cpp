//@sect3{File: step-35.cpp}
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

Copyright Stephan Kramer, Johannes Hagemann, Lutz Künneke, Jan Lebert 2014
*/

//@brief &nbsp; <!-- doxygen fix -->

//CUDA
#include <cuda_runtime.h>

// SciPAL
#include <base/GPUInfo.h>

#include <lac/blas++.h>

//omp
#ifdef HAS_OPENMP
#include <omp.h>
#endif

//Our stuff
#include <step-35/SimParams.h>
#include <step-35/ADMMParams.h>

// #include <step-35/step-35.hh>

#include <step-35/smre_problem.hh>

namespace step35 {

// @ sect3{Class: SimulationHandler}
// A simpler version of the step27::SimulationHandler class.
// It manages the parameter I/O and triggers the actual simulation.
template <typename NumberType, typename BW>
class SimulationManager {

public:
    // The CTor does quite a lot, actually. It first reads the parameters and then
    // initializes the GPU, i.e. figures out how devices there are and which one to use.
    SimulationManager(int argc, char *argv[]);

    void run_admm();

protected:
    void output(int iter, const AugmentedLagrangian<NumberType, BW>& sip) const;

    QString master_prm_filepath;
    QString prm_inverse_problem_filepath;
    QString prm_admm_filepath;


    SimParams params;

    AugmentedLagrangianParams cost_fctl_params;

    ADMMParams admm_params;

    // A simulation has to know where to find its parameters.
    QString  prm_path;
};

}



template <typename NumberType, typename BW>
void step35::SimulationManager<NumberType, BW>::output(int iter,
                                                       const AugmentedLagrangian<NumberType, BW>& sip) const
{
    TiffImage x_d_tmp(sip.params->image_as_read, false);

    if (iter==0) // (this->params.output_rec_image)
    {
        QString x_d_fname = QString ("rhs_d-%1.tif").arg(iter, 5, 10, QChar('0'));

        AssertThrow( sip.measured_data().size() > 0, dealii::ExcMessage("no rhs data available"));


       sip.measured_data().push_to(x_d_tmp.data);
       ImageIO::write_image(x_d_fname.toStdString(),
                            x_d_tmp // , 1.0
                            /*dummy value for Gaussian noise*/ // , params.anscombe
                            );
    }


    if (true) // (this->params.output_rec_image)
    {
        QString x_d_fname = QString ("x_d-%1.tif").arg(iter, 5, 10, QChar('0'));

        AssertThrow( sip.measured_data().size() > 0, dealii::ExcMessage("no rhs data available"));


       sip.signal_estimator().push_to(x_d_tmp.data);
       ImageIO::write_image(x_d_fname.toStdString(),
                            x_d_tmp // , 1.0
                            /*dummy value for Gaussian noise*/ // , params.anscombe
                            );
    }

    if (true)
    {
        QString x_d_fname = QString ("lag_mult_d-%1.tif").arg(iter, 5, 10, QChar('0'));

        sip.lag_mult().push_to(x_d_tmp.data);

        ImageIO::write_image(x_d_fname.toStdString(),
                             x_d_tmp // , 1.0
                             /*dummy value for Gaussian noise*/ // , params.anscombe
                             );

        // TO DO : output

    }

      if (true)
    {
        QString x_d_fname = QString ("e_d-%1.tif").arg(iter, 5, 10, QChar('0'));

        // TO DO : output

        sip.noise().push_to(x_d_tmp.data);

        ImageIO::write_image(x_d_fname.toStdString(),
                             x_d_tmp // , 1.0
                             /*dummy value for Gaussian noise*/ // , params.anscombe
                             );

      }

}



// @sect4{Function: main}
//
// The main function does the handling of the initial parameters
// by means of test file conforming to the dealii::ParameterHandler class.
// Once these basic parameters are read the program knows where to look up further
// parameter files.
// The path to the initial parameter file is either given
// as command line argument argv[1] or it is asumed that it is in a subdirectory prm in the
// current working directory.
//
int main(int argc, char **argv)
{
    typedef cublas BW;

    step35::SimulationManager<float, BW> sim_handler(argc, argv);

    sim_handler.run_admm();
}




template <typename NumberType, typename BW>
 step35::SimulationManager<NumberType, BW>::SimulationManager(int argc, char *argv[])
{
    // The management of the parameter I/O is done by
    dealii::ParameterHandler prm_handler;

    // To read parameters from a file we first have to declare them.
    // Basically, the parameter
    // file must have the same name as the binary. The extension has
    // to be ".prm". What has been read will be dumped into a log file.
    SimParams::declare(prm_handler);

    // Next, we get the current working directory ...
    QDir cwd = QDir::current();

    // i.e. where the program has been started.
    QDir launch_dir = cwd;

    std::cout << "CWD : " << launch_dir.absolutePath().toStdString().c_str() << std::endl;

    // By default, the master parameter file has the same name as the binary
    // and is supposed to be in a subdirectory prm of that directory,
    // where the program has been started.
    // As it may contain paths to the run directory and where results are to be stored
    // it has to be processed before all other parameter files.
    QString master_prm_file_basename;

    QString prog_name =  QFileInfo(argv[0]).baseName();

    QFileInfo tmp;

    // Depending on the number of cmd line arguments we either search for
    // prm file in the directory where the program was started
    if (argc == 1)
    {
        this->prm_path = (launch_dir.absolutePath() + "/prm/");

        tmp = QFileInfo(argv[0]);

        // Since we do not require that in the default case thr prm file actually exists
        // we still make the prm directory.
        cwd.setPath("./prm");
        if (!cwd.exists() )
            launch_dir.mkpath( cwd.absolutePath() );
        cwd.setPath(launch_dir.absolutePath());
    }
    // ... or in the filename pased as first cmd line argument. The path can be relative or absolute.
    else  if (argc == 2)
    {
        // Whatever gets passed as first command line argument is considered as path
        // to a parameter file.
        std::cout << "Given parameter file : " << argv[1] << std::endl;

        // We convert the sequence of characters into something more meaningful.
        tmp = QFileInfo(argv[1]);

        // Before we proceed, let us figure out whether the given parameter file exists.
        // Note: If the file is a symlink that points to a non-existing file,
        // false is returned as well.
        if(!tmp.exists())
        {
            std::cerr << "The following parameter file does not exist:\n"
                      << argv[1] << std::endl;

            qFatal("Cannot proceed without proper path to parameter file");
        }

        // Next, we subdivide the given filename into its path and filename
        // so that the corresponding subdirectories can be created.
        this->prm_path = tmp.absolutePath();
    }
    else
        std::cout << "Usage: " << master_prm_file_basename.toStdString().c_str() << " path-to-master-prm-file"
                  << "\nThe argument is optional and may be either a relative or absolute path." << std::endl;

    master_prm_file_basename = tmp.baseName();

    master_prm_filepath = this->prm_path + "/" + master_prm_file_basename + ".prm";

    this->prm_inverse_problem_filepath = (this->prm_path  + "/" + prog_name + "_inverse_problem.prm");
    this->prm_admm_filepath = (this->prm_path + "/" + prog_name + "_admm.prm");

    // For the convenience of the user we give an overview about where the prm files are sought for.
    std::cout << "Parameter file path used : "
              << this->prm_path.toStdString().c_str()
              << std::endl;

    std::cout << "Parameter file (master) : " << master_prm_filepath.toStdString().c_str()  << std::endl;
    std::cout << "Parameter file (problem data) : " << this->prm_inverse_problem_filepath.toStdString().c_str()  << std::endl;
    std::cout << "Parameter file (ADMM solver) : " << this->prm_admm_filepath.toStdString().c_str()  << std::endl;

    // Now, we can read them.
    prm_handler.read_input (master_prm_filepath.toStdString());


    this->params.get(prm_handler);

    // Create the toplevel run directory ...
    cwd.setPath(this->params.run_dir.absolutePath());

    // and the log dir as a child
    this->params.prm_log_dir.setPath(this->params.run_dir.absolutePath() + "/log/");

    // The following lets a directory make its own path.
    if (!cwd.exists())
        cwd.mkpath( "." );


    // After the run directory we create the log directory.
    if (!this->params.prm_log_dir.exists())
        this->params.prm_log_dir.mkpath(".");

    // Now, change to the run directory ...
    QDir::setCurrent(cwd.absolutePath());


    // and write what has been actually read
    // into log file. Basically, this is just another parameter file
    // and thus could be used again as input to another run after stripping the .log suffix.
        // When assembling the name of the log file including the absolute path we have to
        // strip any path-like components from the name of the master prm file.
    std::ofstream log_master_prm( (this->params.prm_log_dir.absolutePath() + QDir::separator()
                                   + master_prm_file_basename + ".prm.log").toStdString().c_str() );
    prm_handler.print_parameters (log_master_prm,
                                  dealii::ParameterHandler::Text);


    // Before we can create the next prm file we have to clear any previous content
    // from the prm handler.
    prm_handler.clear();

    // The parameters of the inverse problem and the solver are read in an analogous way.
    AugmentedLagrangianParams::declare(prm_handler);

    // Now, we can read them.
    prm_handler.read_input (this->prm_inverse_problem_filepath.toStdString());

    cost_fctl_params.get(prm_handler);

    std::ofstream log_problem_prm( (this->params.prm_log_dir.absolutePath() + QDir::separator()
                                   + QFileInfo(this->prm_inverse_problem_filepath).fileName() + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_problem_prm,
                                  dealii::ParameterHandler::Text);



    prm_handler.clear();

    ADMMParams::declare(prm_handler);

    prm_handler.read_input (this->prm_admm_filepath.toStdString());


    admm_params.get(prm_handler);

    std::ofstream log_admm_prm( (this->params.prm_log_dir.absolutePath() + QDir::separator()
                                   + QFileInfo(this->prm_admm_filepath).fileName() + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_admm_prm,
                                  dealii::ParameterHandler::Text);






    // After reading the parameter files we figure out
    // how many CUDA devices are available.
    int n_CUDA_devices;

    cudaGetDeviceCount(&n_CUDA_devices);
    std::cout
            << "N available CUDA devices : "
            <<  n_CUDA_devices
             << std::endl;

    // This command is used to set the GPU on which
    // we want to run our computations.
    // Currently, we just use the default. In a future verison we make this selectable from
    // the main prm file.
    int DevNo = 0;
    cudaSetDevice(DevNo);

    // The data related to the GPUs available in the computer
    // on which this program gets executed is collected in an auxiliary
    // structure provided by the library.
    // Currently, this information is just printed as information
    // about the GPUs.
    SciPAL::GPUInfo gpu_info(DevNo);

    // Before we can instantiate the simulation we have to get the GPU info.
    gpu_info.get();
}



 template <typename NumberType, typename BW>
 void step35::SimulationManager<NumberType, BW>::run_admm()
{

    // At the very beginning we initialize the blas library.
    BW::Init();

    AugmentedLagrangian<NumberType, BW> sip(this->cost_fctl_params);

    ADMMStepper<NumberType, BW> admm_stepper (this->admm_params);

     bool converged = false;

    int s = 0;
    while (!converged)
    {
        admm_stepper.step(sip);

        if (s%admm_params.report_interval == 0)
            this->output(s, sip);

        s++;
        // dummy
        if (s > admm_params.solver_control.max_steps())
        converged = true;
    }





    std::cout << "simulation goes here " << std::endl;

    // do something in between

    BW::Shutdown();
}


