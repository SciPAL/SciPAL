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


#ifndef STEPTEMPLATE_HH
#define STEPTEMPLATE_HH

// This is standard C++.
#include <iostream>
#include <vector>

// Driver class for GPU part of th program.
#include <cuda_driver_step-template.h>
#include <cuda_driver_step-template.hh>

// The parameter class for your simulation.
#include <SimParams.h>


//deal.II includes
#include <deal.II/lac/vector.h>

// SciPal includes
//
// This an auxiliary structure which collects data related to the GPUs
// available in the computer on which this program gets executed.
// For details about the cuda<something> functions have a look at
// the CUDA reference manual.
#include <base/GPUInfo.h>

namespace steptemplate {

// @sect3{Class: MyFancySimulation}
//
// To make this test facility extendible, we implement
// a class for a simple user interface. Its primary tasks are
// - management of run-time parameters by a simple text-based parameter file
// - setting device parameters according to the user's parameters
// - preprocessing and output of results
class MyFancySimulation {

public:

    MyFancySimulation(int argc, char *argv[], SciPAL::GPUInfo &g);

    void run();

private:
    SciPAL::GPUInfo & gpuinfo;
    dealii::Vector<float> __solution;

protected:
    SimParams params;

    // A simulation has to know where to find its parameters.
    QString  prm_path;
};

}


// @sect4{Constructor: MyFancySimulation}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
steptemplate::MyFancySimulation::MyFancySimulation(int argc,
                                                   char *argv[],
                                                   SciPAL::GPUInfo &g)
    : gpuinfo(g)
{
    // Before setting up the simulation we
    // figure out how many GPUs are available
    cudaGetDeviceCount(&gpuinfo.n_CUDA_devices);
    std::cout
            << "N available CUDA devices : "
            << gpuinfo.n_CUDA_devices << std::endl;

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
    // and is supposed to be in a subdirectory prm of that directory,
    // where the program has been started.
    std::string master_prm_filename;
    if (argc == 1)
    {
        this->prm_path = (launch_dir.absolutePath().toStdString() + "/prm/").c_str();

        QFileInfo tmp(argv[0]);
        master_prm_filename  = tmp.baseName().toStdString();

        master_prm_filename = this->prm_path.toStdString() + master_prm_filename + ".prm";

        cwd.setPath("./prm");
        // At this point the sobdirectory prm may not yet exist. This is fixed after the
        // name of the parameter files are set up.
    }
    else
    {
        // Whatever gets passed as first command line argument is considered as path
        // to a parameter file.
        std::cout << "Given parameter file : " << argv[1] << std::endl;

        // We convert the sequence of characters into something more meaningful.
        QFileInfo tmp(argv[1]);

        // Before we proceed, let us figure out whether the given parameter file exists.
        // Note: If the file is a symlink that points to a non existing file,
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
        cwd.setPath(prm_path);
        cwd.makeAbsolute();
        master_prm_filename = tmp.fileName().toStdString();

        std::cout << "Parameter file path : "
                  << tmp.absolutePath().toStdString().c_str()
                  << std::endl;
    }

    std::cout << "Parameter file : " << master_prm_filename  << std::endl;

    // Before the parameter file can be read, we have to make sure that
    // its directory exists. In case of the default parameter file
    // the directory will be created.
    if (!cwd.exists() )
        launch_dir.mkpath( cwd.absolutePath() );

    QDir::setCurrent(cwd.absolutePath());

    prm_handler.read_input (master_prm_filename);

    QDir::setCurrent(launch_dir.absolutePath());

    this->params.get(prm_handler);

    // Create toplevel run directory
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
                                  (QFileInfo(master_prm_filename.c_str()).fileName()
                                   + ".log") ).toStdString();

    std::cout << "log file : " << master_log_file.c_str()
                                   << std::endl;

    std::ofstream log_master_prm( master_log_file.c_str() );
    prm_handler.print_parameters (log_master_prm,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    QDir::setCurrent(this->params.run_dir.absolutePath());

}



// @sect4{Function: run}
//
// Actual call to run function.
void steptemplate::MyFancySimulation::run()
{
        // Initialize dummy data on host ...
    float init_val = 4.17;

   __solution.reinit(100, true);
   __solution = init_val;

   // instantiate an object of the driver class
    steptemplate::CUDADriver testcase(&__solution(0), __solution.size() );

    // ... and run the computation on the GPU (or other dedicated parallel hardware).
    testcase.run();

    std::cout << "Done." << std::endl;
}

#endif // STEPTEMPLATE_HH
