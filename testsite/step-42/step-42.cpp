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


// This is standard C++.
#include <iostream>
#include <vector>

// Driver class for GPU part of th program.
#include <cuda_driver_step-42.h>
#include <cuda_driver_step-42.hh>

// The parameter class for your simulation.
#include <SimParams.h>


//deal.II includes
#include <deal.II/lac/vector.h>


namespace step42 {

// @sect3{Class: GPUInfo}
//
// This an auxiliary structure which collects data related to the GPUs
// available in the computer on which this program gets executed.
// For details about the cuda<something> functions have a look at
// the CUDA reference manual.
struct GPUInfo {

    int n_CUDA_devices;

    int current_device_id;

    cudaDeviceProp prop;


    GPUInfo(int DevNo)
        :
          current_device_id(DevNo)
    {

    }

    // @sect4{Function: get}
    //
    // this function must be called to retrieve the
    // GPU-related information.
    void get()
    {
        static const int KB = 1024;
        static const int MB = KB*KB;

        // Retrieve information about the currently selected GPU.
        std::cout << "current device ID : " << this->current_device_id << std::endl;

        cudaGetDeviceProperties(&prop, this->current_device_id);

        printf("Currently used GPU: %s \n",prop.name);
        printf("Compute Capability: %d.%d \n",prop.major,prop.minor);
        printf("ClockRate: %uMHz \n",prop.clockRate/1000);
        printf("Warpsize: %d \n",prop.warpSize);
        printf("Number of Multiprocessors: %d \n",prop.multiProcessorCount);

        printf("Shared Memory: %luKB\n",prop.sharedMemPerBlock/KB);
        printf("Constant Memory: %luKB \n",prop.totalConstMem/KB);
        printf("Global Memory: %luMB \n",prop.totalGlobalMem/MB);
        printf("the device %s concurrently copy memory between host and device while executing a kernel\n",
               (prop.deviceOverlap? "can": "cannot"));
    }

    // To keep the compiler from automatically generating
    // a copy constructor and an assignment operator we provide
    // dummy implementations and declare them as private.
    // In case one of them is needed the compiler will complain at compile-time
    // and one can think about whether they are really needed, i.e. one has to
    // review one's software design.
private:
    GPUInfo (const GPUInfo & /*other*/) {}

    GPUInfo & operator = (const GPUInfo & /*other*/) { return *this; }

};


// @sect3{Class: MyFancySimulation}
//
// To make this test facility extendible, we implement
// a class for a simple user interface. Its primary tasks are
// - management of run-time parameters by a simple text-based parameter file
// - setting device parameters according to the user's parameters
// - preprocessing and output of results
class MyFancySimulation {

public:

    MyFancySimulation(int argc, char *argv[], GPUInfo &g);

    void run();

private:
    GPUInfo & gpuinfo;

protected:
    SimParams params;
};

}


// @sect4{Constructor: MyFancySimulation}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
step42::MyFancySimulation::MyFancySimulation(int argc,
                                                   char *argv[],
                                                   step42::GPUInfo &g)
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
    std::string prm_filename;
    if (argc == 1)
    {
        prm_filename  = argv[0];
        prm_filename += ".prm";

        cwd.setPath("./prm");
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

            qFatal("Cannot proceed without proper path to paramter file");
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

    // Before the parameter file can be read, we have to make sure that
    // its directory exists. In case of the default parameter file
    // the directory will be created.
    if (!cwd.exists() )
        launch_dir.mkpath( cwd.absolutePath() );

    QDir::setCurrent(cwd.absolutePath());

    prm_handler.read_input (prm_filename);

    QDir::setCurrent(launch_dir.absolutePath());

    this->params.get(prm_handler);

    // Create toplevel run directory
    cwd.setPath(this->params.run_dir.absolutePath());

    // The following lets a directory make its own path.
    if (!cwd.exists())
        cwd.mkpath( "." );

    // Now, change to the run dir
    QDir::setCurrent(cwd.absolutePath());

    cwd.setPath("./log");
    cwd.makeAbsolute();
    if (!cwd.exists())
        cwd.mkpath(".");

    // Create the log directory and write what has been actually read
    // into log file. Basically, this is just another parameter file
    // and can thus be used again as input to another run after stripping the .log suffix.
    QDir::setCurrent(cwd.absolutePath());

    prm_filename += ".log";
    std::ofstream log_out_text(("./" + QString(prm_filename.c_str()).split("/").last()).toStdString().c_str());
    prm_handler.print_parameters (log_out_text,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    QDir::setCurrent(this->params.run_dir.absolutePath());
}



// @sect4{Function: run}
//
// Actual call to run function.
void step42::MyFancySimulation::run()
{

   // instantiate an object of the driver class
    step42::CUDADriver testcase;

    // ... and run the computation on the GPU (or other dedicated parallel hardware).
    testcase.gemm_tests();

//    testcase.gemv_tests();

//    testcase.complex_tests();

//    testcase.cusolver_demonstration();
    
//    testcase.feature_demonstration();

//    testcase.lin_combo();

    testcase.views();

//    testcase.stacks_of_LAOs();
//    testcase.operator_precedence();

    std::cout << "Done." << std::endl;
}



// @sect3{Function: main}
//
// As usual, the main function is pretty boring.
int main(int argc, char *argv[])
{
    using namespace step42;

    // At the beginning we figure out
    // how many CUDA devices are available.
    int n_CUDA_devices;

    cudaGetDeviceCount(&n_CUDA_devices);
    std::cout
            << "N available CUDA devices : "
            <<  n_CUDA_devices
            << std::endl;

    // This command is used to set the GPU on which
    // we want to run our computations.
    // For a list of GPUs execute nvidia-smi.
    int DevNo = 0;
    cudaSetDevice(DevNo);
    GPUInfo gpu_info(DevNo);

    // Before we can instantiate the simulation we have to get the GPU info.
    gpu_info.get();

    MyFancySimulation machma(argc, argv, gpu_info);

    machma.run();

}
