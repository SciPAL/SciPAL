// Dummy pragma to trick doxygen
#pragma
// <h2>step-28.cpp</h2>
// @sect3{File: step-28.cpp}

// TODO: secondary intro


// This is standard C++.
#include <iostream>
#include <vector>


// Driver class for GPU part of the program.
#include <step-28/cuda_driver_step-28.h>
#include <step-28/cuda_driver_step-28.hh>

#ifdef USE_STEP_32
#include <step-32/cuda_driver_step-32.h>
#include <step-32/cuda_driver_step-32.hh>
#endif
// For choosing on which device to run
#include  <step-27/Architecture.h>

// Include the actual dielectric relaxation spectroscopy simulation,
// courtesy of Stephan Kramer
#include <step-28/drs_simulation.hh>
#include <step-28/assemble_system.hh>
// #include <step-28/assemble_system_multithreaded.hh>

// For the timing
#include <timing.h>

// Include the SimulationHandler class
#include <step-27/simulation_handler.h>

#include <step-27/SimParams.h>

//deal.II includes
#include <deal.II/lac/vector.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>


//for Qt
#include <QtCore>
#include <QDebug>
#include <QVector>
#include <QTextStream>
#include <QThread>
#include <QMutex>


// @sect4{Function: main}
//
// The main function is pretty boring.
int main(int argc, char *argv[])
{
//    using namespace step28;

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

    // Then we set up the drs_simulation
    step27::SimulationManager< step28::DRSProblem<3>, step27::SimParams > simulation(argc, argv/*, gpu_info*/);

    // and run it
    return simulation.run();


}


// Initialize the static mutex
// shared by all instanced of the CellWorker class.
// This cannot be done in a header file.
// template<int dim> QMutex Step16::CellWorker<dim>::mutex;

