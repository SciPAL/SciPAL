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
#include <cstdio>
#include <ctime>

//Our stuff
#include <step-35/step-35.hh>
//Where the preprocessor directives are defined
#include "preprocessor_directives.h"

// @sect4{Function: main}
// We have to start somewhere
//argv[1] path to parameter file (optional)
//argv[2] Device ID of the GPU (optional, default: 0)
int main(int argc, char **argv) {
    cublas::Init();
    srand(time(NULL));
    using namespace step35;

//    ImageDecomposition_1024px iDecomp;
//    std::cout << iDecomp.print();


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
    if (argc == 3) {
        std::cout << "Using device " << argv[2] << std::endl;
        DevNo = atoi(argv[2]);
    }

    cudaSetDevice(DevNo);
    SciPAL::GPUInfo gpu_info(DevNo);

    // Before we can instantiate the simulation we have to get the GPU info.
    gpu_info.get();
    // The simulation object gets the command line arguments for further processing,
    // especially reading the runtime parameters, and the information about the GPUs
    // which are available.
#ifdef DOUBLE_PRECISION //set in preprocessor_directives.h
    ADMM<double> admm_instance(argc, argv, gpu_info);
#else
    ADMM<float> admm_instance(argc, argv, gpu_info);
#endif
    // Start the main routine
    admm_instance.run();

    cublas::Shutdown();
    return EXIT_SUCCESS;
}
