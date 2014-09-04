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


#include <step-template/step-template.hh>

// @sect3{Function: main}
//
// As usual, the main function is pretty boring.
int main(int argc, char *argv[])
{
    using namespace steptemplate;

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
    SciPAL::GPUInfo gpu_info(DevNo);

    // Before we can instantiate the simulation we have to get the GPU info.
    gpu_info.get();

    // The simulation object gets the command line arguments for further processing,
    // especially reading the runtime parameters, and the information about the GPUs
    // which are available.
    MyFancySimulation simulation(argc, argv, gpu_info);

    simulation.run();
}
