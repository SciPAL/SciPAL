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


#ifndef STEPMINI_HH
#define STEPMINI_HH

// This is standard C++.
#include <iostream>

// The declaration of the interface to the CUDA-backend
// is contained in the following header.
#include <cuda.h>
#include <cuda_runtime_api.h>

// SciPal includes
//
// This an auxiliary structure which collects data related to the GPUs
// available in the computer on which this program gets executed.
// For details about the cuda<something> functions have a look at
// the CUDA reference manual.
#include <base/GPUInfo.h>
#include <lac/release/cublas_wrapper.hh>
#include <lac/release/blas_wrapper.hh>
#include <lac/development/cublas_Matrix.h>

namespace stepmini {

typedef double Number;
typedef cublas BW;

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
};

}


// @sect4{Constructor: MyFancySimulation}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
stepmini::MyFancySimulation::MyFancySimulation(int argc,
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

    BW::Init();
}



// @sect4{Function: run}
//
// Actual call to run function.
void stepmini::MyFancySimulation::run()
{
    // some dummy vectors
    const unsigned int n_rows = 4;
    const unsigned int n_cols = 3;
    const unsigned int n_elements = n_rows * n_cols;

    std::vector<Number>
            a(n_elements, 1.),
            b(n_elements, 2.),
            c(n_rows * n_rows, 1.23);

    for (unsigned int i = 0; i < b.size(); i++ )
        b[i] = i+1;

    SciPAL::Matrix<Number, BW>
            A(n_rows, n_cols, a),
            B(n_cols, n_rows, b),
            C(n_rows, n_rows, c);

    Number alpha = 1.;
    Number beta = 1.;

    std::cout << "A : " << std::endl;
    A.print();

    std::cout << "B : " << std::endl;
    B.print();

    std::cout << "C : " << std::endl;
    C.print();

    //gemm test
    std::cout << " ============ C = " << alpha << " * A * B + "  << beta  << " * C======" << std::endl;
    C = alpha * A * B + beta * C;
    std::cout << "C : " << std::endl; C.print();

    std::cout << "Done." << std::endl;
}

#endif // STEPMINI_HH
