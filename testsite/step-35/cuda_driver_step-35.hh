//@sect3{File: cuda_driver_step-35.hh}
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

//Helper classes for CUDADriver class
#ifndef CUDA_DRIVER_STEP_35_HH
#define CUDA_DRIVER_STEP_35_HH

//CUDA
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>

//std stuff
#include <list>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <algorithm>    // std::min

//CUDA cuFFT
#include <cufft.h>

//Our stuff
#include "cuda_driver_step-35.h"
#include "cuda_kernel_wrapper_step-35.cu.h"
#include "patch.h"
#include "cuda_helper.h"
#include "preprocessor_directives.h"
#include "ADMMParams.h"

#include <deal.II/lac/vector.h>

namespace step35 {



}
#endif // CUDA_DRIVER_STEP_35_HH
