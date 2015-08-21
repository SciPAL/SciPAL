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


// To outfox QTCreator's syntax highlighting and especially
// nvcc we put all cuda-related code into files with names
// ending on .cu.c and include them here.
// Then, in the project file we only have to take of one source
// file. This reduces the amount of maintenance.
#include <src/cuda/scipal_kernels.cu>

#include "cuda_kernel_step-4.cu.c"

//autoinstantiations from SciPAL
#include <step-4/autoInstantiations.h>
