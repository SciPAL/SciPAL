//@sect3{File: cuda_helper.h}
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

#ifndef CUDAHELPER_H
#define CUDAHELPER_H

#include<stdlib.h>

//@sect4{CUDA helper functions}
//CUDA helper functions  for error handling

//@sect5{checkCudaErrors}
//@brief print the proper CUDA error strings if a CUDA host call returns an error
#define checkCudaErrors(err)    __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors( cudaError err, const char *file, const int line ) {
    if( cudaSuccess != err) {
        printf("%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString(err));
        std::abort();
    }
}

//@sect5{getLastCudaError}
//@brief print the proper error string of a previous kernel call when calling cudaGetLastError
#define getLastCudaError(msg)   __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError( const char *errorMessage, const char *file, const int line ) {
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        printf("%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n", file, line, errorMessage, (int)err, cudaGetErrorString(err));
        std::abort();
    }
}
#endif // CUDAHELPER_H
