//@sect3{File: preprocessor_directives.h}
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

#ifndef PREPROCESSOR_DIRECTIVES_H
#define PREPROCESSOR_DIRECTIVES_H

//@sect4{Preprocessor directives}
//Should we use double or single precision? \n
//Double precision is only supported for the exact Dykstra method,
//not for the approximative Dykstra method \n
//If DOUBLE_PRECISION is not defined we use floats

//#define DOUBLE_PRECISION


//If defined, print total dykstra kernel execution time to "kernel_time" for speedup plots \n

//#define TIME_KERNELS

#endif // PREPROCESSOR_DIRECTIVES_H
