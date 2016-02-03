/*********************************************************************
 * Copyright © 2011-2014,
 * Marwan Abdellah: <abdellah.marwan@gmail.com>
 *
 * This library (cufftShift) is free software; you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 ********************************************************************/

#ifndef CUFFTSHIFT_1D_H
#define CUFFTSHIFT_1D_H

#include <cufftShiftShared.h>

namespace cufftShift
{
    void cufftShift_3D(cufftReal* input, cufftReal* output, int NX, int NY, int NZ);
    void cufftShift_3D(cufftDoubleReal* input, cufftDoubleReal* output, int NX, int NY, int NZ);
    void cufftShift_3D(cufftComplex* input, cufftComplex* output, int NX, int NY, int NZ);
    void cufftShift_3D(cufftDoubleComplex* input, cufftDoubleComplex* output, int NX, int NY, int NZ);

    void cufftShift_3D(cufftReal* data, int NX, int NY, int NZ);
    void cufftShift_3D(cufftDoubleReal* data, int NX, int NY, int NZ);
    void cufftShift_3D(cufftComplex* data, int NX, int NY, int NZ);
    void cufftShift_3D(cufftDoubleComplex* data, int NX, int NY, int NZ);
}

#endif // CUFFTSHIFT_1D_H
