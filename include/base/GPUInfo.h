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


#ifndef GPUINFO_H
#define GPUINFO_H
#include <iostream>
namespace SciPAL {

// @sect3{Class: GPUInfo}
//
//! This an auxiliary structure which collects data related to the GPUs
//! available in the computer on which this program gets executed.
//! For details about the cuda<something> functions have a look at
//! the CUDA reference manual.
struct GPUInfo {

    int n_CUDA_devices;

    int current_device_id;

    cudaDeviceProp prop;

    //! @param d : ID of the device for which the information is to be retrieved.
    GPUInfo(int d=0)
        :
          current_device_id(d)
    {}


    //! This function allows to reset the device id. It automatically retrieves
    //! the technical infos about the new device.
    //! @param new_value : ID of the device for which the information is to be retrieved.
    void set_device_id(int new_value)
    {
        this->current_device_id = new_value;
        this->get();
    }

    // @sect4{Function: get}
    //
    //! This function must be called to retrieve the
    //! GPU-related information.
    void get()
    {
        static const int KB = 1024;
        static const int MB = KB*KB;

        //! Retrieve information about the currently selected GPU.
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

    //! To keep the compiler from automatically generating
    //! a copy constructor and an assignment operator we provide
    //! dummy implementations and declare them as private.
    //! In case one of them is needed the compiler will complain at compile-time
    //! and one can think about whether they are really needed, i.e. one has to
    //! review one's software design.
private:
    GPUInfo (const GPUInfo & /*other*/) {}

    GPUInfo & operator = (const GPUInfo & /*other*/) { return *this; }

};

} // namespace SciPAL END

#endif // GPUINFO_H
