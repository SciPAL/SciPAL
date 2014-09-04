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


#ifndef CUDATIMER_H
#define CUDATIMER_H

//!! STL header
#include <string>
#include <iostream>
#include <iomanip>

//!! plain C
//!!#include <stdio.h>

//!! utility library from CUDA SDK
//#include <cutil_inline.h>

#include <helper_cuda.h>
#include <helper_timer.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

/*! Wrapper class converting the timer functions
    provided by the cuda utility library into somethin OO-like.
    There should be plenty of use cases in the tutorial programs.
*/
class CUDATimer {

public:
        //!!! Create and start stop watch.
    CUDATimer() : __time(NULL), __is_running(false)
    {
        sdkCreateTimer(&__time);
        sdkStartTimer(&__time);
        //cutilCheckError(cutCreateTimer(&__time));
        reset();
    }

        //!!! Destroy the watch.
    ~CUDATimer()
    {
        sdkDeleteTimer(&__time);
        // ?? cutilCheckError(cutDeleteTimer(__time));
    }

	//!!! Stop the watch and record the time elapsed since reset() has been called.
    void stop()
    {
	if (__is_running)
	{
	    __is_running = false;
        sdkStopTimer(&__time);
            // ?? cutilCheckError(cutStopTimer(__time));
        __interval = .001 * sdkGetTimerValue(&__time);
           // ??  __interval = .001*cutGetTimerValue(__time);
	}
    }

        //!!! Turn back the clock.
    void reset()
    {
        if(__is_running)
            sdkResetTimer(&__time);
        // ??    cutilCheckError(cutStopTimer(__time));
        //!! cutilCheckError(cutDeleteTimer(__time));
        //!! cutilCheckError(cutCreateTimer(&__time));
        __interval = 0.;
        __is_running = true;
        // ?? cutilCheckError(cutStartTimer(__time));

    }

        //!!! Stop the watch and send the time that has elapsed so far to stdout.
    void print_elapsed (std::string message)
    {
	if (__is_running)
	{
        sdkStopTimer(&__time);
        // ?? cutilCheckError(cutStopTimer(__time));
        __interval = .001 * sdkGetTimerValue(&__time);
        // ??     __interval = .001*cutGetTimerValue(__time);
        }
		std::cout << message.c_str() <<  " " << __interval << std::endl;
		//!!std::cout << message.c_str() <<  std::fixed << " " << __interval << std::endl;
		//!!std::cout << message.c_str() << std::fixed << std::setprecision(9) <<" " << __interval << std::endl;
    }
	double really_elapsed() {  	if (__is_running)
		{
            sdkStopTimer(&__time);
            // ?? cutilCheckError(cutStopTimer(__time));
            __interval = .001 * sdkGetTimerValue(&__time);
            // ?? 	__interval = .001*cutGetTimerValue(__time);
			}
		return __interval; }


        //!!! This is a convenience function returning the time interval
		//!!! that has been measured by the last call of print_elapsed().
	double elapsed() {  
		return __interval; }

    double operator() () const { return __interval; }

private:
    StopWatchInterface *__time;
    volatile double       __interval;
    volatile bool         __is_running;
};

#endif //!! CUDATIMER_H
