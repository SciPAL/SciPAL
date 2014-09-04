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


#ifndef DEVICEPARAMS_H
#define DEVICEPARAMS_H

#include <deal.II/base/parameter_handler.h>

namespace step2 {

// @sect3{struct: DeviceParams}
//
// Similar to the parameters of the simulation's user interface
// we collect all runtime parameters for tuning the GPU behavior in a
// dedicated structure.
struct DeviceParams {

    DeviceParams() {}


    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    bool L2_caching_only;

    bool maximize_L1_cache;

    bool use_textures;


private:
    DeviceParams (const DeviceParams & /*other*/) {}

    DeviceParams & operator = (const DeviceParams & /*other*/) { return *this; }
};

}

#endif // DEVICEPARAMS_H
