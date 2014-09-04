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


#include <step-2/DeviceParams.h>

// @sect4{Function: declare}
//
void
step2::DeviceParams::declare(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Device flags");

    prm.declare_entry("Disable L1 caching for global memory",
                      "false",
                      dealii::Patterns::Bool(),
                      "This allows to dedicate the local L1 caches completely "
                      "to buffering register spilling effects.");

    prm.declare_entry("Maximize L1 cache",
                      "false",
                      dealii::Patterns::Bool(),
                      "This allows to rather use 48 KB as L1 cache.");

    prm.declare_entry("Wrap matrix into texture",
                      "true",
                      dealii::Patterns::Bool(),
                      "This allows to access the matrix via texture and "
                      "not by direct memory accesses. On older hardware "
                      "textures can exploit the 2D locality of matrix "
                      "entries. On Fermi- or Kepler-based GPUs this should not make "
                      "much sense.");

    prm.leave_subsection();
}


// @sect4{Function: get}
//
void
step2::DeviceParams::get(dealii::ParameterHandler &prm)
{
    prm.enter_subsection("Device flags");


    this->L2_caching_only   = prm.get_bool("Disable L1 caching for global memory");

    this->maximize_L1_cache = prm.get_bool("Maximize L1 cache");

    this->use_textures      = prm.get_bool("Wrap matrix into texture");


    prm.leave_subsection();
}
