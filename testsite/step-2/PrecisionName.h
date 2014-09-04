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


#ifndef PRECISIONNAME_H
#define PRECISIONNAME_H

#include <string>

// @sect3{struct: PrecisionName}
//
// Auxiliary-structure to write test precisions into the output files
template<typename T> struct PrecisionName
{
    static std::string name();
};

template<>
inline std::string PrecisionName<float>::name()  { return "float"; }

template<>
inline std::string PrecisionName<double>::name() { return "double"; }

template<typename T>
inline std::string PrecisionName<T>::name()      { return "unknown_T"; }


#endif // PRECISIONNAME_H
