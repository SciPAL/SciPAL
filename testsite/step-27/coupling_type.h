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

Copyright S. C. Kramer, J. Hagemann 2010 - 2014, N. Wyderka, R. Czechowski 2014
*/

#ifndef COUPLING_TYPE_H
#define COUPLING_TYPE_H

// @sect4{enum: CouplingType}
//
// Enum to choose the FEM BEM coupling type. @p symmetric will choose
// the symmetric coupling as described in the introduction to this
// documentation. @p nonsymmetric selectes the non-symmetric coupling.
namespace step27 {
enum CouplingType { symmetric, nonsymmetric };
}


#endif // COUPLING_TYPE_H
