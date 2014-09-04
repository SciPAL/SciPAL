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


#ifndef SPARSEMATRIXACCESSOR_H
#define SPARSEMATRIXACCESSOR_H

//!#include  <lac/sparse_matrix_mod.h>
#include <lac/cublas_BlancMatrix.h>

template<typename Number>
class SparseMatrixAccessor : protected BlancMatrix<Number>
{
public:
    //! Hier muss noch ein Copy-Konstruktor hin, der den Accessor aus der
    //! dealii-Matrix erzeugt

    SparseMatrixAccessor(const BlancMatrix<Number> &src);

    Number* vals;

    Number* get_values();

    unsigned int nnz();

};



#endif //! SPARSEMATRIXACCESSOR_H

