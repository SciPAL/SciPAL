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


#ifndef blas_plus_plus_H
#define blas_plus_plus_H


//! deal.II-Komponenten
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>



//! Include implementations

#include <lac/blas_wrapper.hh>
#include <lac/cublas_wrapper.hh>

#include <lac/FullMatrixAccessor.h>

#include <lac/Array.h>

#include <lac/cublas_Vector.h>
#include <lac/cublas_SubVectorView.h>

#include <lac/cublas_Matrix.h>
#include <lac/cublas_SubMatrixView.h>


template<typename T, typename BW>
struct blas_pp {

    typedef T value_type;

    typedef ::FullMatrixAccessor<T> FullMatrixAccessor;

    typedef BW blas_wrapper_type;


    typedef SciPAL::Array<T, BW> Array;

    typedef SciPAL::Matrix<T, BW> Matrix;

    typedef SciPAL::SubMatrixView<T, BW> SubMatrix;

    typedef SciPAL::Vector<T, BW> Vector;


    typedef SciPAL::ColVectorView<T, Matrix> MatrixSubCol;

    typedef SciPAL::RowVectorView<T, Matrix> MatrixSubRow;

    typedef SciPAL::ColVectorView<T, Vector> SubColVector;

    typedef SciPAL::RowVectorView<T, Vector> SubRowVector;

    typedef SciPAL::SubVectorView<T, Vector>    SubVectorBase;
};


#endif //! blas_plus_plus_H
