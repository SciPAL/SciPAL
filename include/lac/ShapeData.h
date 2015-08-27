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


#ifndef  SCIPAL_SHAPEDATA_H
#define  SCIPAL_SHAPEDATA_H
#include <cuda_runtime_api.h>
namespace SciPAL {


//! Objects of this type have to be passed into CUDA kernels.
//! It collects the attributes which would have to be part of
//! the Shape class.
//! The purpose is to separate the shape data from the maintenace
//! functions like reinit() so that a CUDA kernel gets a plain-C structure
//! free of any member functions.
template <typename T>
struct ShapeData {

    //we will need this for the devexpr...
    //but it leads to ambiguities
    typedef T value_type;
    typedef const ShapeData<T> ConstHandle;

    //!  Pointer to the first element of the linear array containing the elements of the LAO.
    T * data_ptr;

    //! Number of rows. If this is set to 1 you get a row vector.
    size_t n_rows;

    //! Nnumber of columns. If this is set to 1 you get a column vector.
    size_t n_cols;

    //! If this value is, e.g., 3, then only every 3rd element of the original data array gets accessed.
    size_t stride;

    //! The length of a column in memory.
    //! This value can be used to optimize memory access of CUDA warps.
    //! We assume column-major storage in order to be compatible with CUBLAS.
    size_t leading_dim;


    size_t  r_begin_active;
    size_t  r_end_active;
    size_t  c_begin_active;
    size_t  c_end_active;
    size_t  n_elements_active;
    T *  view_begin;

    //! The default constructor generates an empty shape.
    __host__ __device__
    ShapeData()
    {
        this->data_ptr = 0;
        this->view_begin = 0;
        this->n_rows = 0;
        this->r_begin_active = 0;
        this->r_end_active = 0;
        this->n_cols = 0;
        this->c_begin_active = 0;
        this->c_end_active = 0;
        this->leading_dim = 0;
        this-> stride = 0;
    }
    __host__ __device__
    ShapeData& operator = (const ShapeData& other)
    {
        this->data_ptr = other.data_ptr;
        this->view_begin = other.view_begin;

        this->r_begin_active = other.r_begin_active;
        this->r_end_active = other.r_end_active;
        this->c_begin_active = other.c_begin_active;
        this->c_end_active = other.c_end_active;
        this->n_elements_active = other.n_elements_active;
        this->n_rows = other.n_rows;
        this->n_cols = other.n_cols;

        this->leading_dim = other.leading_dim;
        this->stride = other.stride;

        return *this;
    }
    __host__ __device__
    ShapeData(const ShapeData& other)
    {
        *this = other;
    }

//      op()(i,j)
//     op()(i)



};

}

#endif //  SCIPAL_SHAPEDATA_H
