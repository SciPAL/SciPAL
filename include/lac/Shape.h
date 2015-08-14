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


#ifndef SCIPAL_SHAPE_H
#define SCIPAL_SHAPE_H

#include <algorithm>
#include <lac/ShapeData.h>

namespace SciPAL {

//!Shape Layouts
enum LAOType {matrix, vector};
//!specialization of layout
//!TODO use that to choose memory representation
enum sub_Layout{general, symm, banded, rowvector, columnvector};

//! The <i>shape</i> of a linear algebra object (LAO) determines whether the LAO can be considered
//! as matrix or vector. By separating the shape from the data of a LAO it is easy to <i>reshape</i>
//! a LAO if necessary. The attributes are public for easy use in CUDA kernels.
template <typename T, LAOType L>
class Shape : public /*protected*/ ShapeData<T>
{
public:

    //! Data type of the LAO elements.
    typedef T value_type;

    //! is it a matrix or vector
    static const Layout layout = L;

    typedef Shape<T, L> Type;

    typedef  SciPAL::ShapeData<T> DevType;

    //! The default constructor generates an empty shape.
    Shape()
    {
        this->data_ptr = 0;
        this->n_rows = 0;
        this->n_cols = 0;
        this->leading_dim = 0;
        this-> stride = 0;
    }

    //! @param data : pointer to the first element of the linear array containing the elements of the LAO.
    //! @param n_rows : number of rows. If this is set to 1 you get a row vector.
    //! @param n_cols : number of columns. If this is set to 1 you get a column vector.
    //! @param stride : If this value is, e.g., 3, then only every 3rd element of the original data array gets accessed.
    //! @param leading_dim : the length of a column in memory. This value can be used to optimize memory access of CUDA warps. We assume column-major storage in order to be compatible with CUBLAS.
    Shape(T * data,
          size_t n_rows, size_t n_cols, size_t leading_dim,
          size_t stride = 1)
    {
        this->reinit(data, n_rows, n_cols, leading_dim, stride);
    }


    //! By assigning one shape to another we get two shapes looking at the same LAO.
    //! @param other: source shape to copy.
    Shape& operator = (/*const*/ Shape& other)
    {
#ifdef DEBUG
    std::cout << "line :" << __LINE__ << ", Shape<T>  " << __FUNCTION__<< "\n"  << std::endl;
#endif
        this->data_ptr = other.data_ptr;

        this->n_rows = other.n_rows;
        this->n_cols = other.n_cols;

        this->leading_dim = other.leading_dim;
        this->stride = other.stride;

        return *this;
    }

    //! Reset the shape data. This function typically gets called in the constructors
    //! and reinit() functions of the Matrix and Vector classes.
    void reinit(T * data,
                size_t  r_begin_active,
                size_t  r_end_active,
                size_t  c_begin_active,
                size_t  c_end_active,
                size_t leading_dim, size_t stride)
    {
        this->data_ptr = data;
        this->view_begin = data + (c_begin_active * leading_dim + r_begin_active);

        this->r_begin_active = r_begin_active;
        this->r_end_active = r_end_active;
        this->c_begin_active = c_begin_active;
        this->c_end_active = c_end_active;

        this->leading_dim = leading_dim;
        this->stride = stride;
    }

    //! Returns active rows
    size_t n_rows_active() const {
        return r_end_active - r_begin_active;
    }

    //! Returns active columns
    size_t n_cols_active() const {
        return c_end_active - c_begin_active;
    }

    //! Returns the number of elements, this is the same nuber the function n_elements of the
    //! array class returns. But note that we can have the situation where only the Shape
    //! iformation but the Array information is available.
    size_t size() const {
        return n_rows_active() * n_cols_active();
    }


protected:
    //! Swap the attributes of two shapes.
    //! This is reasonable only for LAOs which have the same total number of elements.
    //! Whether this is true has to checked externally.
    void swap(Shape& other)
    {
        std::swap(this->data_ptr, other.data_ptr);

        std::swap(this->n_rows, other.n_rows);
        std::swap(this->n_cols, other.n_cols);

        std::swap(this->leading_dim, other.leading_dim);
        std::swap(this->stride, other.stride);
    }

private:
    //! Copy constructor. The new object points to the same LAO as @p other.
    //! @param other: source shape to copy.
    Shape(const Shape& other)
    {
        *this = other;
    }


};

}
#endif // SCIPAL_SHAPE_H
