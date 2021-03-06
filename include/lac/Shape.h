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
#include <base/ForewardDeclarations.h>

namespace SciPAL {

//!Shape Layouts
enum LAOType {matrix, vector};
//!specialization of layout
//!TODO use that to choose memory representation
enum sub_Layout{general, symm, banded, rowvector, columnvector};

//! The <i>shape</i> of a linear algebra object (LAO) determines whether the LAO can be considered
//! as matrix or vector. By separating the shape from the data of a LAO it is easy to <i>reshape</i>
//! a LAO if necessary. The attributes are public for easy use in CUDA kernels.
template <typename T, typename BW, LAOType LT>
class Shape : public ShapeData<T>
{
    friend class Array<T, BW>;
public:

    //! Data type of the LAO elements.
    typedef T value_type;

    //! is it a matrix or vector
    static const LAOType layout = LT;

    typedef Shape<T, BW, LT> Type;
    //FIX ME unnecessary ? but needed for compiling
    typedef Type MyShape;

    typedef  SciPAL::ShapeData<T> DevType;

    typedef BW blas_wrapper_type;

    //! The default constructor generates an empty shape.
    Shape():ShapeData<T>()
    {}

    //! @param data : pointer to the first element of the linear array containing the elements of the LAO.
    //! @param n_rows : number of rows. If this is set to 1 you get a row vector.
    //! @param n_cols : number of columns. If this is set to 1 you get a column vector.
    //! @param stride : If this value is, e.g., 3, then only every 3rd element of the original data array gets accessed.
    //! @param leading_dim : the length of a column in memory. This value can be used to optimize memory access of CUDA warps. We assume column-major storage in order to be compatible with CUBLAS.

    Shape(size_t  r_begin_active,
          size_t  r_end_active,
          size_t  c_begin_active,
          size_t  c_end_active,
          size_t  stride = 1)
    {
        this->reinit(r_begin_active, r_end_active,
                     c_begin_active, c_end_active,
                     stride);
    }

    //! In contrast to the constructor above, this one reuses the memory
    //! information from @p other. But alters the elements the Shape's information.
    //! This is particular intresting when we build Views, which don't have own
    //! memory.
    template<LAOType LT2>
    Shape(const Shape<T, BW, LT2>& other,
          size_t  r_begin_active,
          size_t  r_end_active,
          size_t  c_begin_active,
          size_t  c_end_active,
          size_t  stride = 1)
    {
        this->reinit_attr(r_begin_active,
                          r_end_active,
                          c_begin_active,
                          c_end_active,
                          stride);

        this->data_ptr = other.data();
        this->leading_dim = other.leading_dim;

        this->view_begin = this->data_ptr +
                (this->c_begin_active * this->leading_dim + this->r_begin_active);
    }

    //! By assigning one shape to another we get two shapes looking at the same LAO.
    //! @param other: source shape to copy.
    template<LAOType LT2>
    Shape& operator = (const Shape<T, BW, LT2>& other)
    {
#ifdef DEBUG
        std::cout << "line :" << __LINE__ << ", Shape<T>  " << __FUNCTION__<< "\n"  << std::endl;
#endif
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

    void  reinit_attr(size_t  r_begin_active,
                      size_t  r_end_active,
                      size_t  c_begin_active,
                      size_t  c_end_active,
                      size_t  stride=1)
     {
         this->r_begin_active = r_begin_active;
         this->r_end_active = r_end_active;
         this->c_begin_active = c_begin_active;
         this->c_end_active = c_end_active;
         this->n_rows = this->n_rows_active();
         this->n_cols = this->n_cols_active();

        this->n_elements_active = this->size();

         this->stride = stride;

        //don't remove this, even if it looks as if it's done twice...
        this->view_begin = this->data_ptr +
                (this->c_begin_active * this->leading_dim + this->r_begin_active);
     }

    //! Reset the shape data. This function typically gets called in the constructors
    //! and reinit() functions of the Matrix and Vector classes.
    void reinit(size_t  r_begin_active,
                size_t  r_end_active,
                size_t  c_begin_active,
                size_t  c_end_active,
                size_t  stride)
    {
        //first get sizes than allocate memory
        storage.reinit(r_end_active-r_begin_active, c_end_active-c_begin_active);

        this->reinit_attr( r_begin_active,
                           r_end_active,
                           c_begin_active,
                           c_end_active,
                           stride);

        this->data_ptr = storage.val();
        this->leading_dim = storage.leading_dim();

        this->view_begin = this->data_ptr +
                (this->c_begin_active * this->leading_dim + this->r_begin_active);
    }

    ~Shape()
    {
        this->data_ptr = 0;
        this->view_begin = 0;
    }

    //! The following functions er used to access attributes of ShapeData
    //! Returns active rows
    size_t n_rows_active() const {
        return this->r_end_active - this->r_begin_active;
    }

    //! Returns active columns
    size_t n_cols_active() const {
        return this->c_end_active - this->c_begin_active;
    }

    size_t size() const {
        return n_rows_active() * n_cols_active();
    }

    //! Return ptr to data
    T* data()
    {
        return this->view_begin;
    }

    T* data() const
    {
        return this->view_begin;
    }

    //! Copy constructor. The new object points to the same LAO as @p other.
    //! @param other: source shape to copy.

    Shape(const Shape& other)
    {
        *this = other;
    }

protected:
    //! Swap the attributes of two shapes.
    //! This is reasonable only for LAOs which have the same total number of elements.
    //! Whether this is true has to checked externally.
    void swap(Shape& other)
    {
        std::swap(this->data_ptr, other.data_ptr);
        std::swap(this->view_begin, other.view_begin);

        std::swap(this->n_rows, other.n_rows);
        std::swap(this->r_begin_active, other.r_begin_active);
        std::swap(this->r_end_active, other.r_end_active);
        std::swap(this->n_cols, other.n_cols);
        std::swap(this->c_begin_active, other.c_begin_active);
        std::swap(this->c_end_active, other.c_end_active);
        std::swap(this->n_elements_active, other.n_elements_active);

        std::swap(this->leading_dim, other.leading_dim);
        std::swap(this->stride, other.stride);
    }

//    void swap_memory


private:
    Array<T, BW> storage;

};

}
#endif // SCIPAL_SHAPE_H
