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


#ifndef PRECONDITIONJACOBIADAPTER_H
#define PRECONDITIONJACOBIADAPTER_H

#include <lac/precondition.h>
#include "cuda_driver_step-8.h"


//Purpose of this class: adapt dealii-ish preconditioners to the structure used in solve_qmrcgstab

template <typename Number,
typename MATRIX,
typename Base=dealii::PreconditionJacobi<MATRIX> >
        class PreconditionAdapter
            :
            public Base
{
public:

    PreconditionAdapter(const MATRIX &_A);

    PreconditionAdapter() : Base(), A(0) {}

    //omega : weight factor for diagonal entries
    void initialize(MATRIX &B, double omega = 1.2);

    template<typename VECTOR>
    void vmult (VECTOR &dst, const VECTOR &src) const;

    template<typename VECTOR>
    void postprocessor (VECTOR &dst, const VECTOR &src) const;

protected:
    const MATRIX* A;
};


template <typename Number,
typename MATRIX,
typename Base>
void PreconditionAdapter<Number, MATRIX, Base>::initialize(MATRIX &B, double omega)
{

    A = &B;
    Base::initialize(B, omega);

}

template <typename Number,
typename MATRIX,
typename Base>
PreconditionAdapter<Number, MATRIX, Base>::PreconditionAdapter(const MATRIX &_A)
    :
    Base(),
    A(&_A)
{

}



template <typename Number,
typename MATRIX,
typename Base>
template<typename VECTOR>
void PreconditionAdapter<Number, MATRIX, Base>::vmult (VECTOR &dst, const VECTOR &src) const
{
    VECTOR x;
    x = src;
    //! p__ = A*p_

    A->vmult(x,src);
    this->Base::vmult(dst,x);
}

template <typename Number,
typename MATRIX,
typename Base>
template<typename VECTOR>
void PreconditionAdapter<Number, MATRIX, Base>::postprocessor (VECTOR &dst, const VECTOR &src) const
{

}


#endif // PRECONDITIONJACOBIADAPTER_H
