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


#ifndef EXPRESSION_TEMPLATE_H
#define EXPRESSION_TEMPLATE_H

#include <lac/Expr.h>


//!
//! @param L : Typ des linken Operanden
//! @param Op : Operation in dem Ausdruck
//! @param R : rechter Operand.
//!
template<typename L, typename Op, typename R>
struct X_read_read
        :
        public SciPAL::Expr< X_read_read<L, Op, R> > {

    const L & l;
    const R & r;

    X_read_read(const L & _l, const R & _r) : l(_l), r(_r) {}

    template<typename T>
    void apply(T & dst) const { Op::apply(dst, l, r); }
};


template<typename L, typename Op, typename R>
struct X_read_write
        :
        public SciPAL::Expr< X_read_write<L, Op, R> >{

    const L & l;
    R & r;

    X_read_write(const L & _l, R & _r) : l(_l), r(_r) {}

};


template<typename L, typename Op, typename R>
struct X_write_read
        :
      public SciPAL::Expr< X_write_read<L, Op, R> >{

    L & l;
    const R & r;

    X_write_read(L & _l, const R & _r) : l(_l), r(_r) {}

};


template<typename L, typename Op, typename R>
struct X_write_write
        :
       public SciPAL::Expr< X_write_write<L, Op, R> > {

    L & l;
    R & r;

    X_write_write(L & _l, R & _r) : l(_l), r(_r) {}

};

#endif //! EXPRESSION_TEMPLATE_H
