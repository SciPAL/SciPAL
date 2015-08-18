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


#ifndef EXPR_H
#define EXPR_H

namespace SciPAL {

//differentiating a BinaryExpr and a UnaryExpr can be hard because of the different number of template parameters.
//this enum allows to distinguish them, in an easy manner.
enum EType {binE, unE, leafE};

// @sect3{Struct: Expr}
//
//This is a base class for a CRTP(Curious Recurring Template Programming)
//template <template<typename> class E>
template <class E>
struct Expr
{
    typedef E Type;

    const E &  operator~ () const {
        return static_cast<const E&>(*this);
    }
};

} //End namespace SciPAL
#endif // EXPR_H
