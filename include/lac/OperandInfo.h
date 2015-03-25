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


#ifndef OPERANDINFO_H
#define OPERANDINFO_H

namespace SciPAL {

//forward declaration
template <class E> struct Expr;

template<typename _L,
         typename Operation,
         typename _R>
struct BinaryExpr;

//! Tag for addition. The tag classes are empty structures which get passed into the binary expressions in order to retrieve the operator at a later time during compilation.
struct plus {};


//! Tag for subtraction
struct minus {};


//! Tag for multiplication
struct mult {};


//! Tag for division
struct divide {};


//! Tag for pointwise multiplication
struct pmult {};



//! Tag for pointwise division
struct pdivide {};

}


#endif // OPERANDINFO_H
