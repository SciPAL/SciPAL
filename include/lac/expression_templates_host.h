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


#ifndef EXPRESSION_TEMPLATES_H
#define EXPRESSION_TEMPLATES_H

#include <iostream>

#include <base/PrecisionTraits.h>
#include <base/ArchTraits.h>

// Base class for expressions
#include <lac/Expr.h>
#include <base/ForewardDeclarations.h>
#include <lac/cublas_Matrix.h>
#include <lac/ShapeData.h>
#include <lac/Literal.h>

#include <lac/expression_template_helper.h>

#include <lac/OperandInfo.h>


//references:
// C++ Templates -The Complete Guide: Chapter 18, Chapter 16
// C++ Expression Templates - An Introduction to the Principles of
//Expression Templates
//implementation scheme
//http://www.angelikalanger.com/Articles/Cuj/ExpressionTemplates/TemplateBasedExpression.gif

#include <base/ForewardDeclarations.h>

namespace SciPAL {

// @sect3{Structs: binary and expressions}
//
// This class holds a binary part of the expression. This includes left and right
// operand and the operator.
template<typename _L,
         typename Operation,
         typename _R>
struct BinaryExpr
        :
        public Expr< BinaryExpr<_L, Operation, _R> >
{
    // Redefine template parameters as types, such they are accessible later from the outside
    typedef typename _L::Type L; //<- needed in apply function template
    typedef typename _R::Type R;
    typedef Operation OpTag;
    static const EType I_am = binE;
    // The following two typedefs are for the recursive use of the BinaryExpr class
    typedef BinaryExpr<_L, Operation, _R> Type;
    typedef DevBinaryExpr<typename _L::DevType, Operation,typename _R::DevType> DevType;

    typedef const Type& ConstHandle;

    // TODO: explain ovalue type: why deduced from _L and not _R?
    typedef typename _L::Type::value_type value_type;

    //we get the wrapper type from the _R because L can be a Literal, Literals do not have a blaswrapper.
    typedef typename _L::Type::blas_wrapper_type blas_wrapper_type;

    // For Literals the handle is a value and for matrices or vectors it is a reference.
    // these are the actual attributes of the binary expression
    typename _L::ConstHandle l;
    typename _R::ConstHandle r;

    Operation op;

    //The ~ operator is uesed here again on _l and _r to remove the outer Expr<>
    //envelope. This shortens the types of the generated expressions.
    BinaryExpr(const _L& _l, const _R& _r) : l(_l), r(_r) {}
    BinaryExpr(const BinaryExpr<_L, Operation, _R> & Ax): l(Ax.l), r(Ax.r){}


    operator const BinaryExpr<typename GetMyType<_L>::Type, Operation, typename GetMyType<_R>::Type>() const
    {
        return BinaryExpr<typename GetMyType<_L>::Type,
                Operation,
                typename GetMyType<_R>::Type>(l, r);
    }
};


// @sect3{Structs: Unary Expressions}
//
template<typename _L,
         typename Operation
         >
struct UnaryExpr
        :
        public Expr<UnaryExpr<_L,  Operation>
        >
{
//    typedef typename GetMyType<_L>::Type L; //<- needed in apply function template
    typedef _L L;

    typedef Operation OpTag;
    static const EType I_am = unE;
    typedef UnaryExpr<_L, Operation> Type;
    typedef DevUnaryExpr<typename _L::DevType, Operation> DevType;

    typedef const Type& ConstHandle;
    typedef typename L::value_type value_type;

    typename L::ConstHandle l;
    typedef typename _L::Type::blas_wrapper_type blas_wrapper_type;

    //I had to comment this in order for the Setter function to work.
    //Where does it break down? -jh
    //typedef typename _L::Type::blas_wrapper_type blas_wrapper_type;

    Operation op;

    UnaryExpr(const _L & _l) : l(~_l)
    {}
public:
    typename _L::ConstHandle get_l()
    {
        return this->l;
    }
};



// @sect3{General Operators}
//


template <typename T1, typename T2 >
inline
const BinaryExpr< T1, SciPAL::plus, T2 >
operator + (const Expr<T1> &e1, const Expr<T2> &e2)
{
    return BinaryExpr< T1, SciPAL::plus, T2 >(~e1, ~e2);
}

template <typename T1, typename T2 >
inline
const BinaryExpr< T1, SciPAL::minus, T2 >
operator - (const Expr<T1> &e1, const Expr<T2> &e2)
{
    return BinaryExpr< T1, SciPAL::minus, T2 >(~e1, ~e2);
}

//for the case that we want to multiply a real-valued literal with a complex LAO
template <typename T, typename BW, template <typename, typename> class LAO>
inline
const BinaryExpr<Literal<T, BW>, SciPAL::mult, LAO<T, BW> >
operator *(const typename PrecisionTraits<T, BW::arch>::NumberType e1,
           const  LAO<typename SciPAL::CudaComplex<typename PrecisionTraits<T, BW::arch>::NumberType>, BW>& e2)
{
    return BinaryExpr<Literal<T, BW>, SciPAL::mult, LAO<T, BW> >(~Literal<T, BW>(e1), ~e2);
}

template <typename T, typename BW, template <typename, typename> class LAO >
inline
const BinaryExpr<Literal<T, BW>, SciPAL::mult, LAO<T, BW> >
operator *(const typename LAO<T, BW>::value_type e1,
           const  LAO<T, BW>& e2)
{
    return BinaryExpr<Literal<T, BW>, SciPAL::mult, LAO<T, BW> >(~Literal<T, BW>(e1), ~e2);
}

template <typename T1, typename T2 >
inline
const BinaryExpr<T1, mult, T2>
operator * (const Expr<T1>& e1, const Expr<T2>& e2)
{
    return BinaryExpr<T1, SciPAL::mult, T2 >(~e1, ~e2);
}

//pointwise multiplication and divide operators
template <typename T1, typename T2 >
inline
const BinaryExpr<T1, SciPAL::pmult, T2 >
operator && (const Expr<T1>& e1, const Expr<T2>& e2)
{
    return BinaryExpr<T1, SciPAL::pmult, T2 >(~e1, ~e2);
}


template <typename T1, typename T2 >
inline
const BinaryExpr<T1, SciPAL::pdivide, T2 >
operator || (const Expr<T1>& e1, const Expr<T2>& e2)
{
    return BinaryExpr<T1, SciPAL::pdivide, T2 >(~e1, ~e2);
}


} // END namespace SciPAL

#endif // EXPRESSION_TEMPLATES_H
