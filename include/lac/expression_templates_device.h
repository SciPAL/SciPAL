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


#ifndef DEVICEEXPR_H
#define DEVICEEXPR_H

#include <lac/expression_templates_host.h>

#include <lac/DevLiteral.h>

struct blas;

struct cublas;

namespace SciPAL {

//! aux structure to solve the missing Type problem in ShapeData
template<typename T> struct GetMyType;

template<typename T> struct GetMyType<SciPAL::ShapeData<T> >
{
typedef SciPAL::ShapeData<T> Type;
};

template<typename T> struct GetMyType {

    typedef typename T::Type Type;
};

//!definitions of operations for DevBinaryExpr
template<typename OpTag, typename T> struct  binary;

template<typename T>
struct  binary<plus, T>
{
    template <typename L, typename R>
    __host__ __device__ __forceinline__
    static  T eval(const L& a, const R& b) { return a + b; }
};

template<typename T>
struct  binary<mult, T>
{
    template <typename L, typename R>
    __host__ __device__  __forceinline__
    static  T eval(const L& a, const R& b) { return a * b; }
};

template<typename T>
struct  binary<pmult, T>
{
    template <typename L, typename R>
    __host__ __device__  __forceinline__
    static  T eval(const L& a, const R& b) { return a * b; }
};

template<typename T> struct
binary<minus, T>
{
    template <typename L, typename R>
    __host__ __device__  __forceinline__
    static  T eval(const L& a, const R& b) { return a - b; }
};

template<typename T> struct
binary<divide, T>
{
    template <typename L, typename R>
    __host__ __device__  __forceinline__
    static  T eval(const L& a, const R& b) { return a / b; }
    // Here, we could add a sanity check, whether we divide by zero
};

template<typename T> struct
binary<pdivide, T>
{
    template <typename L, typename R>
    __host__ __device__  __forceinline__
    static  T eval(const L& a, const R& b) { return a / b; }
    // Here, we could add a sanity check, whether we divide by zero
};


//do nothing if it is not +-*/
template<typename OpTag, typename T> struct  binary {};


//! This structure provides the functions for evaluating the expression tree in array arithmetic.
//!
template <typename T> struct ExprTree;

template <typename T>
struct ExprTree {

    //! Specialization for shapes of LAOs. Get element @p i from the data array.
    __host__
    __device__
    __forceinline__
    static T eval(const SciPAL::ShapeData<T> & Ax, int i)
    {
        // printf("%s ", __PRETTY_FUNCTION__); printf("shape data : Ex[%d] : %g\n", i, Ax.data_ptr[i] );
        return Ax.data_ptr[i];
    }

    //! Specialization for literals. Retrieve the value.
    template<template <typename> class Lit>
    __host__
    __device__
    __forceinline__
    static T eval(const Lit<T> & Ax, int i)
    {
       // printf("%s ", __PRETTY_FUNCTION__); printf("literal : Ex[%d] : %g\n", i, Ax.get_val() );
        return Ax.get_val();
    }

    //! This function evaluates a binary expression by evaluating its operands.
    //! Thus, this is recursive. The actual evaluation at this level is defered to the binary structure
    //! which depending on the specialization of @p Op either adds, multiplies,
    //! subtracts or divides the operands.
    template <typename L, typename Op, typename R>
    __host__
    __device__
    __forceinline__
    static T eval(const DevBinaryExpr<L/*e.g.: DevLiteral<T>*/, Op, R> & Ax, int i)
    {
        return Ax[i];
    }

    //! In the specialization for unary expressions the operator @p Op is applied to the elements of
    //! the operand @p Ax.l of the expression @p Ax.
    template <typename L, typename Op>
    __host__
    __device__
    __forceinline__
    static T eval (const DevUnaryExpr<L, Op> &Ax, int i)
    {
        return Ax.op(eval(Ax.l, i));
    }

};

//! DevBinary Expression
template <typename _L, typename Operator, typename _R>
struct DevBinaryExpr
{
    typedef Operator OpTag;

    // ShapeData can not provide "Type" due to the inheritance of Vector:Shape
    // and Matrix:Shape, this would cause double definition errors.
    typedef typename GetMyType<_L>::Type L;
    typedef typename GetMyType<_R>::Type R;

    typedef typename SciPAL::DevBinaryExpr<_L, Operator, _R> Type;

    typedef const Type ConstHandle;

    typedef typename L::value_type value_type;

    typename L::ConstHandle l;
    typename R::ConstHandle r;

    template<typename L_host, typename R_host>
    DevBinaryExpr(const BinaryExpr<L_host, Operator, R_host> &Ax) : l(Ax.l), r(Ax.r) {}

    __host__
    __device__
    __forceinline__
    DevBinaryExpr(const DevBinaryExpr<_L, Operator, _R> & Ax): l(Ax.l), r(Ax.r){}

    //! Evaluate the expression for the ith element.
    //! This operator is defined with shape-agnostic array arithmetic in mind.
    __host__
    __device__
    __forceinline__
    value_type operator [] (int i) const
    {
        // Evaluate operands before the expression at this level.
        // This should simplify the analysis of errors ar compile-time.
        // The 'E' typedef is needed, when the eval function is not static. Then @p e is an object of type @p E.
        typedef ExprTree<value_type> /*E;
        E*/ e;
        value_type left  = e::eval(l, i);
        value_type right = e::eval(r, i);

        return binary<Operator, value_type>::eval(left, right);
    }

};


//!DevUnary Expressions
template <typename _L, typename Operation >
struct DevUnaryExpr
{
    typedef Operation OpTag;
    Operation op;

    // ShapeData does not provide "Type" therefore, we have to add it.
    typedef typename GetMyType<_L>::Type L;

    typedef typename SciPAL::DevUnaryExpr<_L, Operation> Type;

    typedef const Type& ConstHandle;

    typedef typename L::value_type value_type;

    typename L::ConstHandle l;


    // L_host is needed to demote e.g. a Vector to its ShapeData
    template<typename L_host>
    __host__ __device__
    DevUnaryExpr (const UnaryExpr<L_host, Operation> &Ax)
        : l(Ax.l) {}
};


}//End namespace SciPAL
#endif // DEVICEEXPR_H
