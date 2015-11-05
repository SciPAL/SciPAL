#ifndef FOREWARDDECLARATIONS_H
#define FOREWARDDECLARATIONS_H

#include <base/ParallelArch.h>
#include <complex>
#include <QString>
struct blas;
struct cublas;

namespace SciPAL
{
//!expr base class
template <class E> struct Expr;
template <typename T> struct CudaComplex;

//!base object for data
//!Shape Layouts
enum LAOType {matrix, vector};
//!specialization of layout
//!TODO use that to choose memory representation
enum sub_Layout{general, symm, banded, rowvector, columnvector};

template <typename> class ShapeData;
template <typename , typename, LAOType> class Shape;

//! aux structure to solve the missing Type problem in ShapeData
template<typename... T> struct GetMyType;


//!linear algebra objects = LAO
template<typename, typename> class Array;
template<typename, typename> class Matrix;
template<typename, typename> class SubMatrixView;
template<typename, typename> class Vector;
template<typename, typename, typename> class SubVectorView;
template<typename, typename, typename> class ColVectorView;
//    template<typename, typename> class ColVectorViewM;
//wrapper for simple numeric factors
template<typename T, typename BW> struct Literal;
template<typename T> struct DevLiteral;

template<typename T, ParallelArch arch> struct Kernels;

//!unary functions to modify LAOs
struct expr_transpose{};
struct expr_adjoint{};
struct expr_diag{};
template<typename T> struct Setter;


//!declare device expression types
template <typename _L, typename Operator, typename _R> struct DevBinaryExpr;
template <typename _L, typename Operation > struct DevUnaryExpr;

//!declare host expression types
template <typename _L, typename Operator, typename _R> struct BinaryExpr;
template <typename _L, typename Operation > struct UnaryExpr;

//!definitions of operations for DevBinaryExpr
template<typename OpTag, typename T> struct  binary;

//! This structure provides the functions for evaluating the expression tree in array arithmetic.
//!
template <typename T> struct ExprTree;
} // end namespace SciPAL

template<typename T> struct ComplexAdapter;

template <ParallelArch T> struct archTraits;
template<typename T, ParallelArch arch> struct PrecisionTraits;
template<typename T, ParallelArch arch> struct ComplexTraits;


typedef float2 cuComplex;
typedef double2 cuDoubleComplex;

template<typename> class FullMatrixAccessor;

void print_expr_info(QString expr_name);

#endif // FOREWARDDECLARATIONS_H

