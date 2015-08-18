#ifndef FOREWARDDECLARATIONS_H
#define FOREWARDDECLARATIONS_H

#include <base/ParallelArch.h>
#include <complex>

struct blas;
struct cublas;

namespace SciPAL
{
    template <typename T> struct CudaComplex;
    template<typename, typename> class Array;
    template<typename, typename> class Matrix;
    template<typename, typename> class SubMatrixView;
    template<typename, typename> class Vector;
    template<typename, typename> class VectorView;
    template<typename, typename> class ColVectorView;
    template<typename T, typename BW> struct Literal;
    template<typename T> struct DevLiteral;
    template<typename T> struct Setter;
    template<typename T, ParallelArch arch> struct Kernels;

    struct expr_transpose{};
    struct expr_adjoint{};
    struct expr_diag{};

    template <class E> struct Expr;

    //declare device expression types
    template <typename _L, typename Operator, typename _R> struct DevBinaryExpr;
    template <typename _L, typename Operation > struct DevUnaryExpr;

    //declare host expression types
    template <typename _L, typename Operator, typename _R> struct BinaryExpr;
    template <typename _L, typename Operation > struct UnaryExpr;

    //! aux structure to solve the missing Type problem in ShapeData
    template<typename T> struct GetMyType;

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

struct vmu_view;
struct vmu;

#endif // FOREWARDDECLARATIONS_H

