#ifndef EXPRESSION_TEMPLATE_HELPER_H
#define EXPRESSION_TEMPLATE_HELPER_H

#include <base/ForewardDeclarations.h>
#include <QString>
#include <QStringList>
namespace SciPAL {

//! @sect3{struct: ExprChooser}
//! aux structure to determine the correct type of device expression.
//! If we have something like sin(a+b) it must be a DevUnaryExpr.
//! If we just have a+b it is DevBinaryExpr.
//! We can distinguish them by an static constant of type EType in the host
//! expression.
template<EType ExprType, typename E> struct ExprChooser;

template<typename E>
struct ExprChooser<unE, E>{

    typedef DevUnaryExpr<typename E::L::DevType,
    typename E::OpTag>      DevEType;

    typedef UnaryExpr<typename E::L::Type,
    typename E::OpTag>        HostEType;
};

template<typename E>
struct ExprChooser<binE, E>{

    typedef DevBinaryExpr<typename E::L::DevType,
    typename E::OpTag,
    typename E::R::DevType>   DevEType;

    typedef BinaryExpr<typename E::L::Type,
    typename E::OpTag,
    typename E::R::Type>   HostEType;

};

template<typename E>
struct ExprChooser<leafE, E>{

    typedef typename E::DevType  DevEType;

    typedef typename E::Type  HostEType;
};



//! @sect3{struct: GetMyType}
//! Struct to process type information in construction of BinE and UnE
//!
//! Special case 1: needed for leafs in device Expr
template<typename T0>
struct GetMyType<SciPAL::ShapeData<T0> >
{
    typedef SciPAL::ShapeData<T0> Type;
};

template<typename T, typename BW>
struct GetMyType<Matrix<T, BW> >
{
    typedef typename Matrix<T, BW>::MyShape Type;
};

template<typename T, typename BW>
struct GetMyType<Vector<T, BW> >
{
    typedef typename Vector<T, BW>::MyShape Type;
};


//! general case
template<typename T0>
struct GetMyType<T0> {
    typedef typename T0::Type Type;
};

//! LAO case
template<typename T0, typename T1>
struct GetMyType<T0, T1> {
    typedef SciPAL::ShapeData<T0> Type;
};

template<typename T0, typename... T>
struct GetMyType<T0, T...> {
    typedef typename GetMyType<T0, T...>::Type Type;
};



//! function to print information about expression evaluation
#ifdef DEBUG
void print_expr_info(QString expr_name)
{
QStringList names = expr_name.split("[with");
std::cout<<"evaluating expression :"<<names.at(1).toStdString()<<"\n"<<std::endl;
}
#endif
}

#endif // EXPRESSION_TEMPLATE_HELPER_H
