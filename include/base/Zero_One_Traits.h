#ifndef ICH_WEIS_NICHT_H
#define ICH_WEIS_NICHT_H
namespace  SciPAL {


template<typename T> struct One
{
    typedef T Type;
public:
    T operator()(bool plus=true);
};



template<typename T>
inline T One<T>::operator ()(bool plus) {  return (plus ? 1. : -1); }

//!-------------------------------------------------------------------------

template<typename T> struct Zero
{
    typedef T Type;
public:
    T operator()();
};


template<>
inline cuComplex Zero<cuComplex>::operator ()()
{
    Type result;
    result.x = 0; result.y = 0.; return result;
}

template<>
inline cuDoubleComplex Zero<cuDoubleComplex>::operator ()()
{
    Type result;
    result.x = 0; result.y = 0.; return result;
}

//template<>
//std::complex<double> Zero<std::complex<double> >::operator ()()
//{
//    return std::complex<double>(0.0, 0.0);
//}

//template<>
//std::complex<float> Zero<std::complex<float> >::operator ()()
//{
//    return std::complex<float>(0.0, 0.0);
//}

template<typename T>
inline T Zero<T>::operator ()() {  return 0.; }


template<>
inline cuComplex One<cuComplex>::operator ()(bool plus)
{
    Type result; result.x = (plus ? 1. : -1.); result.y = 0.; return result;
}

template<>
inline cuDoubleComplex One<cuDoubleComplex>::operator ()(bool plus)
{
    Type result; result.x = (plus ? 1. : -1.); result.y = 0.; return result;
}

//template<>
//std::complex<double> One<std::complex<double> >::operator ()(bool plus)
//{
//    return std::complex<double>(plus ? 1. : -1., 0.0);
//}

//template<>
//std::complex<float> One<std::complex<float> >::operator ()(bool plus)
//{
//    return std::complex<float>(plus ? 1. : -1., 0.0);
//}
} //SciPAL end
#endif // ICH_WEIS_NICHT_H
