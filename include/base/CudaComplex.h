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


//Written by J. Hagemann, August/Septemper 2012
//actually a rebuild of std::complex for CUDA
//Update March 2013: Bugfixing, thanks to S. Meretzke for testing and pointing
//out some errors
#ifndef CUDACOMPLEX_H
#define CUDACOMPLEX_H
//SciPAL precision traits include
#include<cuComplex.h>
#include<base/PrecisionTraits.h>


// CudaComplex dosen't support mixed precision arithmetics.
namespace SciPAL {

template<typename T>
struct CudaComplex : public PrecisionTraits<T, gpu_cuda>::ComplexType
{
    typedef typename PrecisionTraits<T, gpu_cuda>::ComplexType ComplexType;
    typedef typename PrecisionTraits<T, gpu_cuda>::NumberType  NumberType;


    __host__ __device__ __forceinline__ CudaComplex(const NumberType re = NumberType(),
                                    const NumberType im = NumberType())
     { this->x = re, this->y = im;}

    __host__ __device__ __forceinline__ CudaComplex(const CudaComplex& a){*this = a;}

    __host__ __device__ __forceinline__ CudaComplex(const ComplexType& a)
    {this->x = a.x; this->y = a.y;}

    __host__ __device__ __forceinline__ CudaComplex(const std::complex<T>& a)
    {this->x = std::real(a); this->y = std::imag(a);}


    //return real and imaginary part
    __host__ __device__ __forceinline__ NumberType& real()
    {return this->x;}

    __host__ __device__ __forceinline__ NumberType& imag()
    {return this->y;}

    __host__ __device__ __forceinline__ const NumberType& real() const
    {return this->x;}

    __host__ __device__ __forceinline__ const NumberType& imag() const
    {return this->y;}

    //set real and imaginary part
    __host__ __device__  __forceinline__ void real(NumberType val)
    { this->x = val; }

    __host__ __device__  __forceinline__ void imag(NumberType val)
    { this->y = val; }

    //add value to real part
    __host__ __device__  __forceinline__ CudaComplex&
    operator+=(const NumberType& val)
    {
        this->x += val;
        return *this;
    }

    //substract...
    __host__ __device__ __forceinline__ CudaComplex&
    operator-=(const NumberType& val)
    {
        this->x -= val;
        return *this;
    }

    //multiply
    __host__ __device__ __forceinline__ CudaComplex& operator*=(const NumberType t)
    {
        this->x *= t;
        this->y *= t;
        return *this;
    }

    //divide
    __host__ __device__  __forceinline__ CudaComplex& operator/=(const NumberType t)
    {
        this->x /= t;
        this->y /= t;
        return *this;
    }

    //divide
    __host__ __device__
    __forceinline__
    CudaComplex& operator()(const NumberType re, const NumberType im)
    {
        this->x = re;
        this->y = im;
        return *this;
    }

    //a few declarations
    __host__ __device__ /*__forceinline__*/
                CudaComplex& operator=(const CudaComplex& a);
    __host__ __device__ __forceinline__ CudaComplex&
                operator-=(const CudaComplex& val);
    __host__ __device__ __forceinline__ CudaComplex&
                operator/=(const CudaComplex& val);
    __host__ __device__ __forceinline__ CudaComplex&
                operator*=(const CudaComplex& val);
    __host__ __device__ __forceinline__ CudaComplex<T>&
                operator+=(const CudaComplex<T>& val);

    operator std::complex<NumberType>() {
            return std::complex<NumberType>(this->x, this->y);
        }

private:

}; //end struct



__host__ __device__ static __forceinline__
CudaComplex<double> conj_impl(const CudaComplex<double>& a)
{CudaComplex<double> tmp; tmp/*.Number*/ = cuConj(a/*.Number*/); return tmp;}

__host__ __device__ static __forceinline__
CudaComplex<float> conj_impl(const CudaComplex<float>& a)
{CudaComplex<float> tmp; tmp/*.Number*/ = cuConjf(a/*.Number*/); return tmp;}


template<typename T>
__host__ __device__ __forceinline__
typename PrecisionTraits<T, gpu_cuda>::ComplexType toNumberType2(CudaComplex<T> a)
{
    typename PrecisionTraits<T, gpu_cuda>::ComplexType tmp;
    tmp.x = a.real();
    tmp.y = a.imag();
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
typename PrecisionTraits<T, gpu_cuda>::ComplexType toNumberType2(T a)
{
    typename PrecisionTraits<T, gpu_cuda>::ComplexType tmp;
    tmp.x = a;
    tmp.y = T(0.0);
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
typename PrecisionTraits<T, gpu_cuda>::ComplexType toNumberType2(T a, T b)
{
    typename PrecisionTraits<T, gpu_cuda>::ComplexType tmp;
    tmp.x = a;
    tmp.y = b;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>& CudaComplex<T>::operator=(const CudaComplex<T>& a)
{
    this->real() = a.real();
    this->imag() = a.imag();
    return *this;
}

template<typename T>
__host__ __device__ __forceinline__
 CudaComplex<T>& CudaComplex<T>::operator+=(const CudaComplex<T>& __z)
  {
    this->real() += __z.real();
    this->imag() += __z.imag();
    return *this;
  }

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>& CudaComplex<T>::operator-=(const CudaComplex<T>& __z)
  {
    this->real() -= __z.real();
    this->imag() -= __z.imag();
    return *this;
  }

// This is a grammar school implementation.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>& CudaComplex<T>::operator*=(const CudaComplex<T>& __z)
  {
      const T __r = this->real() * __z.real() - this->imag() * __z.imag();
      this->imag() = this->real() * __z.imag() + this->imag() * __z.real();
      this->real() = __r;
      return *this;
  }

// This is a grammar school implementation.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>& CudaComplex<T>::operator/=(const CudaComplex<T>& __z)
  {

      T s = ::abs(__z.real()) + ::abs(__z.imag());
      T oos = T(1.0) / s;
      T ars = this->real() * oos;
      T ais = this->imag() * oos;
      T brs = __z.real() * oos;
      T bis = __z.imag() * oos;
      s = (brs * brs) + (bis * bis);
      oos = T(1.0) / s;
      this->real() = ((ars * brs) + (ais * bis)) * oos;
      this->imag() = ((ais * brs) - (ars * bis)) * oos;
      return *this;
  }

//Operators:
//operator + for different cases
template<typename T>
__host__ __device__  __forceinline__
CudaComplex<T> operator+(const CudaComplex<T>& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(a);
    tmp += val;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator+(const CudaComplex<T>& a, const T& val)
{
    CudaComplex<T> tmp(a);
    tmp += val;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator+( const T& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(a);
    tmp += val;
    return tmp;
}


//operator - for different cases
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator-( const CudaComplex<T>& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(val);
    tmp -= a;
    return tmp;
}


template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator-(const CudaComplex<T>& a, const T& val)
{
    CudaComplex<T> tmp(a);
    tmp -= val;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator-( const T& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(val, -a.imag());
    tmp -= a.real();
    return tmp;
}

//operator * for different cases
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator*( const CudaComplex<T>& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(a);
    tmp *= val;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator*(const CudaComplex<T>& a, const T& val)
{
    CudaComplex<T> tmp(a);
    tmp *= val;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator*( const T& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(a);
    tmp *= val;
    return tmp;
}

//operator / for different cases
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator/( const CudaComplex<T>& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(val);
    tmp /= a;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator/(const CudaComplex<T>& val, const T& a)
{
    CudaComplex<T> tmp(val);
    tmp /= a;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator/( const T& val, const CudaComplex<T>& a)
{
    CudaComplex<T> tmp(val, 0.0);
    tmp /= a;
    return tmp;
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>
conj(const CudaComplex<T>& __z)
{ return conj_impl(__z);}

// abs(__z):  Returns the magnitude of __z.


__host__ __device__ __forceinline__
float abs_impl(const CudaComplex<float>& a )
{
    return cuCabsf(a/*.Number*/);
}

__host__ __device__ __forceinline__
double abs_impl(const CudaComplex<double>& a )
{
    return sqrt(a.real()*a.real() + a.imag()*a.imag());
}

__host__ __device__ __forceinline__
float abs_impl(const float a )
{
    return fabsf(a);
}

__host__ __device__ __forceinline__
double abs_impl(/*const*/ double a )
{
    return fabs(a);
}


template<typename T>
__host__ __device__ __forceinline__
T abs(const CudaComplex<T>& __z)
{
    return  abs_impl(__z);
}

//template<typename T>
//__host__ __device__ __forceinline__
//T abs(const T& __z)
//{
//    return  abs_impl(__z);
//}

// arg(__z): Returns the phase angle of __z.
template<typename T>
__host__ __device__ static __forceinline__
T arg(const CudaComplex<T>& __z)
{ return  atan2(__z.imag(), __z.real()); }

template<typename T>
__host__ __device__ __forceinline__
T norm(const CudaComplex<T>& __z)
{
    return abs(__z)*abs(__z);
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> polar(const T& __rho, const T& __theta)
{ return CudaComplex<T>(__rho * std::cos(__theta), __rho * std::sin(__theta)); }

//Transcendentals
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> cos(const CudaComplex<T>& __z)
{
    const T __x = __z.real();
    const T __y = __z.imag();
    return CudaComplex<T>(cos(__x) * cosh(__y), -sin(__x) * sinh(__y));
}
// cosh(__z): Returns the hyperbolic cosine of __z.
template<typename T>
__host__ __device__  __forceinline__
T cosh(const T __z)
{
    return ::cosh(__z);
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> cosh(const CudaComplex<T>& __z)
{
    const T __x = __z.real();
    const T __y = __z.imag();
    return CudaComplex<T>(cosh(__x) * cos(__y), sinh(__x) * sin(__y));
}
// exp(__z): Returns the complex base e exponential of x
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> exp(const CudaComplex<T>& __z)
{ return polar(::exp(__z.real()), __z.imag()); }

// log(__z): Returns the natural complex logarithm of __z.
//           The branch cut is along the negative axis.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> log(const CudaComplex<T>& __z)
{ return CudaComplex<T>(log(abs(__z)), arg(__z)); }

template<typename T>
__host__ __device__  __forceinline__
CudaComplex<T> sin(const CudaComplex<T>& __z)
{
    const T __x = __z.real();
    const T __y = __z.imag();
    return CudaComplex<T>(::sin(__x) * ::cosh(__y), ::cos(__x) * ::sinh(__y));
}

// sinh(__z): Returns the hyperbolic sine of __z.
template<typename T>
__host__ __device__  __forceinline__
T sinh(const T __z)
{
    return ::sinh(__z);
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> sinh(const CudaComplex<T>& __z)
{
    const T __x = __z.real();
    const T __y = __z.imag();
    return CudaComplex<T>(::sinh(__x) * ::cos(__y), ::cosh(__x) * ::sin(__y));
}

// sqrt(__z): Returns the complex square root of __z.
//            The branch cut is on the negative axis.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> sqrt(const CudaComplex<T>&__z)
{
    T __x = __z.real();
    T __y = __z.imag();

    if (__x == T())
    {
        T __t = ::sqrt(fabs(__y) / 2); // TO DO: fix
        return CudaComplex<T>(__t, __y < T() ? -__t : __t);
    }
    else
    {
        T __t = ::sqrt(2 * (abs(__z) + fabs(__x)));  // TO DO: fix
        T __u = __t / 2;
        return __x > T()
                ? CudaComplex<T>(__u, __y / __t)
                : CudaComplex<T>(fabs(__y) / __t, __y < T() ? -__u : __u);  // TO DO: fix
    }
}

// tan(__z):  Return the complex tangent of __z.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> tan(const CudaComplex<T>& __z)
{ return sin(__z) / cos(__z); }

// tanh(__z): Returns the hyperbolic tangent of __z.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> tanh(const CudaComplex<T>& __z)
{ return sinh(__z) / cosh(__z); }

// pow(__x, __y): Returns the complex power base of __x
//                raised to the __y-th power.  The branch
//                cut is on the negative axis.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> pow(const CudaComplex<T>& __z, int __n)
{
    return __n < 0
            ? T(1.0)/pow(__z, -__n)
            : pow(__z, __n);
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> pow(const CudaComplex<T>& __x, const T& __y)
{
    if (__x.imag() == T() && __x.real() > T())
        return pow(__x.real(), __y);

    CudaComplex<T> __t = log(__x);
    return polar(exp(__y * __t.real()), __y * __t.imag());
}

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> pow(const CudaComplex<T>& __x, const CudaComplex<T>& __y)
{ return __x == T() ? T() : exp(__y * log(__x)); }

template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> pow(const T& __x, const CudaComplex<T>& __y)
{
    return __x > T() ? polar(pow(__x, __y.real()),
                             __y.imag() * log(__x))
                     : pow(CudaComplex<T>(__x), __y);
}


///  Return @a x.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T> operator+(const CudaComplex<T>& __x)
{ return __x; }

///  Return complex negation of @a x.
template<typename T>
__host__ __device__ __forceinline__
CudaComplex<T>
operator-(const CudaComplex<T>& __x)
{  return CudaComplex<T>(-__x.real(), -__x.imag()); }

//Comparisonoperators

//@{
///  Return true if @a x is equal to @a y.
template<typename T>
__host__ __device__ __forceinline__ bool
operator==(const CudaComplex<T>& __x, const CudaComplex<T>& __y)
{ return __x.real() == __y.real() && __x.imag() == __y.imag(); }

template<typename T>
__host__ __device__ __forceinline__ bool
operator==(const CudaComplex<T>& __x, const T& __y)
{ return __x.real() == __y && __x.imag() == T(); }

template<typename T>
__host__ __device__ __forceinline__ bool
operator==(const T& __x, const CudaComplex<T>& __y)
{ return __x == __y.real() && T() == __y.imag(); }
//@}

//@{
///  Return false if @a x is equal to @a y.
template<typename T>
__host__ __device__ __forceinline__ bool
operator!=(const CudaComplex<T>& __x, const CudaComplex<T>& __y)
{ return __x.real() != __y.real() || __x.imag() != __y.imag(); }

template<typename T>
__host__ __device__ __forceinline__ bool
operator!=(const CudaComplex<T>& __x, const T& __y)
{ return __x.real() != __y || __x.imag() != T(); }

template<typename T>
__host__ __device__ __forceinline__ bool
operator!=(const T& __x, const CudaComplex<T>& __y)
{ return __x != __y.real() || T() != __y.imag(); }
//@}


} //end namespace SciPAL
#endif // CUDACOMPLEX_H

