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


#ifndef LITERAL_H
#define LITERAL_H
#include <lac/Expr.h>
#include<lac/BlasVectorOperations.h>

namespace SciPAL {

//forward declarations
template <typename T> class DevLiteral;
template <typename T> struct Setter;

// @sect3{Struct: Literal}
//
// This class encapsulates simple numbers so that they can be used in expression tree leafs.
// We put into the same file as @p Expr because the only reason to introduce a class for literals
// is that tempalte magic in the ET business does nto understand primitive built-n types as they are not classes.
template <typename T>
class Literal : public Expr<Literal<T> > {
public:
    typedef T value_type;

    typedef Literal<T> Type;
    typedef DevLiteral<T> DevType;
    typedef const Type ConstHandle;

    __host__ __device__
    Literal(const T _val) : val(_val) {}

    template<typename T2>
    __host__ __device__
    Literal(const T2 _val) : val(T(_val)) {}

    operator value_type () const {
    return val;
    }

     __host__ __device__ __forceinline__
     T get_val() const { return val; }

     Literal<T> operator =(T new_val) {
         val = new_val;
         return *this;
     }

     Literal<T> operator =(Literal<T> new_val) {
         val = new_val.get_val();
         return *this;
     }

     template <typename X>
     SciPAL::Literal<T> & operator = (const SciPAL::Expr<X> & e);


private:
     T val;
};
template <typename T>
template <typename X>
Literal<T> & Literal<T>::operator =
(const SciPAL::Expr<X> & e)
{
#ifdef DEBUG
    std::cout << "line :" << __LINE__ << ", Vector<T,BW>  " << __FUNCTION__<< "\n"  << std::endl;
#endif

    SciPAL::LAOOperations::apply<T, cublas>(*this, ~e);
    return *this;
}

}
#endif // LITERAL_H
