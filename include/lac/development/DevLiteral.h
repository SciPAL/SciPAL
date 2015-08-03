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


#ifndef DEVLITERAL_H
#define DEVLITERAL_H

namespace SciPAL {


//! @sect3{Struct: DevLiteral}
//!
//! This class encapsulates simple numbers so that they can be used in expression tree leafs.
//! We put into the same file as @p Expr because the only reason to introduce a class for
//! literals is that the template magic in the ET business does not understand primitive
//! built-in types as they are not classes.
template <typename T>
class DevLiteral
{
public:
    typedef T value_type;

    typedef DevLiteral<T> Type;

    typedef const Type ConstHandle;

    static const bool is_leaf = true;


    __host__ __device__ __forceinline__
    DevLiteral(const T _val) : val(_val) {
         // printf("%s ", __PRETTY_FUNCTION__); printf("literal src : %g : dst : %g\n", _val, val );
//        cudaMemcpyToSymbol(&val, &_val, sizeof(T) );

    }


    __host__ __device__ __forceinline__
    operator const T () const {
        return val;
    }


    __host__ __device__ __forceinline__
    T get_val() const {
      //  printf("literal get_val() : %g\n", val );
        return val;
    }

private:
    //! This attribute stores the actual value of the scalar coefficient used in the original expression.
    const T val;
};

}

#endif // DEVLITERAL_H
