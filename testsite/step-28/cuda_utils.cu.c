#ifndef CUDA_UTILS_CU_C
#define CUDA_UTILS_CU_C

#include "cuda_utils.cu.h"

// @sect5{Function: is_power_of_two}
// Returns whether @p n is a power of two or not.
// This is needed to run the correct kernel. It basically counts the number
// of 1-bits in the unsigned integer @p n. If there is only one bit set to '1'
// (and @p n is not zero), @p n is a power of two.

bool is_power_of_two (uint n)
{
  return (((n & (~n + 1)) == n) && (n != 0));
}


// @sect4{Device Functions}
//
// Prior to the kernels we have to define the device functions that we need.

// @sect5{Device Function: CStyle_index_2D}

// This function maps a row index $row$ and column index $col$
// of a matrix entry to the position $row\cdot row\_length + col$ in a linear array,
// which stores a matrix of row length
// $row\_length$ in row-major (C-style) order.
// @param row : row index
// @param col : column index
// @param row_length : row length

__device__
__forceinline__
uint CStyle_index_2D(uint row, uint col, uint row_length)
{
    return col + row*row_length;
}

// @sect5{Device Function: ColMajor_index_2D}

// This function maps a row index $row$ and column index $col$
// of a matrix entry to the position $row + col\cdot col_length$ in a linear array,
// which stores a matrix of column length
// $col\_length$ in column-major order.
// @param row : row index
// @param col : column index
// @param col_length : column length

__device__
__forceinline__
uint ColMajor_index_2D(uint row, uint col, uint col_length)
{
    return row + col*col_length;
}




// @sect5{Device Function: scalardiff}

// This function takes three points $n$, $x$ and $q$ and returns $n\cdot(x-q)$
// Because of the way these points lie in memory, the sizes of the vectors containing
// the points are needed as well:
// e.g. @p normal has the components: @p normal, (@p normal + @p normal_size),
// (@p normal + 2 * @p normal_size)
// @param normal : pointer to the first component of the normal vector
// @param normal_size : number of points in the normal vector
// @param x: pointer to the first component of the x vector
// @param x_size : number of points in the x vector
// @param q: pointer to the first component of the q vector
// @param q_size : number of points in the q vector

template<int dim>
__forceinline__
__device__ double scalardiff( const double *normal, const uint normal_size,
                              const double *x, const uint x_size,
                              const double *q, const uint q_size)
{
    double val = 0;
    for(uint i = 0; i < dim; ++i){
        val += normal[i*normal_size]*(x[i*x_size]-q[i*q_size]);
    }
    return val;
}


// @sect5{Device Function: single_layer}

// Returns the value of the Green's function $1/|x-q|$, needed to calculate
// the value of the single layer operator
// @param x: pointer to the first component of the x vector
// @param x_size : number of points in the x vector
// @param q: pointer to the first component of the q vector
// @param q_size : number of points in the q vector

template<int dim>
__forceinline__
__device__ double single_layer(const double *x, const uint x_size,
                              const double *q, const uint q_size)
{
    double val = 0;
    for(uint i = 0; i < dim; ++i){
        val += pow((x[i*x_size]-q[i*q_size]),2);
    }
    // @p rsqrt is the reverse square root provided by CUDA
    return rsqrt(val);
}

__device__
__forceinline__
double3 cross(const double3 &a, const double3 &b)
{
    return make_double3(a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x);
}

__device__
__forceinline__
double dot(const double3 &a, const double3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
    //return 1.;
}


__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
            (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}


#endif
