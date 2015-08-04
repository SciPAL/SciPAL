#ifndef CUDA_UTILS_CU_H
#define CUDA_UTILS_CU_H


// TMPFIX

#define uint unsigned int
#define NUM_THREADS 256
#define BLOCK_DIM_Y 16

bool is_power_of_two (uint n);

__device__
__forceinline__
uint ColMajor_index_2D(uint row, uint col, uint col_length);

template<int dim>
__forceinline__
__device__ double scalardiff( const double *normal, const uint normal_size,
                              const double *x, const uint x_size,
                              const double *q, const uint q_size);

template<int dim>
__forceinline__
__device__ double single_layer(const double *x, const uint x_size,
                              const double *q, const uint q_size);

__device__
__forceinline__
double3 cross(const double3 &a, const double3 &b);

__device__
__forceinline__
double dot(const double3 &a, const double3 &b);

__device__ double atomicAdd(double* address, double val);

#endif //CUDA_UTILS_CU_H

