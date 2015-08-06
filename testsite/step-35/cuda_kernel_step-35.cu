//@sect3{File: cuda_kernel_step-35.cu}
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

Copyright  Lutz Künneke, Jan Lebert 2014
*/

//Header containing the declarations of the wrapper functions - not the kernels.
//This header is the interface between the code compiled by nvcc and the one compiled by gcc.

//std
#include <stdio.h>

//CUDA
#include <cuComplex.h>

//Thrust
#include <thrust/functional.h>

#ifndef nUSE_CPP11
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#endif


//SciPAL
#include <base/PrecisionTraits.h>
#include <src/cuda/scipal_kernels.cu>
#include <step-35/autoInstantiations.h>

//Our stuff
#include "cuda_helper.h"

#include <step-35/cuda_kernel_wrapper_step-35.cu.h>

//@sect4{Constant CUDA memory}
//Maximum subset size is 512 (reasonable choices are more like 15-30)
#ifdef DOUBLE_PRECISION
__constant__ double cs[512];
#else
__constant__ float cs[512];
#endif

//@sect5{Function: set_cs}
//@brief Copy $c_s$ to constant memory
//@param cs_h pointer to $c_s$ on host
//@param maxnum number of elements in cs_h
template<typename T>
void step35::Kernels<T>::set_cs(T* cs_h, int maxnum) {
    checkCudaErrors(cudaMemcpyToSymbol(cs, cs_h, sizeof(T)*maxnum));
}

//@sect4{Device Functions}
//
//Prior to the kernels we have to define the device functions that we need.

//@sect5{Struct: SharedMemory}
//Workaround as CUDA doesn't immediately support shared memory arrays in templated functions,
//see <a href="http://docs.nvidia.com/cuda/cuda-samples/#simple-templates">CUDA Simple Templates example</a>.
template <typename T>
struct SharedMemory {
    // Ensure that we won't compile any un-specialized types
    __device__ T *getPointer() {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization
template <>
struct SharedMemory <float> {
    __device__ float *getPointer() {
        extern __shared__ float s_float[];
        return s_float;
    }
};

//Double specialization
template <>
struct SharedMemory <double> {
    __device__ double *getPointer() {
        extern __shared__ double s_double[];
        return s_double;
    }
};

//@sect5{Struct: __abs}
//@brief Unary functions to calculate the absolute value for an element
template <typename T>
struct __abs : public thrust::unary_function<T,T> {
    // Ensure that we won't compile any un-specialized types
    __device__
    T operator()(T x) {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization: use fabsf()
template <>
struct __abs <float> : public thrust::unary_function<float,float> {
    __device__
    float operator()(float x) {
        return fabsf(x);
    }
};

//Double specialization: use fabs()
template <>
struct __abs <double> : public thrust::unary_function<double,double> {
    __device__
    double operator()(double x) {
        return fabs(x);
    }
};

//@sect5{Struct: __element_norm_product}
//@brief Unary functions to calculate $\frac{x\cdot y}{n}$ where $x,y$ are complex numbers and $n$ a constant factor.
template <typename T>
struct __element_norm_product {
    __element_norm_product(T _n) {}
    // Ensure that we won't compile any un-specialized types
    __device__
    T operator()(T x) {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization
template <>
struct __element_norm_product <float> {
    cuFloatComplex normalization;

    __element_norm_product(float _n) {
        normalization = make_cuFloatComplex(1.0/_n, 0);
    }

    __host__ __device__
    cuFloatComplex operator()(const cuFloatComplex& x, const cuFloatComplex& y) const {
        return cuCmulf(cuCmulf(x, y), normalization);
    }
};

//Double specialization
template <>
struct __element_norm_product <double> {
    cuDoubleComplex normalization;

    __element_norm_product(double _n) {
        normalization = make_cuDoubleComplex(1.0/_n, 0);
    }

    __host__ __device__
    cuDoubleComplex operator()(const cuDoubleComplex& x, const cuDoubleComplex& y) const {
        return cuCmul(cuCmul(x, y), normalization);
    }
};

//@sect5{Struct: __mult}
//@brief Unary function to calculate elementwise multiplication with a constant value.
template <typename T>
struct __mult : public thrust::unary_function<T,T> {
    const T factor;

    __mult(T _f) : factor(_f) {}

    __device__
    T operator()(T x) {
        return factor*x;
    }
};

//@sect5{Function: power_of_two}
//@param x exponent
//@return $2^x$ for $x \ge 0$. For $x < 0$ returns 0.
__device__ int pow_of_two(int x) {
    if (x < 0) {
        return 0;
    }
    return 1 << x;
}

//@sect5{Function: mysign}
//Sign function
template<typename T>
__device__ T mysign(T in)
{
    if ( in < (T)0)
        return (T)(-1.0);
    else
        return (T)1;
}

// @sect4{Kernel}
//
// CUDA extends C by <i>Kernel</i>s, which are executed in parallel by
// <i>threads</i>. Furthermore, one can define <i>device functions</i>
// which can only be called from inside a kernel.
// <br>
// Note that unlike device functions kernels cannot be members of class.

//@sect5{Kernel: __dykstra}
//@brief This Kernel provides the Gaussian Dykstra projection.
//@param qg pointer to the q vectors in global memory
//@param A path to the image
//@param cinfo array of pointers where to find n,i,j and numpatch
//@param q_offset array of offsets where to find first $q$ value for a threadblock
//@param ni height of the image
//@param nj width of the image
//NOTICE: only 2D functionality
template<typename T>
__global__ void
__dykstra(T *qg, T *A, int **cinfo, const int * q_offset, const int ni, const int nj) {
    //Find out where I live

    int k=0;
    int border=0;
    int *n=cinfo[4*blockIdx.x];
    int *i=cinfo[4*blockIdx.x+1];
    int *j=cinfo[4*blockIdx.x+2];
    int numpatch=cinfo[4*blockIdx.x+3][0];

    while ( k < numpatch ) {
        if ( border+n[k]*n[k] > threadIdx.x )
            break;
        border+=n[k]*n[k];
        k++;
    }
    if ( threadIdx.x - border >= n[k]*n[k] || k >= numpatch) {
        //This is bad, exit
        return;
    }
    //Allocate shared memory
    SharedMemory<T> shared;
    T* s_mem = shared.getPointer();
    T* q=&(s_mem[0]); //takes first n2 sizeof(T) bytes
    T* f=&(s_mem[blockDim.x]); //takes second third
    T* e=&(s_mem[2*blockDim.x]); //takes last third

    int is = (threadIdx.x-border)/n[k] + i[k];
    int js = (threadIdx.x-border)%n[k] + j[k];

    //Pixel index in the current in the global image
    const int idx = is*nj + js;

    //Fill shared memory with variables we later use
    q[threadIdx.x] = qg[threadIdx.x+q_offset[blockIdx.x]];

    f[threadIdx.x] = A[idx] - q[threadIdx.x];
    e[threadIdx.x] = f[threadIdx.x]*f[threadIdx.x];
    __syncthreads();

    //m is the next power of 2 bigger than n divided by two
    int m=n[k]*n[k];
    int mmax=pow_of_two(floor(log2f((float)m)));
    //Sum over all pixels in one frame
    m=1;
    while (m <= mmax) {
        if (threadIdx.x - border + m < n[k]*n[k]) {
            e[threadIdx.x] += e[threadIdx.x + m];
        }
        m = m << 1;
        __syncthreads();
    }

    //$q = x_{r+1} - x_r$
    q[threadIdx.x] = -f[threadIdx.x];
    if (((T) cs[n[k]-1])*e[border] > 1.0) {
        f[threadIdx.x] = f[threadIdx.x]/(sqrt(((T) cs[n[k]-1])*e[border]));
    }

    //Update global memory (image and q)
    A[idx] = f[threadIdx.x];
    qg[threadIdx.x+q_offset[blockIdx.x]] = q[threadIdx.x] + f[threadIdx.x];
}

//@sect5{Kernel: __incomplete_dykstra}
//@brief CUDA adapted approximate Dykstra
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
__global__ void
__incomplete_dykstra(T *A_image, const int ni, const int nj, const int nk, const int offseti, const int offsetj, const int offsetk, const int smin) {
    //Shared memory
    __shared__ float q[6144];
    __shared__ float f[1024];
    __shared__ float e[1024];

    // Set q = 0
    for (int s=5; s>=0; s--) {
        q[s*1024+threadIdx.x]=0;
    }

    // Temporary variables
    int is, js, ks, m;
    // Pixel index this thread is processing
    int idx;

    // Iteration counter
    int it=0;

    // FIXME: make adjustable from outside
    // Maximum iteration count
    const int itmax=300;

    // FIXME: make adjustable from outside
    // Tolerance for convergence check
    const T tol = 1e-2;
    // Increment
    T delta=2.0*tol;
    while ( delta > tol && it < itmax ) {
        it++;
        delta=0;
        // In each threadblock we apply dykstra's algorithm to
        // subsets with edge lengths 32, 16, 8, 4, 2 and 1.
        //\image html s5.png
        //\image html s4.png
        //\image html s3.png
        //\image html s2.png
        //\image html s1.png
        //\image html s0.png
        // Wait for every step before starting the iteration
        __syncthreads();
        for (int s=5; s>=smin; s--) {
            //Edge length of one chunk = $2^s$
            int ChunkLength = pow_of_two(s);
            //Number of chunks in one threadblock = $2^{5-s}$
            int ChunkNum =pow_of_two(5-s);
            //Number of pixels in one = $2^{2\cdot s}$
            int ChunkSize = pow_of_two(2*s);

            //i = blockInd + row * MaxBlock + col * MaxRow * MaxBlock \n
            //row= ( i / MaxBlock ) % MaxRow                          \n
            //col = ( i / MaxRow ) / MaxBlock

            //Line in global image
            is = ChunkLength*((threadIdx.x/ChunkSize)%ChunkNum) +
                 (threadIdx.x%ChunkSize)%ChunkLength + blockIdx.x*32 + offseti;
            //Column in global image
            js = ChunkLength*((threadIdx.x/ChunkSize)/ChunkNum) +
                 (threadIdx.x%ChunkSize)/ChunkLength + blockIdx.y*32 + offsetj;
            ks=offsetk;
            is=is%ni;
            js=js%nj;
            //Pixel index in the current in the global image
            idx = ks*ni*nj + is*nj + js;

#ifdef DY_DEBUG
            if (idx%891==0)
            printf("A_%d : %f, ", idx, A_image[idx] );
#endif
            //Fill shared memory with variables we later use
            f[threadIdx.x] = A_image[idx] - q[s*1024+threadIdx.x];
            e[threadIdx.x] = f[threadIdx.x]*f[threadIdx.x];

            m=1;
            __syncthreads();

            //Sum over all pixels of one chuck, write result to e[(threadIdx.x/ChunkSize)*ChunkSize]
            while (m <= pow_of_two(2*s-1)) {
                if (threadIdx.x -(threadIdx.x/ChunkSize)*ChunkSize + m < ChunkSize) {
                    e[threadIdx.x] += e[threadIdx.x + m];
                }
                m = m << 1; // m = m*2
                __syncthreads();
            }

            //$q = x_{r+1} - x_r$
            q[s*1024+threadIdx.x] = -f[threadIdx.x];

            if (cs[s]*e[(threadIdx.x/ChunkSize)*ChunkSize] > 1.0) {
                f[threadIdx.x] = f[threadIdx.x]/(sqrt(cs[s]*e[(threadIdx.x/ChunkSize)*ChunkSize]));
            }
            // Update $q$
            q[s*1024+threadIdx.x] += f[threadIdx.x];

            // Templatized abs function
            __abs<T> mabs;
            // Calculate increment
            delta+=mabs(A_image[idx]-f[threadIdx.x]);
            // Update image
            A_image[idx] = f[threadIdx.x];
            __syncthreads();
        }
    }
}

//@sect5{Kernel: __sum}
//@brief Elementwise sum with offset in one array. Writes output to arr1.
//@param arr1 first input array, output is written here.
//@param arr2 second input array
//@param offset offset in rows and columns in arr2 to ignore
//@param nx width of the image
//@param ny height of the image
template<typename T>
__global__ void
__sum(T* arr1, T* arr2, int offset, int nx, int ny, int nk) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int ti = i + offset;
    int tj = j + offset;
    int tk=0;
    if ( nk > 1 )
        tk = k + offset;
    else
        tk = 0;

    //Check if offset is larger than image boundaries
    if ( ti >= nx )
        ti=ti-nx;
    if ( tj >= ny )
        tj=tj-ny;
    if ( tk >= nk )
        tk=tk-nk;

    //Do elementwise sum
    if ( k*nx*ny+i*ny+j < nx*ny*nk )
        arr1[k*nx*ny+i*ny+j] = arr1[k*nx*ny+i*ny+j] + arr2[tk*nx*ny+ti*ny+tj];
}

//@sect5{Kernel: __update_lagrangian}
//@brief Kernel for Lagrangian update
template<typename T>
__global__ void
__update_lagrangian(T* lag1, T* lag2, int offset, int nx, int ny, int nz, T alpha1, T alpha2,
                    T* e, T* im, T* m1, T* x, T* z) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int ti = i + offset;
    int tj = j + offset;
    int tk = k + offset;

    //Check if offset is larger than image boundaries
    while ( ti >= nx )
        ti=ti-nx;
    while ( tj >= ny )
        tj=tj-ny;
    while ( tk >= nz )
        tk=tk-nz;
    if ( k*nx*ny+i*ny+j < nx*ny*nz )
    {
        lag1[k*nx*ny+i*ny+j] += alpha1*(im[k*nx*ny+i*ny+j] - e[k*nx*ny+i*ny+j] - m1[tk*nx*ny+ti*ny+tj]);
        lag2[k*nx*ny+i*ny+j] += alpha2*(x[k*nx*ny+i*ny+j] - z[k*nx*ny+i*ny+j]);
    }
}


//@setc5{Kernel: __tv_regularization}
//@brief kernel for tv regularization
template<typename T>
__global__ void
__tv_regularization(T* x,T* z,T* lag,T lambda,T rho,int nx,int ny,int nz)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int maxit=100;
    T sw=0.3;
    int it=0;
    T tol=1e-5;
    T delta=2.0*tol;
    __abs<T> mabs;
    while ( delta > tol && it < maxit )
    {
        delta=0;
        it++;
        T newval;
        newval=-rho*(x[k*nx*ny+i*ny+j]-z[k*nx*ny+i*ny+j])-lag[k*nx*ny+i*ny+j];
        if ( nx > 1 && i < nx-1 )
            newval=newval-lambda*mysign(z[k*nx*ny+(i+1)*ny+j]-z[k*nx*ny+i*ny+j]);
        if ( nx > 1 && i > 0 )
            newval=newval-lambda*mysign(z[k*nx*ny+(i-1)*ny+j]-z[k*nx*ny+i*ny+j]);
        if ( ny > 1 && j < ny-1 )
            newval=newval-lambda*mysign(z[k*nx*ny+i*ny+j+1]-z[k*nx*ny+i*ny+j]);
        if ( ny > 1 && j > 0  )
            newval=newval-lambda*mysign(z[k*nx*ny+i*ny+j-1]-z[k*nx*ny+i*ny+j]);
        if ( nz > 1 && k < nz-1 )
            newval=newval-lambda*mysign(z[(k+1)*nx*ny+i*ny+j]-z[k*nx*ny+i*ny+j]);
        if ( nz > 1 && k > 0 )
            newval=newval-lambda*mysign(z[(k-1)*nx*ny+i*ny+j]-z[k*nx*ny+i*ny+j]);
        delta=delta+mabs(newval);
        __syncthreads();
        z[k*nx*ny+i*ny+j]=z[k*nx*ny+i*ny+j]-sw*newval;
        if ( z[k*nx*ny+i*ny+j] < 0 )
            z[k*nx*ny+i*ny+j]=0;
        __syncthreads();
    }
}

//@sect5{Kernel: __diff}
//@brief Elementwise substraction with offset in one array. Output is written in arr1.
//@param arr1 first input array, output is written here.
//@param arr2 second input array. arr1[i] - arr[j] will be written in arr1.
//@param offset offset in rows and columns in arr2 to ignore
//@param nx width of the image
//@param ny height of the image
//@param nz number of layers in image stack
template<typename T>
__global__ void
__diff(T* arr1, T* arr2, int offset, int nx, int ny, int nz) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int ti = i + offset;
    int tj = j + offset;
    int tk = k + offset;

    // Check if offset is larger than image boundaries
    while ( ti >= nx )
        ti=ti-nx;
    while ( tj >= ny )
        tj=tj-ny;
    while ( tk >= nz )
        tk=tk-nz;
    if ( k*nx*ny+i*ny+j < nx*ny*nz )
    {
        arr1[k*nx*ny+i*ny+j] = arr1[k*nx*ny+i*ny+j] - arr2[tk*nx*ny+ti*ny+tj];
    }
}

//@sect5{Kernel: __prepare_proj}
//@brief Kernel for Projection preparation
template<typename T>
__global__ void
__prepare_proj(T* e, T* im, T* m1, T* lag, T rho, int sigma, int nx, int ny, int nz) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int ti = i + sigma;
    int tj = j + sigma;
    int tk = k + sigma;

    // Check if offset is larger than image boundaries
    while ( ti >= nx )
        ti=ti-nx;
    while ( tj >= ny )
        tj=tj-ny;
    while ( tk >= nz )
        tk=tk-nz;
    if ( k*nx*ny+i*ny+j < nx*ny*nz )
    {

        e[k*nx*ny+i*ny+j] = im[k*nx*ny+i*ny+j]-m1[tk*nx*ny+ti*ny+tj];//+lag[k*nx*ny+i*ny+j]/rho;

#ifdef DY_DEBUG
            if ((k*nx*ny+i*ny+j)%8891==0)
                 printf("e_%d : %f, im_%d : %f, m1_%d : %f |", k*nx*ny+i*ny+j,  e[k*nx*ny+i*ny+j] , k*nx*ny+i*ny+j, im[k*nx*ny+i*ny+j], tk*nx*ny+ti*ny+tj, m1[tk*nx*ny+ti*ny+tj] );
#endif

    }
}

//@sect5{Kernel: __soft_thresholding}
//@brief Applies a soft threshold to the 1 Norm of a real vector
//@param arr target array
//@param lag lagrangian for soft threholding in ADMM
//@param xvec source array
//@param rho weight of quadratic penalty term in ADMM
//@param gamma penalty for non-vanishing elements in 1 Norm
//@param size number of elements in array
template<typename T>
__global__ void
__soft_thresholding(T* arr,T* lag,T* xvec,T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx] > 0 ) {
            arr[idx]=xvec[idx] +(lag[idx] - gamma)/rho;///rho)/(rho+lagr[idx]);
            if (arr[idx] < 0 )
                arr[idx]=0;
        } else {
            arr[idx]=xvec[idx] +(lag[idx] + gamma)/rho;///rho)/(rho+lagr[idx]);
            if (arr[idx] > 0 )
                arr[idx]=0;
        }
    }
}

//@sect5{Kernel: __soft_thresholding_complex}
//@brief Kernel for comlex soft thresholding to the 0 Norm of a complex vector
template<typename T>
__global__ void
__soft_thresholding_complex(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr,
                            T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx].x > 0 ) {
            arr[idx].x=arr[idx].x - gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].x < 0 )
                arr[idx].x=0;
        } else {
            arr[idx].x=arr[idx].x + gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].x > 0 )
                arr[idx].x=0;
        }
        if ( arr[idx].y > 0 ) {
            arr[idx].y=arr[idx].y - gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].y < 0 )
                arr[idx].y=0;
        } else {
            arr[idx].y=arr[idx].y + gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].y > 0 )
                arr[idx].y=0;
        }
    }
}


//@sect5{Kernel: __sparse_thresholding}
//@brief Kernel for real valued soft thresholding to enforce 0 Norm sparsity in ADMM
//@param arr target array
//@lag Lagrangian
//@param rho quadratic penalty term for ADMM
//@param gamma penalty for non-vanishing entries in 0-Norm
//@parma size size of the array
template<typename T>
__global__ void
__sparse_thresholding(T* arr, T *lag, T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx] > 0 ) {
            arr[idx]=arr[idx]-(gamma/rho)/(rho+lag[idx]);
            if ( arr[idx] < 0 )
                arr[idx]=0;
        } else {
            arr[idx]=arr[idx]+(gamma/rho)/(rho+lag[idx]);
            if ( arr[idx] > 0 )
                arr[idx]=0;
        }
    }
}

//@sect5{Kernel: __pseudo_inverse}
//@brief Kernel to minimize the 2-Norm of arr with a constraint. The solution corresponds to the pseudo-inverse solution of the constraint
template<typename T>
__global__ void
__pseudo_inverse(T* arr, T *lag, T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        arr[idx]=(rho*arr[idx]+lag[idx])/(rho+gamma);
    }
}
//@sect5{Haar Wavelet Transformation}
//@sect6{Kernel: __haar_horizontal}
//@brief Calculates 1D Haar Wavelet Transformation on all rows for one step with step width $w$.
//@see <a href="http://unix4lyfe.org/haar/">Here</a> and <a href="http://www.tilman.de/uni/ws05/scivis/2d-transformation.html">here</a> for a complete Haar Wavelet Transformation algorithm details.
//@param in input
//@param out output
//@param ny height of the image
//@param w step width for one step

//Treat the first $w$ elements of each row as $w/2$ pairs called $(a, b)$. \n
//Calculate $(a + b) / \sqrt{2}$ for each pair, these values will be the first half of the considered of elements in the output array. \n
//Calculate $(a - b) / \sqrt{2}$ for each pair, these values will be the second half of the considered of elements in the output array.\n
//For $w = 8$ and input array \n
//\htmlonly <pre>[  1.000   2.000   3.000   1.000   2.000   3.000   4.000   0.000  ]</pre>\endhtmlonly
//will result into \htmlonly <pre>[  2.121   2.828   3.536   2.828  -0.707   1.414  -0.707   2.828  ]</pre>\endhtmlonly
//Treating this output array as input array with $w=4$ will yield \n
//\htmlonly <pre>[  3.500   4.500  -0.500   0.500  -0.707   1.414  -0.707   2.828  ]</pre>\endhtmlonly
template<typename T>
__global__ void __haar_horizontal(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (x_idx >= (w+1)/2 || y_idx >= w)
        return;

    int in_idx  = 2*x_idx + y_idx*ny;
    int sum_idx =   x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    out[sum_idx]       = (in[in_idx] + in[in_idx+1])*one_over_sqrt_2;
    out[sum_idx + w/2] = (in[in_idx] - in[in_idx+1])*one_over_sqrt_2;
}

//@sect6{Kernel: __inverse_haar_horizontal}
//@brief Inverse 1D row Haar Wavelet Transformation step. Inverts __haar_horizontal
template<typename T>
__global__ void __inverse_haar_horizontal(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (x_idx >= (w+1)/2 || y_idx >= w)
        return;

    int sum_idx      =   x_idx + y_idx*ny;
    int original_idx = 2*x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    // $a = ((a+b)/\sqrt{2}+(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx]   = (in[sum_idx] + in[sum_idx + w/2])*one_over_sqrt_2;
    // $b = ((a+b)/\sqrt{2}-(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx+1] = (in[sum_idx] - in[sum_idx + w/2])*one_over_sqrt_2;
}


//@sect6{Kernel: __haar_vertical}
//@brief Calculates 1D Haar Transformation on all columns for one step with width $w$.
//@see __haar_horizontal
template<typename T>
__global__ void __haar_vertical(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (y_idx >= (w+1)/2 || x_idx >= w)
        return;

    int in_idx1 = x_idx + 2*y_idx*ny;
    int in_idx2 = x_idx + (2*y_idx+1)*ny;
    int sum_idx = x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    out[sum_idx]        = (in[in_idx1]+in[in_idx2])*one_over_sqrt_2;
    out[sum_idx+ny*w/2] = (in[in_idx1]-in[in_idx2])*one_over_sqrt_2;
}

//@sect6{Kernel: __inverse_haar_vertical}
//@brief Inverse 1D column Haar Wavelet Transformation step. Inverts __haar_vertical
template<typename T>
__global__ void __inverse_haar_vertical(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (y_idx >= (w+1)/2 || x_idx >= w)
        return;

    int original_idx1 = x_idx + 2*y_idx*ny;
    int original_idx2 = x_idx + (2*y_idx+1)*ny;
    int sum_idx       = x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    // $a = ((a+b)/\sqrt{2}+(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx1] = (in[sum_idx]+in[sum_idx+ny*w/2])*one_over_sqrt_2;
    // $b = ((a+b)/\sqrt{2}-(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx2] = (in[sum_idx]-in[sum_idx+ny*w/2])*one_over_sqrt_2;
}

//@sect4{Kernel wrapper functions}
//The functions which call the kernels

//@sect5{Haar Wavelet Transformation}

//@sect6{Function: haar}
//@brief Calculates full 2D nonstandard Haar Wavelet decomposition.
//@param in_d pointer to the array on device, will be transformed inplace. Has to be of quadratic size.
//@param tmp_d pointer to an array on device with same size as in_d. Will be filled with garbage.
//@param ny edge length of the matrix, has to be of the power of two

//Alternating row and column 1D Haar Wavelet steps
template<typename T>
void step35::Kernels<T>::haar(T* in_d, T* tmp_d, const int ny) {
    dim3 blocksize(32,32);
    dim3 gridsize;

    int w = ny;
    gridsize.x=(w+blocksize.x-1)/blocksize.x;
    gridsize.y=(w+blocksize.y-1)/blocksize.y;

    while(w>1) {
        __haar_horizontal<T><<<gridsize,blocksize>>>(in_d, tmp_d, ny, w);
        __haar_vertical<T><<<gridsize,blocksize>>>(tmp_d, in_d, ny, w);
        w /= 2;
    }
    cudaDeviceSynchronize();
    getLastCudaError("haar() execution failed\n");
}

//@sect6{Function: inverse_haar}
//@brief Inverse 2D nonstandard Haar Wavelet decomposition algorithm
//@param in_d pointer to the array on device, will be transformed inplace. Has to be of quadratic size.
//@param tmp_d pointer to an array on device with same size as in_d. Will be filled with garbage.
//@param ny edge length of the matrix, has to be of the power of two
//@see haar()
template<typename T>
void step35::Kernels<T>::inverse_haar(T* in_d, T* tmp_d, const int ny) {
    dim3 blocksize(32,32);
    dim3 gridsize;

    int w = 2;
    gridsize.x=(ny+blocksize.x-1)/blocksize.x;
    gridsize.y=(ny+blocksize.y-1)/blocksize.y;

    while(w <= ny) {
        __inverse_haar_horizontal<T><<<gridsize,blocksize>>>(in_d, tmp_d, ny, w);
        __inverse_haar_vertical<T><<<gridsize,blocksize>>>(tmp_d, in_d, ny, w);
        w *= 2;
    }
    cudaDeviceSynchronize();
    getLastCudaError("inverse_haar() execution failed\n");
}

//@sect5{Function: dykstra}
//@brief This is a wrapper for the __dykstra Kernel function.
//@param qg pointer to the q vectors on device
//@param A pointer to the image on device
//@param cinfo array of device pointers where to find n,i,j and numpatch
//@param q_offset array of offsets where to find first $q$ value for each cluster
//@param num_of_cluster numbers of clusters (threadblocks)
//@param ni height of the image
//@param nj width of the image
//@param maximum number of pixels per threadblock
//@param mystream a CUDA stream
template<typename T>
void step35::Kernels<T>::dykstra(T *qg, T *A,  int **cinfo, const int * q_offset, const int num_of_cluster,
                                 const int ni, const int nj, const int numpx, cudaStream_t *mystream) {
    dim3 grid(num_of_cluster,1);
    dim3 blocks(numpx,1);
    __dykstra<T><<<grid, blocks, 3*numpx*sizeof(T), *mystream>>> (qg, A, cinfo, q_offset, ni, nj);
}

//@sect5{Function: dyadyc_dykstra}
//@brief This is a wrapper for the __dyadyc_dykstra Kernel function.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
void step35::Kernels<T>::dyadic_dykstra(T *A_image, const int ni, const int nj, const int nk, const int offseti, const int offsetj, const int offsetk, const int so) {
    int gridsi=ni/32;
    int gridsj=nj/32;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(1024,1);
    __incomplete_dykstra<T><<<grid,blocks>>> (A_image, ni, nj, nk, offseti, offsetj, offsetk, so);
    getLastCudaError("__incomplete_dykstra<<<>>> execution failed\n");
    cudaDeviceSynchronize();
}

//@sect5{Function: sum}
//@brief Elementwise sum with offset in one array. Inplace output to arr1.
//@param arr1 first input array, output is written here.
//@param arr2 second input array
//@param offset offset in rows and columns in arr2 to ignore
//@param ni width of the image
//@param nj height of the image
template<typename T>
void step35::Kernels<T>::sum(T* arr1, T* arr2, int offset, int ni, int nj, int nk) {
    int gridsi=ni/32;
    int gridsj=nj/32;
    int gridsk=nk;//in 2D case defaults to 1
    dim3 grid(gridsi,gridsj,gridsk);
    dim3 blocks(32,32,1);
    __sum<T><<<grid,blocks>>>(arr1, arr2, offset, ni, nj, nk);
    getLastCudaError("__sum<<<>>> execution failed\n");
    cudaDeviceSynchronize();
}

//@sect5{Function: diff}
//@brief Elementwise substraction with offset in one array. Output is written in arr1.
//@param arr1 first input array, output is written here.
//@param arr2 second input array. arr1[i] - arr[j] will be written in arr1.
//@param offset offset in rows and columns in arr2 to ignore
//@param width width of the array
//@param height height of the array
template<typename T>
void step35::Kernels<T>::diff(T* arr1, T* arr2, int offset, int width, int height, int depth) {
    int gridsi=width/32;
    int gridsj=height/32;
    int gridsk=depth;
    dim3 grid(gridsi,gridsj,gridsk);
    dim3 blocks(32,32,1);
    __diff<T><<<grid,blocks>>>(arr1, arr2, offset, width, height, depth);
    getLastCudaError("__diff<<<>>> execution failed\n");
    cudaDeviceSynchronize();
}

// @sect5{Function: element_norm_product}
//@brief Elementwise multiplication of two arrays with complex values and division by the number of elements in each array.
//@param arr1 first input array, output will be calculated inplace here
//@param arr2 second input array
//@param width width of the array
//@param height height of the array
template<typename T>
void step35::Kernels<T>::element_norm_product(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr1,
        typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr2,
        int width, int height, int depth) {
    int size=depth*width*(height/2+1);
    T normalization=width*height*depth;
    __element_norm_product<T> op(normalization);
#ifndef nUSE_CPP11
    thrust::device_ptr<typename PrecisionTraits<T, gpu_cuda>::ComplexType> ptr1 = thrust::device_pointer_cast(arr1);
    thrust::device_ptr<typename PrecisionTraits<T, gpu_cuda>::ComplexType> ptr2 = thrust::device_pointer_cast(arr2);
    thrust::transform(thrust::device, ptr1, ptr1 + size, ptr2, ptr1, op);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("element_norm_product<<<>>> execution failed\n");
}

//@sect5{Function: mult}
//@brief Inplace elementwise multiplication with a constant value.
//@param arr input array
//@param factor factor for the elementwise multiplication
//@param size number of elements in the input array
template<typename T>
void step35::Kernels<T>::mult(T* arr, T factor, int size) {
    __mult<T> op(factor);
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::transform(thrust::device, ptr, ptr + size, ptr, op);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::transform __mult execution failed\n");
}

//@sect5{Function: reset}
//@brief Fills input array with zeros
//@param arr input array
//@param size number of elements in the input array
template<typename T>
void step35::Kernels<T>::reset(T* arr, int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::fill(thrust::device, ptr, ptr + size, 0);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::fill execution failed\n");
}

//@sect5{Function: abs}
//@brief Replaces every element in an array with it the absolute value of the element
//@param arr input array
//@param size number of elements in array
template<typename T>
void step35::Kernels<T>::abs(T* arr, int size) {
    __abs<T> op;
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::transform(thrust::device, ptr, ptr + size, ptr, op);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::transform __abs execution failed\n");
}

//@sect5{Function: update_lagrangian}
//@brief Lagrangian Multipliers Update on device
//@param alpha1 Update step for first constraint
//@param alpha2 Update step for second constraint
//@param lag1 Lagrangian for first constraint
//@param lag2 Lagrangian for second constraint
//@param sigma half width of psf
//@param nx width of image
//@param ny height of image
//@param nz number of images in case of 3d
//@param im pointer to noisy image
//@param m1 estimate convoluted with point spread function
//@param x current estimate
//@param z smoothed estimate
template<typename T>
void step35::Kernels<T>::update_lagrangian(T* lag1, T* lag2, int sigma, int nx, int ny, int nz,
                                           T alpha1, T alpha2, T* e, T* im, T* m1, T* x, T* z) {
    int gridsi=nx/32;
    int gridsj=ny/32;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(32,32);
    __update_lagrangian<T><<<grid,blocks>>>(lag1, lag2, sigma,nx,ny,nz,alpha1,alpha2,e,im,m1,x,z);
    cudaDeviceSynchronize();
    getLastCudaError("__update_lagrangian<<<>>> execution failed\n");
}

//@sect5{Function: prepare_e}
//@brief Put the vector to be projected to *e before Dykstra
//@param rho rho1 from first constraint
//@param e Field of Residuals
//@param im Noisy image
//@param m1 estimate convoluted with psf
//@param lag Lagrangian
template<typename T>
void step35::Kernels<T>::prepare_e(T* e, T* im, T* m1, T* lag, T rho, int sigma, int nx, int ny, int nz) {
    int gridsi=nx/32;
    int gridsj=ny/32;
    int gridsk=nz;
    dim3 grid(gridsi,gridsj,gridsk);
    dim3 blocks(32,32,1);
    __prepare_proj<T><<<grid,blocks>>>(e,im,m1,lag,rho,sigma,nx,ny,nz);
    cudaDeviceSynchronize();
    getLastCudaError("__prepare_e<<<>>> execution failed\n");
}

//@sect5{Function: tv_regularization}
//@brief TODO
template<typename T>
void step35::Kernels<T>::tv_regularization(T* x, T* z, T* lag, T lambda, T rho, int nx, int ny, int nz) {
    int gridsi=nx/32;
    int gridsj=ny/32;
    int gridsk=nz;
    dim3 grid(gridsi,gridsj,gridsk);
    dim3 blocks(32,32,1);
    __tv_regularization<T><<<grid,blocks>>>(x,z,lag,lambda,rho,nx,ny,nz);
    cudaDeviceSynchronize();
    getLastCudaError("__tv_regularization<<<>>> execution failed\n");
}

//@sect5{Function: soft_threshold}
//@brief Applies a soft threshold. Absolute values smaller than a certain threshold will be to zero,
//       larger values will be shifted by the value of the threshold closer to zero.
//@param arr device array to sort, result will be inplace
//@param threshold threshold value
//@param size number of elements in array
template<typename T>
void step35::Kernels<T>::soft_threshold(T* arr, T* lag,T *xvec, T rho, T threshold, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __soft_thresholding<T><<<grid,blocks>>>(arr, lag,xvec,rho, threshold, size);
    cudaDeviceSynchronize();
    getLastCudaError("__soft_threshold<<<>>> execution failed\n");
}

//@sect5{Function: soft_threshold_complex}
//@brief Soft thresholding for complex fields
//@param threshold threshold for both complex and real part
//@param arr inplace thresholding on arr
//@param size size of arr
template<typename T>
void step35::Kernels<T>::soft_threshold_complex(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr,
        T threshold, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __soft_thresholding_complex<T><<<grid,blocks>>>(arr,threshold,size);
    cudaDeviceSynchronize();
    getLastCudaError("__soft_thresholding_complex<<<>>> execution failed\n");
}


//@sect5{Function: sparse}
//@brief Corresponds to soft thresholding in real fields
//@param gamma thresholding parameter
//@param arr array to perform inplace thresholding on
//@param lag Lagrangian
//@param rho parameter rho from ADMM
//@param size size of arr
template<typename T>
void step35::Kernels<T>::sparse(T* arr, T* lag, T rho, T gamma, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __sparse_thresholding<T><<<grid,blocks>>>(arr,lag,rho,gamma,size);
    cudaDeviceSynchronize();
    getLastCudaError("sparse execution failed\n");
}
template<typename T>

//@sect5{Function: pseudo_inverse}
//@brief Pseudoinverse
void step35::Kernels<T>::pseudo_inverse(T* arr, T* lag, T rho, T gamma, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __pseudo_inverse<T><<<grid,blocks>>>(arr,lag,rho,gamma,size);
    cudaDeviceSynchronize();
    getLastCudaError("__pseudo_inverse execution failed\n");
}

//@sect5{Function: sort}
//@brief Sorts array using thrust::sort
//@param arr device array to sort, result will be inplace
//@param size number of elements in array
//@brief Wrapper around thrust::sort
template<typename T>
void step35::Kernels<T>::sort(T* arr,int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::sort(thrust::device, ptr, ptr + size);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::sort execution failed\n");
}

//@sect5{Function: sum_all}
//@return Sum over all elements
//@param arr device array
//@param size number of elements in array
//@brief Wrapper around thrust::reduce
template<typename T>
T step35::Kernels<T>::sum_all(T* arr,int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    T sum = thrust::reduce(thrust::device, ptr, ptr + size, 0, thrust::plus<T>());
#else
    T sum = -1;
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::reduce plus execution failed\n");
    return sum;
}

// Finally, we have to specialize the templates to force the compiler to actually compile something.
// This has to be at the end of file, because all functions have to be declared and their bodies defined before the
// class can be explictly instantiated by the compiler.
template class step35::Kernels<float>;
template class step35::Kernels<double>;
