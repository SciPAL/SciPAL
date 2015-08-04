// @sect3{File: cuda_kernel_step-28.cu.c}

// Include the header containing the declarations of the kernel wrapper functions - not the kernels themselves.
// This header is the interface between the code compiled by nvcc and the one compiled by gcc:
#include <step-28/cuda_kernel_wrapper_step-28.cu.h>
// For @p printf() debugging in CUDA
#include <stdio.h>


#include <step-28/cuda_utils.cu.h>
#include <lac/release/ShapeData.h>

// @sect5{Device Function: SL_to_DL}

// This function takes the calculated value of the single layer operator
// (and some other parameters) and returns the value of the double layer
// operator. It uses takes the equation for the double-layer operator
// given in the Introduction (DL):
// \f{align*}
//     K_{h,ij}&=\frac{1}{2}\delta_{ij}-\sum_{E\subset S_j}\sum_{a\in E}H_{ia}
//     w_{aj}JxW_a
// \f}
// Ignoring the $\frac{1}{2}\delta_{ij}$ for now:
// \f{align*}
// H_{ia}w_{aj}JxW_a &= \vec{n}_a \frac{\vec{x}_i-\vec{x}_a^\prime}
// {|\vec{x}_i-\vec{x}_a^\prime|^3}w_{aj}JxW_a = \vec{n}_a(\vec{x}_i-\vec{x}_a^\prime)
// G_{ia}^3\cdot w_{aj}JxW_a \\
// &= \vec{n}_a(\vec{x}_i-\vec{x}_a^\prime)G_{ia}^2\cdot V_{h,ij}
// \f}
// With the definitions given in the introduction.
// @param normal : pointer to the first component of the normal vector, $\mathrel{\widehat{=}} \vec{n}_a$
// @param normal_size : number of points in the normal vector
// @param x: pointer to the first component of the x vector, $\mathrel{\widehat{=}} \vec{x}_i$
// @param x_size : number of points in the x vector
// @param q: pointer to the first component of the q vector, $\mathrel{\widehat{=}} \vec{x}_a^\prime$
// @param q_size : number of points in the q vector
// @param G_ia : value of the Green's function
// @param SL_ij : value of the single layer operator, $\mathrel{\widehat{=}} V_{h,ij}$

template<int dim>
__device__ double SL_to_DL(const double *normal, const uint normal_size,
                           const double *x, const uint x_size,
                           const double *q, const uint q_size,
                           const double G_ia, const double SL_ij)
{
    return scalardiff<dim>(normal, normal_size, x, x_size, q, q_size)*pow(G_ia,2)*SL_ij;
}


// @sect4{Kernel}
//
// The kernels used in this project to assemble the matrices


// @sect5{Kernel: __assemble_bem_matrices_Po2}
//
// Kernel to assemble the BEM-Matrices if the number of quadrature points
// is a power of two. This greatly simplifies the calculation
template<int dim>
__global__
void
__assemble_bem_matrices_Po2(SciPAL::ShapeData<double> DL_matrix,
                            SciPAL::ShapeData<double> SL_matrix,
                            const SciPAL::ShapeData<double> x, // support points
                            const SciPAL::ShapeData<double> q, // q points
                            const SciPAL::ShapeData<double> n, // normals
                            const SciPAL::ShapeData<double> W  // matrix with trial function values * JxW
                            // old list of args:
                            /* double * DL_matrix, double * SL_matrix,
                                                    const int DL_SL_n_rows,
                                                    const int DL_SL_n_cols,
                                                    const double * x, const int x_size,
                                                    const double * __restrict q, const int q_size,
                                                    const double * __restrict n, const int n_size,
                                                    const double * __restrict W, const int W_n_rows, const int W_n_cols*/
                            )
{
    // Each thread gets one pair (i,a) to calculate.
    const uint a = threadIdx.x;
    const uint i = blockIdx.y * blockDim.y + threadIdx.y;

    // The index of the location in shard memory where we store
    // this thread's contributino to the sums.
    const uint tidy_a_bDy = ColMajor_index_2D(threadIdx.y, a, blockDim.y);

    const uint n_support_points = x.n_rows/dim; // x.n_rows has to become n_...points

    const uint n_q_points = q.n_rows/dim;

    const uint n_normals = n.n_rows/dim;

    const uint num_js = W.n_cols;

    // Next set pointers to the data that does not change during the execution of an instance of
    // the kernel.

    const double * x_i =  (x.data_ptr + i);

    const double * q_a = (q.data_ptr + a);

     const double * n_a = (n.data_ptr + a);

    __shared__ double SL_sum [NUM_THREADS]; // Each thread gets a double val
    __shared__ double DL_sum [NUM_THREADS];

    // Range check
    if( i > n_support_points-1) return;

    // Calculate G_ia
    double G_ia = single_layer<dim>(x_i, n_support_points,
                                    q_a, n_q_points);

    // For loop over j
    for(uint j = 0; j < num_js; ++j){
        // Each thread calculates 'its' single layer value
        double local_SL = G_ia * W.data_ptr[ColMajor_index_2D(a, j, W.n_rows)];

        // and writes it to shared memory
        SL_sum[tidy_a_bDy] = local_SL;

        // Same with the double layer
        DL_sum[tidy_a_bDy] =
                SL_to_DL<dim>( n_a, n_normals,
                               x_i, n_support_points,
                               q_a, n_q_points,
                               G_ia, local_SL);
        // Because we wrote to shared memory, we need to synchronize
        __syncthreads();
        // Now reduce this in each row
        // This is where the power of two is important:
        // This naive reduction only works if blockDim.x is a power of two
        for(int k = blockDim.x/2; k > 0; k /=2)
            if(a < k)
            {
                uint ps_pos = ColMajor_index_2D(threadIdx.y, a + k, blockDim.y);

                SL_sum[tidy_a_bDy] += SL_sum[ps_pos];
                DL_sum[tidy_a_bDy] += DL_sum[ps_pos];
            }

        // Again, because we wrote to shared memory, we need to synchronize
        __syncthreads();
        if(a == 0)
        { // Write to global Matrix
            SL_matrix.data_ptr[ColMajor_index_2D(i, j, SL_matrix.n_rows)] = -SL_sum[tidy_a_bDy];

            DL_matrix.data_ptr[ColMajor_index_2D(i, j, DL_matrix.n_rows)] = +DL_sum[tidy_a_bDy];
        }
    }
}

// @sect5{Kernel: __assemble_bem_matrices}
//
// Kernel to assemble the BEM-Matrices if the number of quadrature points
// is NOT a power of two, but only divisible by four.
template<int dim>
__global__
void
__assemble_bem_matrices(double * DL_matrix, double * SL_matrix,
                        const int DL_SL_n_rows,
                        const int DL_SL_n_cols,
                        const double * x, const int x_size,
                        const double * __restrict q, const int q_size,
                        const double * __restrict n, const int n_size,
                        const double * __restrict W, const int W_n_rows, const int W_n_cols
                        )
{
    // Each thread gets one i. a's are calculated in blocks of 4
    uint i = blockDim.x * blockIdx.x + threadIdx.x;
    uint local_i = threadIdx.x;
    uint thread_a_offset = threadIdx.y;

    uint num_js = W_n_cols;
    uint num_as = q_size;
    uint num_is = x_size;


    // Shared memory for the double layer and the single layer matrices
    __shared__ double SL_sum [NUM_THREADS];
    __shared__ double DL_sum [NUM_THREADS];

    // Range check
    if(i >= num_is )
        return;



    // First, precalculate the needed values of G_ia in each thread
    double G_i [NUM_THREADS/4];
    uint ctr = 0;
    for( uint a = 0; a < num_as; a += BLOCK_DIM_Y)
        G_i[ctr++] = single_layer<dim>( (x + i), x_size,
                                        (q + a + thread_a_offset), q_size);

    double local_SL;
    // Loop over all j's
    for( uint j = 0; j < num_js; ++j){
        // Calculate a's in snippets of BLOCK_DIM_Y (=4)
        ctr = 0;
        // Initialize shared mem to zero:
        SL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)] = 0;
        DL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)] = 0;

        for( uint a = 0; a < num_as; a += BLOCK_DIM_Y){

            local_SL =  G_i[ctr] * W[ColMajor_index_2D(a + thread_a_offset, j, W_n_rows)];

            SL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)]
                    += local_SL;

            DL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)]
                    += SL_to_DL<dim>( (n+a+thread_a_offset), n_size, (x+i), x_size,
                                      (q+a+thread_a_offset), q_size,
                                      G_i[ctr], local_SL);
            ++ctr;
            __syncthreads();
        }

        // now reduce the matrices of partial sums in each row
        for(int k = blockDim.y/2; k > 0; k /=2){
            if(thread_a_offset < k){
                SL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)]
                        += SL_sum[ColMajor_index_2D(local_i,
                                                    thread_a_offset + k, blockDim.x)];

                DL_sum[ColMajor_index_2D(local_i, thread_a_offset, blockDim.x)]
                        += DL_sum[ColMajor_index_2D(local_i,
                                                    thread_a_offset + k, blockDim.x)];
            }
            // Because we wrote to shared memory, we need to synchronize
            __syncthreads();
        }
        // Write to global matrix:
        if(thread_a_offset == 0){
            SL_matrix[ColMajor_index_2D(i, j, DL_SL_n_rows)] =
                    -SL_sum[ColMajor_index_2D(local_i,
                                              0, blockDim.x)];

            DL_matrix[ColMajor_index_2D(i, j, DL_SL_n_rows)] =
                    DL_sum[ColMajor_index_2D(local_i,
                                             0, blockDim.x)];
        }
    }
}

// @sect5{Wrapper Function: assemble_bem_matrices}

// Wrapper function around the kernel calls.\n
// This function decides the block and grid sizes for CUDA and calls the appropriate kernels.
// The parameters are the pointers to device arrays allocated in the CUDA-driver constructor, together
// the sizes of the matrices and vectors. These sizes are needed because device arrays are simple C-arrays
// and so they do not have information on their one size.
template<int dim>
void step28::BEMKernels<dim>::assemble_bem_matrices(SciPAL::ShapeData<double> DL_matrix,
                                                    SciPAL::ShapeData<double> SL_matrix,
                                                    const SciPAL::ShapeData<double> x, // support points
                                                    const SciPAL::ShapeData<double> q, // q points
                                                    const SciPAL::ShapeData<double> n, // normals
                                                    const SciPAL::ShapeData<double> W  // matrix with trial function values * JxW
                                                    // old list of args:
                                                    /*double * DL_matrix, double * SL_matrix,
                                                    const int DL_SL_n_rows,
                                                    const int DL_SL_n_cols,
                                                    const double * x, const int x_size,
                                                    const double * q, const int q_size,
                                                    const double * n, const int n_size,
                                                    const double * W, const int W_n_rows, const int W_n_cols*/)
{

    const uint q_size = q.n_rows/dim;

    const uint x_size = x.n_rows/dim;

    if(q_size > NUM_THREADS){
        fprintf(stderr, "To be implemented: number of quad points (%d() > NUMTHREADS(%d)",
                q_size, NUM_THREADS);
        return;
    }
    //Range of i: 0--x_size-1 \n
    //Range of j: 0--W_n_cols-1 \n
    //Range of a: 0--q_size-1 \n

    // We call different kernels, depending on whether we are dealing with a power of two
    // here. The kernel for a power of two is significantly faster, but only works
    // for powers of two
    if(is_power_of_two(q_size)){

        // If @p q_size (the size of $\mathbf{x}_a$) is a power of two,
        // then we allocate one thread for each entry of $G_{ia}$ that we have to calculate.

        //\htmlimage{visu_Po2.png, 350, Parallization strategy\, here for the single-layer operator $V_{h\,ij}$}

        // Then each thread calculates its entry $G_{ia}$, multiplies it with the
        // corresponding entry of $w_{aj}$ and writes it to a shared memory matrix.
        // Because that matrix has a power of two as its number of columns, the sum can be calculated
        // with a simple reduction, adding the second half-block onto the first half-block and then
        // halfing the blocksize in each step.\n
        // Lastly, there is a @p for() loop over all the @p j's.
        uint blocksize_x = q_size;
        uint blocksize_y = NUM_THREADS / blocksize_x;

        uint n_blocks = (x_size + blocksize_y - 1)/blocksize_y;

        dim3 block_size(blocksize_x, blocksize_y);
        dim3 grid_size (1, n_blocks);

        __assemble_bem_matrices_Po2<dim><<<grid_size,block_size>>>(DL_matrix, SL_matrix,
                                                                   // DL_SL_n_rows,
                                                                   // DL_SL_n_cols,
                                                                   x, // x_size,
                                                                   q, // q_size,
                                                                   n, // n_size,
                                                                   W //, W_n_rows, W_n_cols
                                                                   );
    } else {
        // If @p q_size (the size of $\mathbf{x}_a$) is a not a power of two,
        // the reduction strategy outlined above does not work. So here, instead of
        // allocating one thread for each entry of $G_{ia}$, one block of size
        // @p num_is*4 is allocated and the entries of $G_{ia}$ are calculated in
        // blocks of 4.
        // 4 is chosen here, because for sensible values of the quadrature order, the resulting
        // number of quadrature points (q_size) is always divisible by 4.
        //\htmlimage{visu_allg.png, 350, Parallization strategy\, here for the single-layer operator $V_{h\,ij}$}
        // So the parallelization strategy here is: pre-calculate the needed values of $G_{ia}$
        // in each thread (so the first thread calculates $a=0,4,8\dots$, the second calculates $a=1,5,9\dots$).

        // Then each thread multiplies it's values of $G_{ia}$ with the corresponding values of $w_{aj}$ and
        // adds up those partial sums. Now we only have to add up four values to calculate an entry
        // of e.g. $V_{h,ij}$ which is easily done using the reduction strategy above, since 4 is a power
        // of 2.\n
        // Lastly, we have a @p for() loop over j again, to calculate all values $i,j$.
        uint blocksize_y = BLOCK_DIM_Y;
        uint blocksize_x = NUM_THREADS/blocksize_y;

        uint n_blocks = (x_size + blocksize_x - 1)/blocksize_x;

        dim3 block_size(blocksize_x, blocksize_y);
        dim3 grid_size (n_blocks, 1);

#ifdef wergwerwerwerf
        __assemble_bem_matrices<dim><<<grid_size,block_size>>>(DL_matrix, SL_matrix,
                                                               DL_SL_n_rows,
                                                               DL_SL_n_cols,
                                                               x,x_size,
                                                               q, q_size,
                                                               n, n_size,
                                                               W, W_n_rows, W_n_cols
                                                               );
#endif
    }
    cudaDeviceSynchronize();
}

// Finally, we have to specialize the templates to force the compiler to actually compile something.
// This has to be at the end of file, because all functions have to be declared and their bodies defined
// before the class can be explictly instantiated by the compiler.
//template class step28::BEMKernels<2>;
template class step28::BEMKernels<3>;
