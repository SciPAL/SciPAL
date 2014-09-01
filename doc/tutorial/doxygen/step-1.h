/**
 * @page step_1                 The step-1 tutorial program
@htmlonly
<table class="tutorial" width="100%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
      </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ParallelizationofCholeskyFactorization">Parallelization of Cholesky Factorization</a>
        <li><a href="#ClassKernels">Class: Kernels</a>
      <ul>
        <li><a href="#ClassCholesky">Class: Cholesky</a>
        <li><a href="#ClassLU">Class: LU</a>
      </ul>
        <li><a href="#DeviceFunctions">Device Functions</a>
      <ul>
        <li><a href="#DeviceFunctionlex_index_2D">Device Function: lex_index_2D</a>
        <li><a href="#DeviceFunctionglobal_pos">Device Function: global_pos</a>
        <li><a href="#DeviceFunctioninv_sqrt">Device Function: inv_sqrt</a>
      </ul>
        <li><a href="#CholeskyKernels">Cholesky Kernels</a>
      <ul>
        <li><a href="#Kernel__single_thread">Kernel: __single_thread</a>
        <li><a href="#Kernelfactorize_diag_block">Kernel: factorize_diag_block</a>
        <li><a href="#Kernelstrip_update">Kernel: strip_update</a>
        <li><a href="#Kerneldiag_update">Kernel: diag_update</a>
        <li><a href="#Kernello_update">Kernel: lo_update</a>
      </ul>
        <li><a href="#Wrapperfunctions">Wrapper functions</a>
      <ul>
        <li><a href="#Functionblackbox">Function: blackbox</a>
        <li><a href="#Functionsingle_thread">Function: single_thread</a>
        <li><a href="#Functionfactorize_diag_block">Function: factorize_diag_block</a>
        <li><a href="#Functionstrip_update">Function: strip_update</a>
        <li><a href="#Functiondiag_update">Function: diag_update</a>
        <li><a href="#Functionlo_update">Function: lo_update</a>
      </ul>
        <li><a href="#ClassCUDADriver">Class: CUDADriver</a>
      <ul>
        <li><a href="#ConstructorCUDADriver">Constructor: CUDADriver</a>
        <li><a href="#Functionfactorize">Function: factorize</a>
        <li><a href="#Functionchol_fac">Function: chol_fac</a>
        <li><a href="#Functionsingle_thread_cholesky">Function: single_thread_cholesky</a>
      </ul>
        <li><a href="#ClassSimParams">Class: SimParams</a>
      <ul>
        <li><a href="#Functiondeclare">Function: declare</a>
        <li><a href="#Functionget">Function: get</a>
      </ul>
        <li><a href="#ClassCholeskyTest">Class: CholeskyTest</a>
        <li><a href="#ClassCholesky">Class: Cholesky</a>
      <ul>
        <li><a href="#ConstructorCholeskyTest">Constructor: CholeskyTest</a>
        <li><a href="#Functionsetup_and_assemble_test_matrix">Function: setup_and_assemble_test_matrix</a>
        <li><a href="#Functionrun">Function: run</a>
        <li><a href="#Functioncpu_tiled">Function: cpu_tiled</a>
        <li><a href="#FunctionLLtMult">Function: LLtMult</a>
      </ul>
        <li><a href="#ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#Functionprecision_id">Function: precision_id</a>
        <li><a href="#Functionrun">Function: run</a>
        <li><a href="#Funktionmain">Funktion: main</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
      </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
    <ul>
        <li><a href="#plain-ParallelizationofCholeskyFactorization">Parallelization of Cholesky Factorization</a>
        <li><a href="#plain-ClassKernels">Class: Kernels</a>
      <ul>
        <li><a href="#plain-ClassCholesky">Class: Cholesky</a>
        <li><a href="#plain-ClassLU">Class: LU</a>
      </ul>
        <li><a href="#plain-DeviceFunctions">Device Functions</a>
      <ul>
        <li><a href="#plain-DeviceFunctionlex_index_2D">Device Function: lex_index_2D</a>
        <li><a href="#plain-DeviceFunctionglobal_pos">Device Function: global_pos</a>
        <li><a href="#plain-DeviceFunctioninv_sqrt">Device Function: inv_sqrt</a>
      </ul>
        <li><a href="#plain-CholeskyKernels">Cholesky Kernels</a>
      <ul>
        <li><a href="#plain-Kernel__single_thread">Kernel: __single_thread</a>
        <li><a href="#plain-Kernelfactorize_diag_block">Kernel: factorize_diag_block</a>
        <li><a href="#plain-Kernelstrip_update">Kernel: strip_update</a>
        <li><a href="#plain-Kerneldiag_update">Kernel: diag_update</a>
        <li><a href="#plain-Kernello_update">Kernel: lo_update</a>
      </ul>
        <li><a href="#plain-Wrapperfunctions">Wrapper functions</a>
      <ul>
        <li><a href="#plain-Functionblackbox">Function: blackbox</a>
        <li><a href="#plain-Functionsingle_thread">Function: single_thread</a>
        <li><a href="#plain-Functionfactorize_diag_block">Function: factorize_diag_block</a>
        <li><a href="#plain-Functionstrip_update">Function: strip_update</a>
        <li><a href="#plain-Functiondiag_update">Function: diag_update</a>
        <li><a href="#plain-Functionlo_update">Function: lo_update</a>
      </ul>
        <li><a href="#plain-ClassCUDADriver">Class: CUDADriver</a>
      <ul>
        <li><a href="#plain-ConstructorCUDADriver">Constructor: CUDADriver</a>
        <li><a href="#plain-Functionfactorize">Function: factorize</a>
        <li><a href="#plain-Functionchol_fac">Function: chol_fac</a>
        <li><a href="#plain-Functionsingle_thread_cholesky">Function: single_thread_cholesky</a>
      </ul>
        <li><a href="#plain-ClassSimParams">Class: SimParams</a>
      <ul>
        <li><a href="#plain-Functiondeclare">Function: declare</a>
        <li><a href="#plain-Functionget">Function: get</a>
      </ul>
        <li><a href="#plain-ClassCholeskyTest">Class: CholeskyTest</a>
        <li><a href="#plain-ClassCholesky">Class: Cholesky</a>
      <ul>
        <li><a href="#plain-ConstructorCholeskyTest">Constructor: CholeskyTest</a>
        <li><a href="#plain-Functionsetup_and_assemble_test_matrix">Function: setup_and_assemble_test_matrix</a>
        <li><a href="#plain-Functionrun">Function: run</a>
        <li><a href="#plain-Functioncpu_tiled">Function: cpu_tiled</a>
        <li><a href="#plain-FunctionLLtMult">Function: LLtMult</a>
      </ul>
        <li><a href="#plain-ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#plain-Constructor">Constructor</a>
        <li><a href="#plain-Functionprecision_id">Function: precision_id</a>
        <li><a href="#plain-Functionrun">Function: run</a>
        <li><a href="#plain-Funktionmain">Funktion: main</a>
      </ul>
      </ul>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Introduction"></a><h1>Introduction</h1>

<p>
In this example program we discuss how to implement the Cholesky factorization of dense matrices.
Once a matrix is factorized the corresponding linear algebraic system can be solved by solving two
auxiliary triangular linear systems.
</p>
<p>
The starting point is a set of linear equations
\f{eqnarray}
 \label{LSYS}
\begin{array}{rcccr}
a_{0,0}x_0   & + & \ldots & + & a_{0,n-1}x_{n-1} \\
             &   & \ldots &   &                  \\ 
             &   & \ldots &   &                  \\
             &   & \ldots &   &                  \\ 
a_{n-1,0}x_0 & + & \ldots & + & a_{n-1,n-1}x_{n-1}
\end{array}
& = &
\begin{array}{c}
f_0 \\
 \cdot\\
 \cdot\\
 \cdot\\
f_{n-1}
\end{array}
\f}
with symmetric coefficient matrix $A = (a_{ij})_{i,j=0}^{i,j=n-1}$, $a_{ij}=a_{ji}\,\forall i,j$.
The solution vector is denoted by $x=(x_j)_{j=0}^{j=n-1}$ and the right-hand side by $f=(f_i)_{i=0}^{i=n-1}$.
The factorization process yields a lower triangular matrix $L$, such that
\f{eqnarray}
LL^T & = & A \,. 
\f}
Solving the linear system given is then achieved by solving two triangular systems
of equations by forward and backward substitution.
The first step yields the auxiliary solution $y = L^Tx$.
\f{eqnarray}
Ax ~ = ~ f  \Rightarrow LL^Tx & = & f \\
 y ~ = ~ L^Tx \Rightarrow Ly & = & f \\
\Rightarrow y & = & L^{-1} f \\
\Rightarrow x & = & (L^T)^{-1} y \,. 
\f}
</p>
<p>

<a name="CholeskyFactorizationAlgorithm"></a><h2>Cholesky Factorization - Algorithm</h2>


The entries of $L$ follow from the condition $LL^T = A$, that is
\f{eqnarray}
(LL^T)_{ij} ~ = ~ \sum_{k=0}^{\min\{i,j\}} L_{ik}L^T_{kj} 
& = & 
\sum_{k=0}^{\min\{i,j\}} L_{ik}L_{jk}  ~ = ~ A_{ij} \,. 
\f}
The upper limit of the summation is due to the triangular shape of matrix $L$.
It implies that one has to start computing its elements at the uppermost diagonal element
and then has to proceed down the column.
Only then one can go over to the next column.
Therefore the computation of the diagonal elements is inherently serial. Each diagonal element
needs the entries of $L$ which are in the block above and to the left of it.
<ol>
<li>
For the top row we have
<ol>
    <li> 
        \f{eqnarray}
            L_{00} & = & \sqrt{A_{00}}
        \f}
    </li>
    <li>do in parallel : $j = 1,\ldots, n-1$ :
        \f{eqnarray}
            L_{j0} & = & \frac{1}{\sqrt{L_{00}}} A_{0j}
        \f}
    </li>
</ol>
</li>
<li>
and for the following rows $i = 1,\ldots, n-1$
    <ol>
        <li><br>
            \f{eqnarray}
            \textrm{sum} & = &  \sum_{k=0}^{i-1} L_{ik}^2 \\
            L_{ii} & = & \sqrt{A_{ii} - \textrm{sum}  }
            \f}
        </li>
        <li>do in parallel : $j = i+1,\ldots, n-1$ :
        \f{eqnarray}
        \textrm{sum} & = & \sum_{k=0}^{i-1} L_{ik}L_{jk} \\
    L_{ji} & = & \frac{1}{\sqrt{L_{ii}}} \left( A_{ij} - \textrm{sum} \right)     \,.
        \f}
        </li>
    </ol>
</ol>
This shows, that at least the computation of the off-diagonal elements can be parallelized.
Fortunately that's the part containing most of the the computational costs in
Cholesky factorization.
</p>
<a name="SpecialrequirementsduetoCUDA"></a><h2>Special requirements due to CUDA</h2>

<ul>
<li> Hide latency of memory accesses, especially those to global memory
<li> memory accesses require a certain order to be efficient (bank conflicts)</li>
<li> Kernels should load as little data from memory as possible and should
do as much computations with it as possible</li>
<li> synchronisation - always an issue in parallel programming</li>
<li> Minimize the dependencies on nvcc
</li>
</ul>
<!--
<p>
How to handle these issues is discussed further below in the source code.
</p>

-->
<a name="ProgramStructure"></a><h2>Program Structure</h2>


<p>
The diagram in the figure roughly sketches the overall class layout
and the distribution of the classes over the source files.
On the host side you have a toplevel class which manages the interaction
with the user (data input and output; file step-1.cpp).
The actual Cholesky factorization is distributed over a stack of classes and files.
This stack reflects the
hierarchical structure of the compute environment formed by CPU and GPU.
The front end is a driver class which offers a black-box function for the factorization.
The factorization routine internally delegates the work to the CUDA kernels.
To do this, the access is not direct but via wrapper functions which encapsulate the kernel calls
and the set up of thread grids and blocks. The purpose of this large amount of indirection is to
facilitate porting the program to different parallelization architectures if new, promising ones appear.
The dependence of the host side code on the device-specific one is kept at a minimum.
The bridge between the two is formed by the file cuda_kernel_wrapper_step-1.cu.h.
</p>

\image


\image html step-1-class-design.png "Distribution of classes over source files."

<a name="Literature"></a><h2>Literature</h2>


Cambridge CUDA-Course <a href="http://www.many-core.group.cam.ac.uk/archive/CUDAcourse09/">Lectures 3 and 4</a>

 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * / *
 *    Copyright 2010,2011,2012,2013 Stephan Kramer, as of 2013: Dr. Stephan Kramer
 * 
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 * 
 *        http://www.apache.org/licenses/LICENSE-2.0
 * 
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 * 
 *    This file includes
 * 
 *    - the CUDA kernels for computing a Cholesky factorization
 *      of a real, symmetric (and hopefully positive definite) matrix.
 * 
 *    - device functions needed by the kernels for mapping thread and block indices
 *      to row and column indices of the matrix which is to be factorized or
 *      to the positions in the linear array holding the matrix entries.
 * 
 *    - the definitions of the wrapper functions for the kernels. These wrapper functions
 *      are declared in a separate header file and either allow to call
 *      the kernels individually or to execute the complete factorization
 *      via the 'blackbox' function.
 * 
 *    All kernels and wrapper functions are enclosed in a namespace 'step1'
 * * /
 * #ifndef CUDA_KERNEL_STEP_1_CU_H
 * #define CUDA_KERNEL_STEP_1_CU_H
 * 
 * 
 * @endcode
 * 
 * The most interesting part of this program is the way parallelization
 * is implemented with CUDA. Therefore, we walk through the source code
 * in a bottom-up fashion. It begins with the kernels and ends with the main function.
 * A general feature of matrix operations is that large matrices get tiled
 * into smaller submatrices on which the work is done at the end.
 * Typically, the choice of the tile size depends on the hardware.
 * 

 * 
 * 
 * <a name="ParallelizationofCholeskyFactorization"></a> 
 * <h3>Parallelization of Cholesky Factorization</h3>
 * 
 * 
 * Looking at the algorithm given in the introduction, we see that the main computational
 * effort is caused by computing the auxiliary variable @p sum in the update of the
 * off-diagonal elements.
 * 
 * 
 * For an efficient parallelization we subdivide the matrix into blocks and start with
 * factorizing the upper leftmost diagonal block.
 * Afterwards we compute the blocks in the columns below the diagonal block.
 * Finally, we update the lower right part of the matrix which is still to be factorized
 * by subtracting the auxiliary variable  @p sum.
 * This happens in diag_update() and lo_update().
 * 

 * 
 * For CUDA, since the size of a warp is 32, the most efficient choice for the blocksize
 * is 16 in most cases where double precision is the type of choice for floating point numbers
 * and 32 should be the optimum for single precision.
 * #define DEFAULT_TILE_SIZE 16
 * 

 * 
 * 

 * 
 * 
 * @code
 * namespace step1 {
 * 
 * @endcode
 * 
 * As of CUDA 5.0 it is possible to pass the size of shared memory arrays via
 * template parameters into the kernels. Therefore we use a static constant
 * rather than some macro to set the value for the size of matrix tiles.
 * For CUDA 2.x and 3.x it did not seem to work.
 * For 4.x we did not try.
 * 
 * 
 * The square of this number gives the number of threads in a thread block.
 * 16 has turned out to be a fairly good compromise.
 * 
 * @code
 * static const int DEFAULT_TILE_SIZE = 16;
 * 
 * @endcode
 * 
 * 
 * <a name="ClassKernels"></a> 
 * <h3>Class: Kernels</h3>
 * 
 * 
 * The CUDA-based kernels for the parallel computation are encapsulated
 * into wrapper functions which are all collected into one structure.
 * From a performance point of view one of the more interesting things is the
 * dependence on the precision. Therefore, the kernels are templatized with respect to the
 * number type @p T.
 * This minimizes the amount of work for creating all the implementations
 * for the different number types and precisions.
 * To do this,
 * we provide explicit template-specializations of this class.
 * The number of specializations
 * is kept at a minimum by grouping the wrapper functions
 * for the different factorization methods
 * into different private internal classes. Access to the kernels is
 * via the @p lu and @p cholesky attributes.
 * 
 * @code
 * template<typename T>
 * class Kernels {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassCholesky"></a> 
 * <h4>Class: Cholesky</h4>
 * 
 * 
 * This internal class provides the interface to
 * the actual kernels of the Cholesky factorization.
 * 
 * @code
 *     struct Cholesky {
 * 
 *         void single_thread(T * A, int n_cols, int leading_dim);
 * 
 *         cudaError_t factorize_diag_block(T *A,
 *                                          int n_blocks_done, int n_cols, int leading_dim );
 * 
 *         void strip_update(T *A,
 *                           int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void diag_update(T *A,
 *                          int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void lo_update(T *A,
 *                        int n_blocks_done, int n_blocks, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void blackbox(T *A, int n_cols, int leading_dim);
 * 
 *     };
 * 
 * @endcode
 * 
 * 
 * <a name="ClassLU"></a> 
 * <h4>Class: LU</h4>
 * 
 * 
 * A possible extension of this program would be to add an
 * internal class which provides the interface to
 * LU-specific kernels. This is the reason why the wrapper
 * functions for the Cholesky kernels have been encapsulated
 * in a private, internal class.
 * 

 * 
 * 

 * 
 * 
 * @code
 * public:
 *     Cholesky cholesky;
 * 
 * };
 * 
 * } // namespace step1 END
 * 
 * #endif // CUDA_KERNEL_STEP_1_CU_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * / *
 *      This file includes
 * 
 *    - the CUDA kernels for computing a Cholesky factorization
 *      of a real, symmetric (and hopefully positive definite) matrix.
 * 
 *    - device functions needed by the kernels for mapping thread and block indices
 *      to row and column indices of the matrix which is to be factorized or
 *      to the positions in the linear array holding the matrix entries.
 * 
 *    - the definitions of the wrapper functions for the kernels. These wrapper functions
 *      are declared in a separate header file and either allow to call
 *      the kernels individually or to execute the complete factorization
 *      via the 'blackbox' function.
 * 
 *    All kernels and wrapper functions are enclosed in a namespace 'step1'
 * * /
 * #ifdef USE_CU_C
 * 
 * @endcode
 * 
 * In order to use printf() from within a kernel
 * the good old C-header has to be included.
 * 
 * @code
 * #include <stdio.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * namespace step1 {
 * 
 * @endcode
 * 
 * 
 * <a name="DeviceFunctions"></a> 
 * <h3>Device Functions</h3>
 * 
 * 
 * Before discussing the kernels we take a look at so-called device functions.
 * They execute only on the GPU and up to CUDA 3.2 are automatically inline.
 * Nowadays they may be not inline. To enforce inlining the keyword @p __forceinline__
 * get introduced.
 * 
 * 
 * 
 * <a name="DeviceFunctionlex_index_2D"></a> 
 * <h4>Device Function: lex_index_2D</h4>
 * 
 * 
 * Compute a position $r\cdot row\_length + c$ in a linear array
 * from a row index $r$ and column index $c$
 * of a matrix entry. It is assumed that the matrix is stored row-wise.
 * @param r : row index
 * @param c : column index
 * @param leading_dim : number of matrix entries per row
 * 
 * @code
 * __forceinline__
 * __device__ int lex_index_2D(int r, int c, int leading_dim)
 * {
 *     return c +  r*leading_dim;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="DeviceFunctionglobal_pos"></a> 
 * <h4>Device Function: global_pos</h4>
 * 
 * 
 * Compute a global row or column index given the Block index and thread index.
 * @param t_pos : local row or comlumn index within block
 * @param n_blocks_done : Index of the block in C-counting
 * 
 * @code
 * template<int TILE_SIZE>
 * __forceinline__
 * __device__ int global_pos(int t_pos, int n_blocks_done)
 * {
 *     return t_pos + TILE_SIZE*n_blocks_done;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="DeviceFunctioninv_sqrt"></a> 
 * <h4>Device Function: inv_sqrt</h4>
 * 
 * 
 * CUDA provides single and double-precision versions for the
 * computation of the reciprocal of a square root.
 * In order to use it in a template context we have to provide
 * our own little wrapper function for unifying the name
 * by employing C++'s polymorphism. As for the other device functions we would like
 * to enforce inlining. However, this leads to compiler errors
 * (at least when using the nvcc which ships with CUDA 5).
 * Therefore, we do not inline.
 * The function @p rsqrtf
 * is documented in the CUDA Programming Guide.
 * @param x : real to take the square root of.
 * 
 * @code
 * __device__ float inv_sqrt(float x)
 * {
 *     return rsqrtf(x);
 * }
 * 
 * 
 * __device__ double inv_sqrt(double x)
 * {
 *     return rsqrt(x);
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="CholeskyKernels"></a> 
 * <h3>Cholesky Kernels</h3>
 * 
 * 
 * Kernels run on the GPU and are global functions. They cannot be members
 * of a class but of a namespace. Therefore, in order to group them according
 * to the factorization method they are used for we put them into different naemspaces.
 * 
 * 
 * 
 * @code
 * namespace Chol {
 * 
 * @endcode
 * 
 * 
 * <a name="Kernel__single_thread"></a> 
 * <h4>Kernel: __single_thread</h4>
 * 
 * 
 * This kernel is only started once to
 * compute the factorization of the whole matrix.
 * This shows that the source of the CPU-based Cholesky factorization
 * would also run unchanged on the GPU.
 * In memory the matrix is stored as a linear array whose length is a
 * multiple of @p leading_dim. This is for improved memory troughput. For instance,
 * this allows to make rows to have a length which is a multiple of the cache line length
 * which avoids misaligned accesses.
 * Since for Cholesky factorizations a matrix must be square it does not matter
 * whether it is stored in row- or column-major format.
 * @param A : linear array containing the matrix entries in row-major order.
 * @param n_rows : length of a column
 * @param leading_dim :
 * 
 * @code
 * template<typename T>
 * __global__
 * void
 * __single_thread(T *A, const int n_rows, const int leading_dim)
 * {
 * @endcode
 * 
 * The outer loop runs over the rows <i>L</i>.
 * 
 * @code
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *     {
 * @endcode
 * 
 * Compute diagonal entry of Cholesky factor.
 * 
 * @code
 *         T sum = 0.;
 *         unsigned int idx;
 *         unsigned int idx_c;
 * @endcode
 * 
 * Sum squares of entries computed so far in this row.
 * 
 * @code
 *         for (unsigned int u = 0; u < r; ++u)
 *         {
 *             idx = lex_index_2D(r, u, leading_dim);
 *             sum += A[idx] * A[idx];
 *         }
 *         idx = lex_index_2D(r, r, leading_dim);
 *         A[idx] = sqrt(A[idx] - sum);
 * 
 * @endcode
 * 
 * Off-diagonal entries. Here, we exploit the symmetry of <i>A</i>.
 * The auxiliary variable @p sum corresponds to step 2.1 of the algorithm
 * given in the introduction.
 * 
 * @code
 *         for (unsigned int c = r+1; c < n_rows; ++c)
 *         {
 *             sum = 0.;
 * 
 *             for (unsigned int u = 0; u < r; ++u)
 *             {
 *                 idx_c = lex_index_2D(c, u, leading_dim);
 *                 idx   = lex_index_2D(r, u, leading_dim);
 *                 sum += A[idx_c]*A[idx];
 *             }
 * 
 *             idx_c = lex_index_2D(c, r, leading_dim);
 *             idx   = lex_index_2D(r, c, leading_dim);
 *             A[idx_c]  = A[idx] - sum;
 * 
 *             idx   = lex_index_2D(r, r, leading_dim);
 *             A[idx_c] /= A[idx];
 *         }
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * 
 * 
 * <a name="Kernelfactorize_diag_block"></a> 
 * <h4>Kernel: factorize_diag_block</h4>
 * 
 * 
 * This kernel factorizes a diagonal block assuming that all previous
 * diagonal blocks have already been factored.
 * 

 * 
 * In contrast to a serial implementation we hide the summation
 * of the off-diagonal elements from the factorized part
 * in the usage of the  <i>thread</i> index and in the choice of
 * synchronization points.
 * 
 * 
 * Each instance of this kernel computes one matrix entry of the Cholesky factor $L$.
 * @param A : Linear array containing the elements of matrix to factorize
 * @param n_blocks_done : distance of the block which is to be factorized
 * from the left uppermost diagonal block.
 * @param n_cols : length of a row of $A$.
 * 
 * @code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __factorize_diag_block(T *A, int n_blocks_done,
 *                        int n_cols, int leading_dim)
 * {
 * @endcode
 * 
 * In C, arrays are stored row-wise; thus the $x$ coordinate
 * of @p threadIdx indicates the column index.
 * 
 * @code
 *     int col = threadIdx.x;
 * 
 * @endcode
 * 
 * The $y$ coordinate
 * of @p threadIdx indicates the row.
 * 
 * @code
 *     int row = threadIdx.y;
 * 
 * @endcode
 * 
 * From the thread and block index we have to compute the index
 * of the matrix element this thread has to work on.
 * This is delegated to device functions.
 * 
 * @code
 *     int global_row = global_pos<TILE_SIZE>(row, n_blocks_done);
 *     int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 * @endcode
 * 
 * For matrices whose number of rows is not
 * a multiple of @p TILE_SIZE we have to
 * take care that thread do not work on non-existing matrix entries.
 * 
 * @code
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 * @endcode
 * 
 * Simplify debugging especially of index problems we provide
 * possibility for some output.
 * 
 * @code
 * #ifdef INDEX_BOUND_DEBUG
 *     if (row == 0 && col == 0)
 *     printf("%s:\n------------------------\n", __FUNCTION__);
 *     __syncthreads();
 * 
 *     printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d\n ",
 *            row, col, global_row, global_col, idx);
 * #endif
 * 
 * @endcode
 * 
 * Copy the diagonal block to <i>shared memory</i>
 * so that threads can exchange their results.
 * To avoid memory bank conflicts
 * when consecutively accessing the elements of a column we add one column.
 * This trick is discussed in the <i>CUDA Best Practices Guide</i> in the
 * chapter about multiplying a matrix with its transpose.
 * 
 * @code
 *     __shared__ T L[TILE_SIZE][TILE_SIZE+1];
 * 
 * @endcode
 * 
 * To minimize the number of accesses to global memory we copy
 * the matrix entries of the block to shared memory and synchronize
 * all threads within the block before we go on.
 * 
 * @code
 *     L[row][col]= A[idx];
 *     __syncthreads();
 * 
 *     T fac;
 * 
 * @endcode
 * 
 * Now, we can compute the entries of the Cholesky factors.
 * We have to distinguish between diagonal elements
 * i.e. $ r = c = k$,
 * elements of the uppermost row, and the rest.
 * For matrices whose number of rows is not
 * a multiple of @p TILE_SIZE we have to modify
 * the upper bound of the loop over the diagonal of
 * matrix block.
 * 
 * @code
 *     int k_max = TILE_SIZE;
 * @endcode
 * 
 * Next, figure out whether we are in the rightmost block column of the matrix.
 * Note that @p global_pos(0, bo) cannot exceed @p n_cols due to the way threads get started.
 * 
 * @code
 *     if (n_cols - global_pos<TILE_SIZE>(0, n_blocks_done) < TILE_SIZE)
 *         k_max = n_cols%TILE_SIZE;
 * 
 * #ifdef INDEX_BOUND_DEBUG
 *       printf("k_max : %d\n", k_max);
 * #endif
 * 
 *     for(int k=0; k < k_max; k++)
 *     {
 *         __syncthreads();
 * @endcode
 * 
 * Compute the inverse square root of diagonal element $1/\sqrt{A_{kk}}$
 * using the device function defined above,
 * and store the result in a register.
 * The device function is necessary to encapsulate the precision-dependent
 * function names of the reciprocal square root.
 * Precomputing $1/\sqrt{A_{kk}}$ allows to map the division by $\sqrt{A_{kk}}$
 * of the off-diagonal elements of $L$ to a multiplication which is computationally cheaper.
 * 
 * @code
 *         fac = inv_sqrt(L[k][k]);
 *         __syncthreads();
 * 
 * @endcode
 * 
 * We compute the diagonal element and the row to its right-hand side
 * \f{eqnarray}
 * L_{ck}  =  fac \cdot A_{ck} & = &
 * \left\lbrace
 * \begin{array}{ll}
 * \sqrt{A_{kk}} & c = k \\
 * \frac{A_{ck}}{\sqrt{A_{kk}}} & c > k
 * \end{array}
 * \right.
 * \f}
 * 
 * @code
 *         if ((row==k)&&(col>=k)) L[col][row]=(L[col][row])*fac;
 * 
 * @endcode
 * 
 * Again we have to synchronize.
 * 
 * @code
 *         __syncthreads();
 * 
 * @endcode
 * 
 * Next, compute the lower left triangle.
 * 

 * 
 * The off-diagonal entries follow from
 * \f{eqnarray}
 * s_k & = & L_{ck}L_{rk} \\
 * s_c & = & \sum_{k=c+1}^{TILE\_SIZE}(-s_k) \\
 * L_{rc} & = & \frac{1}{\sqrt{A_{kk}}} \left(A_{rc}
 * + s_c \right)\,.
 * \f}
 * We only have to perform the first of these steps, i.e. computing $s_k$,
 * and one addition in the sum.
 * Hence, $s_c$ is computed incrementally by the @p k loop.
 * Instead of storing  $s_c$, we subtract the individual terms  $s_k$ from $A_{rc}$
 * and store the modified  $A_{rc}$ in the Cholesky factor.
 * 
 * 
 * When the condition
 * @p row == @p k is true, @p L[row][col] contains the final value
 * $A_{rc} + s_c$. The multiplication with @p fac, that is $1/\sqrt{A_{kk}}$,
 * is performed implicitly during this process.
 * 
 * @code
 *         if ((row>=col)&&(col>k)) L[row][col] -= L[col][k]*L[row][k];
 *     }
 * 
 *     __syncthreads();
 * 
 * @endcode
 * 
 * At the end, we copy the result back to global memory.
 * 
 * @code
 *     if (row>=col) A[idx] = L[row][col];
 * 
 * 
 * #ifdef INDEX_BOUND_DEBUG
 *     printf("A_%d : %f\n", idx,  L[row][col]);
 * #endif
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Kernelstrip_update"></a> 
 * <h4>Kernel: strip_update</h4>
 * 
 * 
 * This kernel updates $L_{ij}$ in the columns below
 * diagonal block <i>j</i>.
 * 
 * 
 * @param A : global Matrix
 * @param n_blocks_done : distance from the diagonal block just factorized.
 * @param n_cols : length of a row of $A$.
 * 
 * @code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __strip_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
 * {
 * @endcode
 * 
 * @p boffy und @p boffx are the coordinates of the matrix block
 * this thread works on.
 * 
 * @code
 *     int boffy=n_blocks_done;
 * 
 * @endcode
 * 
 * The "+1" is needed since @p n_blocks_done denotes the position of the left uppermost
 * block.
 * 
 * @code
 *     int boffx = blockIdx.x + boffy + 1;
 * 
 * @endcode
 * 
 * As in the last kernel.
 * 
 * @code
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 * @endcode
 * 
 * Again avoid bank conficts by adding a column
 * 
 * @code
 *     __shared__ T topleft[TILE_SIZE][TILE_SIZE+1];
 *     __shared__ T workingmat[TILE_SIZE][TILE_SIZE+1];
 * 
 * @endcode
 * 
 * Grab the data of the diagonal block just factorized ...
 * 
 * @code
 *     int global_row = global_pos<TILE_SIZE>(row,n_blocks_done);
 *     int global_col = global_pos<TILE_SIZE>(col,n_blocks_done);
 * 
 * @endcode
 * 
 * For matrices whose number of rows is not
 * a multiple of @p TILE_SIZE we have to
 * take care that thread do not work on non-existing matrix entries.
 * 
 * @code
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     topleft[row][col]=A[idx];
 * 
 * @endcode
 * 
 * and the transposed block which is to be processed
 * 
 * @code
 *     global_row = global_pos<TILE_SIZE>(row,boffx);
 *     int idx_w = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     workingmat[col][row] = A[idx_w];
 * 
 *     __syncthreads();
 * 
 * 
 *     int k_max = TILE_SIZE;
 * 
 * @endcode
 * 
 * Do step 2.2 of the algorithm given in the introduction.
 * Each thread works on one column.
 * 
 * @code
 *     if(row==0)
 *         for (int k=0; k < k_max; k++)
 *         {
 *             T sum=0.;
 *             for (int m = 0; m < k; m++)
 *                 sum += topleft[k][m]*workingmat[m][col];
 * 
 *             workingmat[k][col] = (workingmat[k][col] - sum)/topleft[k][k];
 *         }
 * 
 *     __syncthreads();
 * 
 *     A[idx_w] = workingmat[col][row];
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Kerneldiag_update"></a> 
 * <h4>Kernel: diag_update</h4>
 * 
 * 
 * This Kernel computes the contribution of the last factorized diagonal block
 * to the auxiliary variable @p sum in step 2.1.
 * @param A : global matrix
 * @param n_blocks_done : block offset from upper left corner of the matrix.
 * @param n_cols : length of a row of $A$.
 * 
 * @code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __diag_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
 * {
 * @endcode
 * 
 * Finding the global indices and setup of the shared memory is as above.
 * 
 * @code
 *     int boffx = blockIdx.x + n_blocks_done + 1;
 * 
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 *     int global_row = global_pos<TILE_SIZE>(row, boffx);
 *     int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 * @endcode
 * 
 * For matrices whose number of rows is not
 * a multiple of @p TILE_SIZE we have to
 * take care that threads do not work on non-existing matrix entries.
 * 
 * @code
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     __shared__ T left[TILE_SIZE][TILE_SIZE+1];
 * 
 * @endcode
 * 
 * Copy and synchronize.
 * 
 * @code
 *     left[row][col]= A[idx];
 * 
 *     __syncthreads();
 * 
 * @endcode
 * 
 * The thread with index (row, col) computes the corresponding
 * term from step 2.1.
 * 
 * @code
 *     T sum = 0.f;
 * 
 * 
 *     int k_max = TILE_SIZE;
 * 
 *     if(row>=col)
 *     {
 *         for(int kk=0; kk<k_max; kk++) sum += left[row][kk]*left[col][kk];
 * 
 * @endcode
 * 
 * Subtract the result from the global Matrix entry.
 * 
 * @code
 *         global_col = global_pos<TILE_SIZE>(col, boffx);
 *         idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *         A[idx] -= sum;
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Kernello_update"></a> 
 * <h4>Kernel: lo_update</h4>
 * 
 * 
 * This kernel applies the intermediate results produced by strip_update() and diag_update()
 * to the rest of the matrix.
 * @param A : global matrix
 * @param n_blocks_done : Block-offset
 * @param n_blocks : Number of blocks, in which a row of @p A is subdivided.
 * @param n_cols : length of a row of @p A.
 * 
 * @code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __lo_update(T *A, int n_blocks_done, int n_blocks, int n_cols, int leading_dim)
 * {
 * 
 * @endcode
 * 
 * Start with local and global Indices
 * of the entry this thread works on ...
 * 
 * @code
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 *     int boffy=blockIdx.y+n_blocks_done+1;
 *     int boffx=boffy+1;
 * 
 *     __shared__ T left[TILE_SIZE][TILE_SIZE];
 * 
 * @endcode
 * 
 * The extra column is needed only for @p upt.
 * 
 * @code
 *     __shared__ T upt[TILE_SIZE][TILE_SIZE+1];
 * 
 * 
 * @endcode
 * 
 * Start reading data at the lower left.
 * 
 * @code
 *     int global_row_src = global_pos<TILE_SIZE>(row, boffy);
 *     int global_col_src = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 * @endcode
 * 
 * For matrices whose number of rows is not
 * a multiple of @p TILE_SIZE we have to
 * take care that thread do not work on non-existing matrix entries.
 * 
 * @code
 *     if ((global_row_src >= n_cols) || (global_col_src >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row_src, global_col_src, leading_dim);
 * 
 *     upt[row][col]=A[idx];
 *     __syncthreads();
 * 
 * 
 * #ifdef nSCHUR_DEBUG
 *     if (row == 0 && col == 0)
 *     printf("%s block (%d,%d):\n------------------------\n",
 *            __FUNCTION__, blockIdx.x, blockIdx.y);
 *     __syncthreads();
 * 
 * #endif
 * 
 *     for (;boffx<n_blocks;boffx++)
 *     {
 *         int global_row = global_pos<TILE_SIZE>(row, boffx);
 *         idx = lex_index_2D(global_row, global_col_src, leading_dim);
 * 
 * @endcode
 * 
 * Reset shared memory.
 * 
 * @code
 *         left[row][col]= 0.;
 *         left[row][col]=A[idx];
 * 
 * #ifdef SCHUR_DEBUG
 *         printf("loading  left[%d][%d]=A[%d] == %f\n", row, col, idx, A[idx]);
 * #endif
 *         __syncthreads();
 * 
 *         if (global_row < n_cols)
 *         {
 *             T matrixprod=0.f;
 * 
 * @endcode
 * 
 * The thread with index (row, col) computes the corresponding term from step 2.2.
 * 
 * @code
 *             int k_max = TILE_SIZE;
 * 
 *             for (int kk=0;kk<k_max;kk++)
 *             {
 *                 matrixprod+=left[row][kk]*upt[col][kk];
 * 
 * #ifdef SCHUR_DEBUG
 *                 printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f, L_r%d : %f, U_1c : %f \n ",
 *                        row, col, global_row, global_col_src, idx, matrixprod, kk, left[row][kk], upt[col][kk]);
 * #endif
 *             }
 * 
 *             int global_col = global_pos<TILE_SIZE>(col, boffy);
 * 
 *             if (global_col < n_cols)
 *             {
 *                 idx = lex_index_2D(global_row, global_col, leading_dim);
 *                 A[idx] -= matrixprod;
 * 
 * #ifdef SCHUR_DEBUG
 *                 if (row == 0 && col == 0)
 *                     printf("%s:\n------------------------\n", __FUNCTION__);
 *                 __syncthreads();
 * 
 *                 printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f\n ",
 *                        row, col, global_row, global_col, idx, matrixprod);
 * #endif
 *             }
 *         }
 *     }
 * 
 * }
 * 
 * } // namespace Chol END
 * 
 * } // namespace step1 END
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Wrapperfunctions"></a> 
 * <h3>Wrapper functions</h3>
 * 
 * 
 * The wrapper functions mainly repeat the arguments of the kernels. Besides that
 * they manage the set up of the thread blocks and grids. Due to the tiling of
 * matrix all kernels will be executed by grids of 2-dimensional thread blocks.
 * The block size reflects the tiling. The grids are always 1-dimensional
 * since the off-diagonal update operations are effectively independent
 * updates of several rows or columns at once.
 * 
 * 
 * 
 * <a name="Functionblackbox"></a> 
 * <h4>Function: blackbox</h4>
 * 
 * 
 * This function provides a complete GPU-based implementation
 * of the Cholesky factorization and is an example of how to
 * call several kernels from one wrapper function.
 * 
 * @code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::blackbox(T * a_d, int n_cols, int leading_dim)
 * {
 *     cudaError_t error;
 * 
 * @endcode
 * 
 * Compute the number of blocks needed to cover the matrix.
 * 
 * @code
 *     int n_blocks = (n_cols+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);
 * 
 * @endcode
 * 
 * A thread-block should be as large as a matrix block.
 * 
 * @code
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     dim3 logrid;
 * 
 *     for(int i=n_blocks; i>2; --i)
 *     {
 *         logrid.x=1;
 *         logrid.y=i-2;
 * 
 *         dim3 stripgrid(i-1);
 * 
 * @endcode
 * 
 * For the diagonal block we need only one block, thus the grid size is 1.
 * 
 * @code
 *         Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE>
 *                 <<<1, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__strip_update<T, DEFAULT_TILE_SIZE>
 *                 <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__diag_update<T, DEFAULT_TILE_SIZE>
 *                 <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 *         Chol::__lo_update<T, DEFAULT_TILE_SIZE>
 *                 <<< logrid, threads >>>(a_d, n_blocks-i, n_blocks, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 *     }
 * 
 * @endcode
 * 
 * For the last 2x2-Block submatrix @p lo_update() is not needed anymore.
 * 
 * @code
 *     if(n_blocks>1)
 *     {
 *         Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols,
 *                                                    leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *     }
 * 
 * @endcode
 * 
 * Factorize the last diagonal block.
 * 
 * @code
 *     Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-1, n_cols, leading_dim);
 * 
 * @endcode
 * 
 * Make all threads finish before the host is allowed to go on.
 * Remember, kernel starts are asynchronous.
 * 
 * @code
 *     cudaThreadSynchronize();
 * 
 * @endcode
 * 
 * ... check error state of GPU.
 * 
 * @code
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsingle_thread"></a> 
 * <h4>Function: single_thread</h4>
 * 
 * 
 * This function wraps the kernel call.
 * 
 * 
 * @param a_d : Pointer to device memory containing the matrix entries.
 * @param n_rows : the number of rows of the matrix. This implies the number of columns and total number of elements as the matrix behind @p a_d is supposed to be square.
 * 
 * @code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::single_thread(T * a_d, int n_cols, int leading_dim)
 * {
 * @endcode
 * 
 * Start the kernel which is supposed to work only in one thread.
 * 
 * @code
 *     Chol::__single_thread<<<1,1>>>(a_d, n_cols, leading_dim);
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionfactorize_diag_block"></a> 
 * <h4>Function: factorize_diag_block</h4>
 * 
 * 
 * 
 * @code
 * template<typename T>
 * cudaError_t step1::Kernels<T>::Cholesky::factorize_diag_block(T * a_d,
 *                                                               int n_blocks_done,
 *                                                               int n_cols,
 *                                                               int leading_dim)
 * {
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 *     Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks_done,  n_cols, leading_dim);
 *     cudaThreadSynchronize();
 * 
 *     return cudaGetLastError();
 * }
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionstrip_update"></a> 
 * <h4>Function: strip_update</h4>
 * 
 * 
 * 
 * @code
 * template<typename T>
 * void
 * step1::Kernels<T>::Cholesky::strip_update(T *a_d,
 *                                           int n_blocks_done,
 *                                           int n_remaining_blocks, int n_cols,
 *                                           int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 stripgrid(n_remaining_blocks-1);
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
 *                                                  n_blocks_done, n_cols,
 *                                                  leading_dim);
 * 
 * @endcode
 * 
 * Every update must be synchronized, because we can't continue the factorization before
 * having updated the submatrix of A.
 * 
 * @code
 *     cudaThreadSynchronize();
 * 
 * @endcode
 * 
 * The last task is to query the error state of the CUDA context.
 * 
 * @code
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 * 
 * @endcode
 * 
 * In case of an error
 * 
 * @code
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functiondiag_update"></a> 
 * <h4>Function: diag_update</h4>
 * 
 * 
 * 
 * @code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::diag_update(T *a_d,
 *                                               int n_blocks_done,
 *                                               int n_remaining_blocks,
 *                                               int n_cols,
 *                                               int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 stripgrid(n_remaining_blocks-1);
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
 *                                                 n_blocks_done, n_cols,
 *                                                 leading_dim);
 * 
 *     cudaThreadSynchronize();
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionlo_update"></a> 
 * <h4>Function: lo_update</h4>
 * 
 * 
 * 
 * @code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::lo_update(T *a_d,
 *                                             int n_blocks_done,
 *                                             int n_blocks,
 *                                             int n_remaining_blocks ,
 *                                             int n_cols,
 *                                             int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 logrid;
 *     logrid.x=1;
 *     logrid.y=n_remaining_blocks-2;
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__lo_update<T, DEFAULT_TILE_SIZE><<< logrid, threads >>>(a_d,
 *                                              n_blocks_done, n_blocks,  n_cols, leading_dim);
 *     cudaThreadSynchronize();
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * This 2-liner provides all possible template sepcializations
 * for real-valued matrices.
 * 
 * @code
 * template class step1::Kernels<float>;
 * template class step1::Kernels<double>;
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_1_H
 * #define CUDADriver_STEP_1_H
 * 
 * #include <step-1/SimParams.h>
 * 
 * #include <lac/FullMatrixAccessor.h>
 * #include <lac/blas++.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * 
 * 
 * @endcode
 * 
 * Each example program is put into a separate namespace so that it
 * is easier to combine classes from different examples.
 * 
 * @code
 * namespace step1 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassCUDADriver"></a> 
 * <h3>Class: CUDADriver</h3>
 * 
 * 
 * This class is responsible for managing host-device communication,
 * mainly data transfer and invoking kernels.
 * To simplify host-device data transfer, we use our SciPal-library
 * which encapsulates all the details. The implementation of the kernels
 * is inherited privately in order to express that
 * this class is <i>implemented with</i> the Kernels class.
 * 
 * @code
 * template<typename T>
 * class CUDADriver : private Kernels<T> {
 * 
 * public:
 * @endcode
 * 
 * Some typedefs to make life easier.
 * 
 * @code
 *     typedef typename blas_pp<T, cublas>::blas_wrapper_type BW;
 *     typedef typename blas_pp<T, cublas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, cublas>::Matrix Matrix;
 *     typedef typename blas_pp<T, cublas>::SubMatrix SubMatrix;
 *     typedef typename blas_pp<T, cublas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, cublas>::Vector Vector;
 *     typedef typename blas_pp<T, cublas>::SubColVector SubColVector;
 *     typedef typename blas_pp<T, cublas>::SubVectorBase SubVectorBase;
 * 
 *     typedef std::map<std::string,double> TimerName2Value;
 * 
 * 
 * 
 *     CUDADriver(const SimParams &p);
 * 
 * 
 * 
 *     double factorize(FullMatrixAccessor& A);
 * 
 *     double factorize(dealii::FullMatrix<T> &A);
 * 
 *     void chol_fac(FullMatrixAccessor& A, TimerName2Value& times);
 * 
 *     void lu_fac(FullMatrixAccessor& A, TimerName2Value& times);
 * 
 *     void single_thread_cholesky(FullMatrixAccessor& A);
 * 
 * 
 * private:
 *     Matrix A_d;
 * 
 *     const SimParams * params;
 * };
 * 
 * } // namespace step1 END
 * 
 * #endif // CUDADriver_STEP_1_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_1_HH
 * #define CUDA_DRIVER_STEP_1_HH
 * 
 * #include <base/CUDATimer.h>
 * 
 * #include <step-1/cuda_driver_step-1.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * 
 * 
 * #include <QTime>
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorCUDADriver"></a> 
 * <h4>Constructor: CUDADriver</h4>
 * 
 * 
 * In this case the constructor of the driver class does not do much.
 * 
 * @code
 * template<typename T>
 * step1::CUDADriver<T>::CUDADriver(const SimParams &p)
 *     :
 *       params(&p)
 * {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionfactorize"></a> 
 * <h4>Function: factorize</h4>
 * 
 * 
 * This function copies the data from the host to the device, starts the CUDA-based
 * factorization and copies the result back to the host.
 * There are two version of this function. The first one executes
 * the Cholesky decomposition as a whole.
 * @param A : Matrix to factorize
 * @return Number of seconds spent in GPU version of Cholesky factorization. NOT milliseconds.
 * 
 * @code
 * template<typename T>
 * double
 * step1::CUDADriver<T>::factorize(FullMatrixAccessor& A)
 * {
 * @endcode
 * 
 * Copy from host to device
 * 
 * @code
 *     this->A_d = A;
 * 
 *     QTime t;
 *     t.start();
 * @endcode
 * 
 * Call the multi-threaded factorization.
 * To do this, we have to dereference the Matrix object and retrieve
 * the bare pointer to the array containing the matrix entries.
 * 
 * @code
 *     this->cholesky.blackbox(this->A_d.array().val(), this->A_d.n_cols(), this->A_d.leading_dim );
 * 
 * @endcode
 * 
 * Convert milliseconds into seconds. Cf. QT-doc.
 * 
 * @code
 *     double kernel_time = t.elapsed()/1000.;
 * 
 * @endcode
 * 
 * Finally, copy the result from device back to host.
 * 
 * @code
 *     A = this->A_d;
 * 
 *     return kernel_time;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionchol_fac"></a> 
 * <h4>Function: chol_fac</h4>
 * 
 * 
 * The second version allows to
 * measure the performance of the different parts of the algorithm
 * 
 * @code
 * template<typename T>
 * void
 * step1::CUDADriver<T>::chol_fac(FullMatrixAccessor& A, TimerName2Value& times)
 * {
 * @endcode
 * 
 * Copy from host to device
 * 
 * @code
 *     this->A_d = A;
 * 
 * @endcode
 * 
 * For timing measurements we use
 * 
 * @code
 *     QTime t;
 * 
 * @endcode
 * 
 * We dereference the Matrix object and retrieve
 * the bare pointer to the linear array containing the matrix entries.
 * 
 * @code
 *     T* a_d = this->A_d.array().val();
 * 
 *     int n_rows = this->A_d.n_rows();
 * 
 *     int leading_dim = this->A_d.leading_dim;
 * 
 * 
 * @endcode
 * 
 * Compute the number of blocks needed to cover the matrix.
 * 
 * @code
 *     int n_blocks = (A.n_rows()+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);
 * 
 *     times["factorize_diagonal_block"] = 0.;
 *     times["strip_update"] = 0.;
 *     times["diag_update"] = 0.;
 *     times["lo_update"] = 0.;
 * 
 * @endcode
 * 
 * Loop over the virtual diagonal blocks of the matrix.
 * 
 * @code
 *     for(int i = n_blocks; i > 2; --i)
 *     {
 *         t.restart();
 * 
 * 
 *         cudaError_t error = this->cholesky.factorize_diag_block(a_d, n_blocks-i, n_rows, leading_dim);
 * 
 *         times["factorize_diagonal_block"]+=t.elapsed();
 * 
 *         AssertThrow(error == cudaSuccess, dealii::ExcMessage( cudaGetErrorString(error) ) );
 * 
 * 
 *         t.restart();
 *         this->cholesky.strip_update(a_d, n_blocks-i, i, n_rows, leading_dim);
 * 
 *         times["strip_update"]+=t.elapsed();
 * 
 * 
 *         t.restart();
 *         this->cholesky.diag_update(a_d, n_blocks-i, i, n_rows, leading_dim);
 * 
 *         times["diag_update"]+=t.elapsed();
 * 
 * 
 *         t.restart();
 *         this->cholesky.lo_update(a_d, n_blocks-i, n_blocks, i, n_rows, leading_dim);
 * 
 *         times["lo_update"]+=t.elapsed();
 *     }
 * 
 * @endcode
 * 
 * For the last 2x2-Block submatrix @p lo_update() is not needed anymore.
 * 
 * @code
 *     if(n_blocks>1)
 *     {
 *         this->cholesky.factorize_diag_block(a_d, n_blocks-2, n_rows, leading_dim);
 * 
 *         this->cholesky.strip_update(a_d, n_blocks-2, 2, n_rows, leading_dim);
 * 
 *         this->cholesky.diag_update(a_d, n_blocks-2, 2,  n_rows, leading_dim);
 *     }
 * 
 * @endcode
 * 
 * Factorize the last diagonal block.
 * 

 * 
 * 
 * @code
 *     std::cout << "Cholesky decomposition..." << std::endl;
 *     this->cholesky.factorize_diag_block(a_d, n_blocks-1, n_rows, leading_dim);
 * 
 * @endcode
 * 
 * Finally, copy the result from device back to host.
 * 
 * @code
 *     A = this->A_d;
 * }
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsingle_thread_cholesky"></a> 
 * <h4>Function: single_thread_cholesky</h4>
 * 
 * 
 * This factorization method uses only one thread on the GPU. This demonstrates that CUDA can
 * be basically used like C.
 * @param A : Matrix to factorize
 * 
 * @code
 * template<typename T>
 * void
 * step1::CUDADriver<T>::single_thread_cholesky(FullMatrixAccessor& A)
 * {
 * @endcode
 * 
 * Copy from host to device
 * 
 * @code
 *     this->A_d = A;
 * 
 * @endcode
 * 
 * Call the single-threaded factorization.
 * To this end, we have to dereference the Matrix object and retrieve
 * the bare pointer to the array containing the matrix entries.
 * 
 * @code
 *     this->cholesky.single_thread(this->A_d.array().val(),
 *                                  this->A_d.n_rows(), this->A_d.leading_dim );
 * 
 * @endcode
 * 
 * Finally, copy the result from device back to host.
 * 
 * @code
 *     A = this->A_d;
 * }
 * 
 * #endif // CUDA_DRIVER_STEP_1_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef SIM_PARAMETER_H
 * #define SIM_PARAMETER_H
 * 
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <QDir>
 * 
 * namespace step1 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassSimParams"></a> 
 * <h3>Class: SimParams</h3>
 * 
 * 
 * This structure contains all parameters necessary for controling
 * the global test properties, i.e. precision and what BLAS to use.
 * The documentation strings given in the declare function provide
 * a detailed documentation of the individual attributes of this class
 * and are available at run-time as well.
 * 
 * @code
 * struct SimParams {
 * 
 *     SimParams() {}
 * 
 *     static void declare(dealii::ParameterHandler & prm);
 * 
 *     void get(dealii::ParameterHandler & prm);
 * 
 *     int
 *     device,
 *     matrix_low,
 *     matrix_high,
 *     step_size,
 *     average_runs;
 * 
 *     bool use_double;
 * 
 *     QDir run_dir;
 * 
 * private:
 * @endcode
 * 
 * As usual, inhibit automatic generation of copy ctor and assignment operator.
 * 
 * @code
 *     SimParams(const SimParams& / *other* /) {}
 * 
 *     SimParams& operator= (const SimParams& / *other* /)
 *     {
 *         return *this;
 *     }
 * 
 * };
 * }
 * 
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef SIM_PARAMETER_HH
 * #define SIM_PARAMETER_HH
 * 
 * #include <deal.II/base/parameter_handler.h>
 * #include <step-1/SimParams.h>
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functiondeclare"></a> 
 * <h4>Function: declare</h4>
 * 
 * 
 * 
 * @code
 * void
 * step1::SimParams::declare(dealii::ParameterHandler & prm)
 * {
 *     prm.enter_subsection("Simulation basics");
 * 
 *     prm.declare_entry("Run directory", "./test_me",
 *                       dealii::Patterns::Anything(),
 *                       "Specify a directory where results of "
 *                       "the test are to be stored. This can be either "
 *                       "an absolute path or path relative to the directory "
 *                       "where the program has been started. The default is "
 *                       "subdir called test_me-<date> where <date> will be replaced "
 *                       "by the date at which the program has been started. "
 *                       "this simplifies keeping the projects directory clean "
 *                       "");
 * 
 *     prm.leave_subsection();
 * 
 * 
 *     prm.enter_subsection("CUDA parameters");
 * 
 * 
 *     prm.declare_entry("Device", "0",
 *                       dealii::Patterns::Integer(),
 *                       "which CUDA-enabled GPU should be used");
 * 
 *     prm.declare_entry("Shared Memory", "true",
 *                       dealii::Patterns::Bool(),
 *                       "Whether shared (true) or L1 (false) memory should be used.");
 * 
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Testcase parameters");
 * 
 *     prm.declare_entry("Double-Precision","true",
 *                       dealii::Patterns::Bool(),
 *                       "Decide between double (true) or float (false) precision");
 * 
 *     prm.declare_entry("Matrix size - lower limit", "256",
 *                       dealii::Patterns::Integer(),
 *                       "Start value of the range of matrix sizes tested by the simulation");
 * 
 * 
 *     prm.declare_entry("Matrix size - upper limit", "513",
 *                       dealii::Patterns::Integer(),
 *                       "End value of the range of matrix sizes tested by the simulation");
 * 
 *     prm.declare_entry("Matrix size - step size", "1024",
 *                       dealii::Patterns::Integer(),
 *                       "Increment for the size of the test matrices");
 * 
 *     prm.declare_entry("Average - runs", "10",
 *                       dealii::Patterns::Integer(),
 *                       "Number of runs being used for averaging");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionget"></a> 
 * <h4>Function: get</h4>
 * 
 * 
 * 
 * @code
 * void
 * step1::SimParams::get(dealii::ParameterHandler & prm)
 * {
 *     prm.enter_subsection("Simulation basics");
 * 
 *     run_dir.setPath(prm.get("Run directory").c_str());
 *     run_dir.makeAbsolute();
 * 
 *     if (!run_dir.exists())
 *         run_dir.mkpath(".");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("CUDA parameters");
 * 
 *     device        = prm.get_integer("Device");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Testcase parameters");
 * 
 * 
 *     use_double   = prm.get_bool("Double-Precision");
 * 
 *     matrix_low   = prm.get_integer("Matrix size - lower limit");
 * 
 *     matrix_high  = prm.get_integer("Matrix size - upper limit");
 * 
 *     step_size    = prm.get_integer("Matrix size - step size");
 * 
 *     average_runs = prm.get_integer("Average - runs");
 * 
 * 
 *     prm.leave_subsection();
 * 
 * 
 * }
 * 
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * @endcode
 * 
 * STL header
 * 
 * @code
 * #include <iostream>
 * #include <vector>
 * 
 * @endcode
 * 
 * QT
 * 
 * @code
 * #include <QThread>
 * #include <QTime>
 * 
 * @endcode
 * 
 * deal.II-components
 * 
 * @code
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/lac/full_matrix.h>
 * 
 * #include <QDir>
 * 
 * 
 * @endcode
 * 
 * Drivers for the GPU part.
 * Include all other header files needed.
 * 
 * @code
 * #include <step-1/SimParams.h>
 * 
 * #include <step-1/cuda_driver_step-1.h>
 * #include <step-1/cuda_driver_step-1.hh>
 * 
 * @endcode
 * 
 * Headers from SciPal.
 * 
 * @code
 * #include <lac/FullMatrixAccessor.h>
 * #include <lac/MatrixCreator.h>
 * 
 * 
 * namespace step1 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassCholeskyTest"></a> 
 * <h3>Class: CholeskyTest</h3>
 * 
 * 
 * This class is responsible for executing the Cholesky factorization in a separate thread.
 * To this end we have to inherit from QThread and overwrite the run()-Function.
 * @param Number : type of matrix entries.
 * 
 * @code
 * template<typename Number>
 * class CholeskyTest : public QThread
 * {
 * @endcode
 * 
 * We only want to print small matrices to the console.
 * 
 * @code
 *     static const unsigned int max_display_size = 20;
 * 
 * public:
 * 
 *     CholeskyTest(int n_r,
 *                  dealii::ConvergenceTable& s_table,
 *                  const SimParams &_params);
 * 
 * protected:
 * 
 *     void setup_and_assemble_test_matrix();
 * 
 *     void factorize();
 * 
 *     void check_results();
 * 
 *     void run();
 * 
 * 
 * 
 *     int n_rows;
 * 
 *     dealii::ConvergenceTable & speedup_table;
 * 
 *     dealii::FullMatrix<Number> A, L, L_T;
 * 
 * private:
 * 
 *     const SimParams * params;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClassCholesky"></a> 
 * <h3>Class: Cholesky</h3>
 * 
 * 
 * For performance comparisons we also need a CPU-based Cholesky and LU factorization.
 * 
 * @code
 * class Cholesky {
 * 
 * public:
 *     template<typename T> static void cpu(std::vector<std::vector<T> > & A);
 * 
 *     template<typename T> static void cpu_tiled(T* A, int tile_size);
 * 
 *     template<typename T> static void LLtMult(T * A, const T * L, int n_rows);
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorCholeskyTest"></a> 
 * <h4>Constructor: CholeskyTest</h4>
 * 
 * 
 * 
 * @code
 * template<typename Number>
 * CholeskyTest<Number>::CholeskyTest(int n_r,
 *                                    dealii::ConvergenceTable& s_table,
 *                                    const SimParams &_params)
 *     :
 *       n_rows(n_r),
 *       speedup_table(s_table),
 *       params(&_params)
 * {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionsetup_and_assemble_test_matrix"></a> 
 * <h4>Function: setup_and_assemble_test_matrix</h4>
 * 
 * 
 * 
 * @code
 * template<typename Number>
 * void CholeskyTest<Number>::setup_and_assemble_test_matrix()
 * {
 *     this->A.reinit(n_rows, n_rows);
 *     this->L.reinit(n_rows, n_rows);
 *     this->L_T.reinit(n_rows, n_rows);
 * 
 *     QTime t;
 * 
 *     this->speedup_table.add_value("n rows", n_rows);
 * 
 *     std::cout << "Initial Cholesky factor deal.II :" << std::endl;
 *     std::cout << "----------------------------------------------------"
 *               << std::endl;
 * 
 * @endcode
 * 
 * For debugging purposes it is useful to
 * have the possibility of factorizing a symmetric and orthogonal matrix like
 * Hadamard matrices.
 * 
 * @code
 * #ifdef USE_HADAMARD
 *     MatrixCreator::extended_hadamard(n_rows, A);
 * #else
 *     t.start();
 * 
 * @endcode
 * 
 * Because of the design of the Tmmult-function of the class FullMatrix provided by deal.II
 * we have to compute the transpose of the reference
 * Cholesky factor $L^T_{rc} = (r+2)(c+2),
 * \quad r = 0,\ldots , n_{rows}-1\,, c = r,\ldots , n_{rows}-1$
 * 
 * @code
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *         for (unsigned int c = r; c <n_rows; ++c)
 *             L(r,c) = 1e-0*(r+2)*(c+2);
 * 
 * 
 *     if(false)
 *     qDebug("Time for CPU-based setup of Cholesky factor : %f s",
 *            t.elapsed()/1000.);
 * 
 *     if ( L.n_rows() < max_display_size)
 *         L.print(std::cout, 10, 5);
 *     if ( L.n_rows() < max_display_size)
 *         std::cout << std::endl;std::cout << std::endl;
 * 
 *     {
 * @endcode
 * 
 * At this place it is instrcutive to measure the time
 * needed for copying matrices. To do this, we reinitialize the timer for measuring the copy time.
 * 
 * @code
 *         t.restart();
 *         L_T.copy_from(L);
 *         if(false)
 *         qDebug("Time for CPU-based copying of Cholesky factor : %f s",
 *                t.elapsed()/1000.);
 * 
 * @endcode
 * 
 * Then we construct the matrix which is to be factorized by computing $A = L\cdot L^T$.
 * 
 * @code
 *         t.restart();
 *         L_T.Tmmult(A, L, false);
 *         if(false)
 *         qDebug("Time for CPU-based multiplication of Cholesky factor"
 *                " : %f s",
 *                t.elapsed()/1000.);
 * @endcode
 * 
 * The Matrix @p A is backed up, as its recomputation is pretty expensive.
 * 
 * @code
 *         L_T.copy_from(A);
 *     }
 * 
 *     if ( A.n_rows() < max_display_size)
 *     {
 *         std::cout << "Matrix to factorize :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         A.print(std::cout, 10, 5);
 * 
 *         std::cout << std::endl;std::cout << std::endl;
 *     }
 * #endif
 * 
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun"></a> 
 * <h4>Function: run</h4>
 * 
 * 
 * Initialization and factorization.
 * 
 * @code
 * template<typename Number>
 * void CholeskyTest<Number>::run()
 * {
 * @endcode
 * 
 * Prepare the test data.
 * 
 * @code
 *     this->setup_and_assemble_test_matrix();
 * 
 * @endcode
 * 
 * Timer
 * 
 * @code
 *     QTime t;
 *     double cpu_time, gpu_time;
 * 
 * @endcode
 * 
 * Execute CPU-based Cholesky factorization.
 * 
 * @code
 *     FullMatrixAccessor<Number> A_h_cpu(A, true);
 *     {
 *         t.restart();
 * 
 * @endcode
 * 
 * We need the temporary object @p A_h to get access to the value
 * array of @p A.
 * 

 * 
 * 
 * @code
 *         Cholesky::cpu_tiled<Number>(A_h_cpu.val(), A_h_cpu.n_rows() );
 * 
 * 
 *         cpu_time =  t.elapsed()/1000.;
 *         if (false)
 *         qDebug("Time for CPU-based Cholesky factorization : %f s",
 *                cpu_time);
 * 
 * @endcode
 * 
 * Save timing results in a deal.II convergence table.
 * 
 * @code
 *         this->speedup_table.add_value("CPU factorization", cpu_time);
 *         this->speedup_table.set_precision("CPU factorization", 10);
 * 
 * 
 *         std::cout << "CPU-factorized Matrix (lower triangle contains "
 *                   << "transposed Cholesky factor) :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         if ( A.n_rows() < max_display_size)
 *             A_h_cpu.print(); //std::cout, 10, 5);
 * 
 *     }
 * 
 *     std::cout << std::endl;std::cout << std::endl;
 * 
 * 
 *     Assert( A.n_rows() == A.n_cols(),
 *             dealii::ExcMessage("Matrix not square! Cholesky impossible"));
 * 
 * 
 *     double kernel_time = 0;
 * 
 * @endcode
 * 
 * Compute the Cholesky factorization on the GPU.
 * 
 * @code
 *     CUDADriver<Number> run(*params);
 *     t.restart();
 *     {
 *         FullMatrixAccessor<Number> A_h(A, true);
 *         kernel_time = run.factorize(A_h);
 * 
 *         gpu_time =  t.elapsed()/1000.;
 * 
 * @endcode
 * 
 * For the impatient user dump the results to screen.
 * 
 * @code
 *         if (false)
 *         qDebug("Time for GPU-based Cholesky factorization %f"
 *                " including data transfer : %f s\n"
 *                "speed up factor factorization : %f netto : %f n_rows : %d\n",
 *                kernel_time,
 *                gpu_time,
 *                cpu_time/kernel_time,
 *                cpu_time/gpu_time,
 *                n_rows);
 * 
 * 
 *         std::cout << "GPU-factorized Matrix (lower triangle contains "
 *                   << "transposed Cholesky factor) :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         if ( A_h.n_rows() < max_display_size)
 *             A_h.print(); //std::cout, 10, 5);
 * 
 * @endcode
 * 
 * Both factorizations should lead to the same result:
 * the original matrix still remains in the upper and the Cholesky
 * factor has appeared in the lower triangle. Thus, taking the norm
 * of their difference should yield a numerical zero.
 * 
 * @code
 *         A_h_cpu -= A_h;
 * 
 *         std::cout << "difference of factorized matrices "
 *                         << " :" << std::endl;
 *               std::cout << "----------------------------------------------------"
 *                         << std::endl;
 * 
 *               if ( A_h_cpu.n_rows() < max_display_size)
 *                   A_h_cpu.print(); //std::cout, 10, 5);
 * 
 *               double F_norm =  A_h_cpu.frobenius_norm();
 * 
 *         std::cout << "||A_cpu - A_d ||_F = " << F_norm << "\n"
 *                   << "||A_cpu - A_d ||_F/n_el = " << F_norm/A_h_cpu.n_elements() << "\n"
 *                   << "||A_cpu - A_d ||_F/||A_d||_F = " << F_norm/A_h.frobenius_norm() << std::endl;
 * 
 *     }
 * 
 *     this->speedup_table.add_value("pure GPU fac", kernel_time);
 *     this->speedup_table.set_precision("pure GPU fac", 10);
 * 
 *     this->speedup_table.add_value("GPU fac incl data transfer", gpu_time);
 *     this->speedup_table.set_precision("GPU fac incl data transfer", 10);
 * 
 * @endcode
 * 
 * Timing of individual components of the factorization in order
 * to compare manual version against CUBLAS-based variant
 * 
 * @code
 *     FullMatrixAccessor<Number> A_h(A, true);
 *     FullMatrixAccessor<Number> A_original = A_h;
 * 
 *     {
 *         typename CUDADriver<Number>::TimerName2Value times;
 * 
 *         run.chol_fac(A_h, times);
 * 
 * 
 *         typename CUDADriver<Number>::TimerName2Value::const_iterator
 *                 e=times.begin(),
 *                 end_t=times.end();
 * 
 *         for( ; e != end_t ; ++e)
 *         {
 *             this->speedup_table.add_value(e->first, e->second);
 *             this->speedup_table.set_precision(e->first,10);
 *         }
 *     }
 *     return;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functioncpu_tiled"></a> 
 * <h4>Function: cpu_tiled</h4>
 * 
 * 
 * @p T must be float or double.
 * 
 * 
 * @param A : Matrix to factorize. factorization is in-place.
 * Entries must be in the upper triangle including the diagonal.
 * 
 * 
 * The Cholesky factor will be stored in the lower triangle and overwrites the diagonal.
 * This function is identical to the Chol::__singe_thread() kernel.
 * 
 * @code
 * template<typename T>
 * void Cholesky::cpu_tiled(T* A,
 *                          int tile_size)
 * {
 * 
 *     for (int r = 0; r < tile_size; ++r)
 *     {
 * @endcode
 * 
 * Compute diagonal entry.
 * 
 * @code
 *         T sum = 0.;
 *         int idx;
 *         int idx_c;
 * 
 *         for (int u = 0; u < r; ++u)
 *         {
 *             idx = r*tile_size + u;
 *             sum += A[idx] * A[idx];
 *         }
 *         idx = r*tile_size + r;
 *         A[idx] = sqrt(A[idx] - sum);
 * 
 *         for (int c = r+1; c < tile_size; ++c)
 *         {
 *             T tmp = 0.;
 * 
 *             for (int u = 0; u < r; ++u)
 *             {
 *                 idx_c = c*tile_size + u;
 *                 idx   = r*tile_size + u;
 *                 tmp += A[idx_c]*A[idx];
 *             }
 * 
 *             idx_c = c*tile_size + r;
 *             idx   = r*tile_size + c;
 *             A[idx_c]  = A[idx] - tmp;
 *             A[idx_c] /= A[r*tile_size + r];
 *         }
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="FunctionLLtMult"></a> 
 * <h4>Function: LLtMult</h4>
 * 
 * 
 * Computes the original matrix.
 * @param A : Pointer to value array of matrix
 * @param L : Pointer to value array of Cholesky factor
 * @param n_rows : Number of rows
 * 
 * @code
 * template<typename T>
 * void Cholesky::LLtMult(T * A, const T * L, int n_rows)
 * {
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *         for (unsigned int c = 0; c <=r; ++c)
 *         {
 *             unsigned int idx = c + (r*(r+1))/2;
 *             unsigned int k_max = std::min(r,c);
 * 
 *             A[idx] = 0.;
 * 
 *             for (unsigned int k = 0; k < k_max; ++k)
 *             {
 *                 unsigned int idx_k   = k + (r*(r+1))/2;
 *                 unsigned int idx_k_T = k + (c*(c+1))/2;
 * 
 *                 A[idx] += L[idx_k]*L[idx_k_T];
 *             }
 *         }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClassMyFancySimulation"></a> 
 * <h3>Class: MyFancySimulation</h3>
 * 
 * 
 * The final class which drives the simulation.
 * This class is primarily intended
 * to manage the user's input.
 * 
 * @code
 * template<typename Number>
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     MyFancySimulation(SimParams &p);
 * 
 *     void run();
 * 
 *     static std::string precision_id();
 * 
 * private:
 *     const SimParams * params;
 * 
 * };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 
 * 
 * The constructor of the simulation class sets the pointer to the runtime parameters
 * and which GPU ("device") to use.
 * 
 * @code
 * template <typename Number>
 * step1::MyFancySimulation<Number>::MyFancySimulation(SimParams &p)
 *     :
 *       params(&p)
 * {
 *     cudaSetDevice(params->device); 
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionprecision_id"></a> 
 * <h4>Function: precision_id</h4>
 * 
 * 
 * Returns a string for identifying the precision.
 * 
 * @code
 * template<>
 * std::string MyFancySimulation<float>::precision_id()
 * {
 *     return "float";
 * }
 * 
 * template<>
 * std::string MyFancySimulation<double>::precision_id()
 * {
 *     return "double";
 * }
 * 
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun"></a> 
 * <h4>Function: run</h4>
 * 
 * 
 * Compute the factorization for different matrix sizes.
 * 
 * @code
 * template<typename Number>
 * void step1::MyFancySimulation<Number>::run()
 * {   
 * 
 * @endcode
 * 
 * Before we start the test we setup the naem of the file
 * which at the end contains the runtimes.
 * 
 * @code
 *     std::ostringstream filename;
 * 
 *     filename << "chol_fac_times_" << params->matrix_low << "_" << precision_id().c_str() << ".dat";
 * 
 * 
 * @endcode
 * 
 * The results of the factorization are stored in a convergence table.
 * 
 * @code
 *     dealii::ConvergenceTable factorization_times;
 * 
 * @endcode
 * 
 * Loop over all matrix sizes in the specified range.
 * 
 * @code
 *     for (int n = params->matrix_low; n < params->matrix_high; n+=params->step_size)
 *     {
 *         CholeskyTest<Number> driver(n, factorization_times, *params);
 * 
 *         driver.start();
 * @endcode
 * 
 * For debugging purposes it is sometimes useful to disable
 * the inheritance from QThread and to call the run()-function directly.
 * 
 * @code
 *         / * driver.run();* /
 *         driver.wait();
 * 
 * @endcode
 * 
 * To avoid data loss we save the results after each factorizatrion.
 * 
 * @code
 *         std::ofstream out(filename.str().c_str());
 *         factorization_times.write_text(out);
 *     }
 * 
 *     std::cout << "Done." << std::endl;
 * }
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Funktionmain"></a> 
 * <h4>Funktion: main</h4>
 * 
 * 
 * Instantiate and execute the simulation.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *     using namespace step1;
 * 
 *     SimParams params;
 * 
 * @endcode
 * 
 * First we declare the parameters to expect ...
 * 
 * @code
 *     dealii::ParameterHandler prm_handler;
 * 
 * @endcode
 * 
 * Get the current working directory
 * 
 * @code
 *     QDir cwd = QDir::current();
 * 
 * @endcode
 * 
 * and backup the location where the program has been started.
 * Here, we assume the we use QTCreators shadow-biuld mechanism
 * whoch puts the build directory at the same level as the directory <i>@ref step_1 "step-1"</i>
 * containing the source code.
 * 
 * @code
 *     const QDir launch_dir = cwd;
 *     cwd.setPath("../step-1");
 * 
 * @endcode
 * 
 * By default, the parameter file has the same name as the binary
 * and is supposed to be in a subdirectory prm of the directory,
 * where the program has been started.
 * 
 * @code
 *     std::string prm_filename;
 *     if (argc == 1)
 *     {
 *         std::string tmp = argv[0];
 *         int found=tmp.find_last_of('/');
 *         prm_filename = tmp.substr(found+1);
 *         prm_filename += "-Decomp.prm";
 * 
 *         cwd.setPath("./prm");
 *     }
 *     else
 *     {
 *         QFileInfo tmp(argv[1]);
 * 
 * @endcode
 * 
 * Subdivide the given filename into its path and filename
 * so that the corresponding subdirectories can be created.
 * 
 * @code
 *         QString prm_path = tmp.absolutePath();
 *         cwd.setPath(prm_path);
 *         cwd.makeAbsolute();
 *         prm_filename = tmp.fileName().toStdString();
 * 
 *         std::cout << "chosen prm file : " << tmp.absoluteFilePath().toStdString().c_str() << std::endl;
 *     }
 * 
 * @endcode
 * 
 * Before the parameter file can be read, we have to make sure that
 * its directory exists
 * 
 * @code
 *     if (!cwd.exists() )
 *         launch_dir.mkpath( cwd.absolutePath() );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 *     SimParams::declare(prm_handler);
 *     prm_handler.read_input (prm_filename);
 * 
 *     QDir::setCurrent(launch_dir.absolutePath());
 * 
 *     params.get(prm_handler);
 * 
 * @endcode
 * 
 * Create the toplevel run directory.
 * 
 * @code
 *     cwd.setPath(params.run_dir.absolutePath());
 * @endcode
 * 
 * The following lets a directory make its own path.
 * 
 * @code
 *     if (!cwd.exists())
 *         cwd.mkpath( "." );
 * 
 * @endcode
 * 
 * Now, change to the run dir
 * 
 * @code
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     cwd.setPath("./log");
 *     cwd.makeAbsolute();
 *     if (!cwd.exists())
 *         cwd.mkpath(".");
 * 
 * @endcode
 * 
 * Create the log directory and write what has been actually read
 * into log file. Basically, this is just another parameter file
 * and can thus be used again as input to another run.
 * 
 * @code
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(prm_filename.c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 * @endcode
 * 
 * At this point the toplevel run dir must exist.
 * Thus, we can change to it without any further sanity test.
 * 
 * @code
 *     QDir::setCurrent(params.run_dir.absolutePath());
 * 
 * 
 * @endcode
 * 
 * Now, run the comparison of GPU vs. CPU for the selected precision.
 * 
 * @code
 *     if (!params.use_double) {
 * 
 *         MyFancySimulation<float> machma_float(params);
 * 
 *         machma_float.run();
 *     }
 *     else {
 * 
 *         MyFancySimulation<double> machma_double(params);
 * 
 *         machma_double.run();
 *     }
 * }
 * 
 * 
 * 
 * @endcode
<a name="Results"></a><h1>Results</h1>

<p>
To generate figures comparable to the ones shown below start step-1
with the parameter file <i>step-1-Cholesky.prm</i> from the <i>prm</i> subdirectory.
A gnuplot script <i>chol_fac_plot_results.gp</i> for generating the graphical output
is in the <i>scripts</i> subdirectory. The plots are generated by switching to the
directory where the results file is stored and starting gnuplot from that directory with the plot script as argument.
Depending on the particular parameters you will have
to modify the name of the data file read by that script.
</p>
<p>
The initial comparisons for the Cholesky factorization were done on a Macbook Pro with Core2Duo processor
running at 2.53 GHz
and NVidia's Gforce 8600m GT chip.
Neglecting the time needed for host-to-device memory transfer the break even is
at roughly 170 rows.
Including transfer times this rise to 700 rows.
The Cholesky factorization does not take into account any information about
the size of the matrix entries.
Therefore the choice of the testmatrix does not really matter as long as only
timing tests are performed.
Due to the fact that GPUs have a different rounding behavior, one should also
test whether for a given
matrix and right-hand side the solution obtained from
both variants is the same up to numerical accuracy.


<p>
\htmlimage{cholesky_macbookpro_c2d_2.53GHz_8600m_GT.png, 1024, Comparison of GPU vs. CPU performance for various matrix sizes.}
</p>

</p>
<p>
The next more serious comparison is run on a system consisting of
a quad-core Xeon with 2.26 GHz and a GeForce GTX 460.

\htmlimage{chol_fac_times_256_float.png, 1024, Execution times on quad-core Xeon 2.26 GHz and gtx 460.}
\htmlimage{ chol_fac_speedup.png, 1024, Speedup of GPU over CPU on quad-core Xeon 2.26 GHz and gtx 460.}
\htmlimage{chol_fac_gpu_individual_components.png, 1024, Individual timing of the different kernels of the Cholesky factorization on gtx 460.}

</p>
 * <a name="PlainProg"></a>
 * <h1> The plain program</h1>
 * 
 * (If you are looking at a locally installed CUDA HPC Praktikum version, then the
 * program can be found at <i>
 *  .. /.. /testsite / /step-1 /step-cu.cc
 * </i>. Otherwise, this is only
 * the path on some remote server.)
 @code

 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * / *
 *    Copyright 2010,2011,2012,2013 Stephan Kramer, as of 2013: Dr. Stephan Kramer
 * 
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 * 
 *        http://www.apache.org/licenses/LICENSE-2.0
 * 
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 * 
 *    This file includes
 * 
 *    - the CUDA kernels for computing a Cholesky factorization
 *      of a real, symmetric (and hopefully positive definite) matrix.
 * 
 *    - device functions needed by the kernels for mapping thread and block indices
 *      to row and column indices of the matrix which is to be factorized or
 *      to the positions in the linear array holding the matrix entries.
 * 
 *    - the definitions of the wrapper functions for the kernels. These wrapper functions
 *      are declared in a separate header file and either allow to call
 *      the kernels individually or to execute the complete factorization
 *      via the 'blackbox' function.
 * 
 *    All kernels and wrapper functions are enclosed in a namespace 'step1'
 * * /
 * #ifndef CUDA_KERNEL_STEP_1_CU_H
 * #define CUDA_KERNEL_STEP_1_CU_H
 * 
 * 
 * 
@endcode
 <a name="plain-ParallelizationofCholeskyFactorization"></a>
@code
 * 
 * 
 * 
 * namespace step1 {
 * 
 * static const int DEFAULT_TILE_SIZE = 16;
 * 
@endcode
 <a name="plain-ClassKernels"></a>
@code
 * template<typename T>
 * class Kernels {
 * 
@endcode
 <a name="plain-ClassCholesky"></a>
@code
 *     struct Cholesky {
 * 
 *         void single_thread(T * A, int n_cols, int leading_dim);
 * 
 *         cudaError_t factorize_diag_block(T *A,
 *                                          int n_blocks_done, int n_cols, int leading_dim );
 * 
 *         void strip_update(T *A,
 *                           int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void diag_update(T *A,
 *                          int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void lo_update(T *A,
 *                        int n_blocks_done, int n_blocks, int n_remaining_blocks, int n_cols, int leading_dim);
 * 
 *         void blackbox(T *A, int n_cols, int leading_dim);
 * 
 *     };
 * 
@endcode
 <a name="plain-ClassLU"></a>
@code
 * 
 * 
 * public:
 *     Cholesky cholesky;
 * 
 * };
 * 
 * } // namespace step1 END
 * 
 * #endif // CUDA_KERNEL_STEP_1_CU_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * / *
 *      This file includes
 * 
 *    - the CUDA kernels for computing a Cholesky factorization
 *      of a real, symmetric (and hopefully positive definite) matrix.
 * 
 *    - device functions needed by the kernels for mapping thread and block indices
 *      to row and column indices of the matrix which is to be factorized or
 *      to the positions in the linear array holding the matrix entries.
 * 
 *    - the definitions of the wrapper functions for the kernels. These wrapper functions
 *      are declared in a separate header file and either allow to call
 *      the kernels individually or to execute the complete factorization
 *      via the 'blackbox' function.
 * 
 *    All kernels and wrapper functions are enclosed in a namespace 'step1'
 * * /
 * #ifdef USE_CU_C
 * 
 * #include <stdio.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * namespace step1 {
 * 
@endcode
 <a name="plain-DeviceFunctions"></a>
@code
@endcode
 <a name="plain-DeviceFunctionlex_index_2D"></a>
@code
 * __forceinline__
 * __device__ int lex_index_2D(int r, int c, int leading_dim)
 * {
 *     return c +  r*leading_dim;
 * }
 * 
 * 
@endcode
 <a name="plain-DeviceFunctionglobal_pos"></a>
@code
 * template<int TILE_SIZE>
 * __forceinline__
 * __device__ int global_pos(int t_pos, int n_blocks_done)
 * {
 *     return t_pos + TILE_SIZE*n_blocks_done;
 * }
 * 
 * 
@endcode
 <a name="plain-DeviceFunctioninv_sqrt"></a>
@code
 * __device__ float inv_sqrt(float x)
 * {
 *     return rsqrtf(x);
 * }
 * 
 * 
 * __device__ double inv_sqrt(double x)
 * {
 *     return rsqrt(x);
 * }
 * 
@endcode
 <a name="plain-CholeskyKernels"></a>
@code
 * namespace Chol {
 * 
@endcode
 <a name="plain-Kernel__single_thread"></a>
@code
 * template<typename T>
 * __global__
 * void
 * __single_thread(T *A, const int n_rows, const int leading_dim)
 * {
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *     {
 *         T sum = 0.;
 *         unsigned int idx;
 *         unsigned int idx_c;
 *         for (unsigned int u = 0; u < r; ++u)
 *         {
 *             idx = lex_index_2D(r, u, leading_dim);
 *             sum += A[idx] * A[idx];
 *         }
 *         idx = lex_index_2D(r, r, leading_dim);
 *         A[idx] = sqrt(A[idx] - sum);
 * 
 *         for (unsigned int c = r+1; c < n_rows; ++c)
 *         {
 *             sum = 0.;
 * 
 *             for (unsigned int u = 0; u < r; ++u)
 *             {
 *                 idx_c = lex_index_2D(c, u, leading_dim);
 *                 idx   = lex_index_2D(r, u, leading_dim);
 *                 sum += A[idx_c]*A[idx];
 *             }
 * 
 *             idx_c = lex_index_2D(c, r, leading_dim);
 *             idx   = lex_index_2D(r, c, leading_dim);
 *             A[idx_c]  = A[idx] - sum;
 * 
 *             idx   = lex_index_2D(r, r, leading_dim);
 *             A[idx_c] /= A[idx];
 *         }
 *     }
 * }
 * 
 * 
@endcode
 <a name="plain-Kernelfactorize_diag_block"></a>
@code
 * 
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __factorize_diag_block(T *A, int n_blocks_done,
 *                        int n_cols, int leading_dim)
 * {
 *     int col = threadIdx.x;
 * 
 *     int row = threadIdx.y;
 * 
 *     int global_row = global_pos<TILE_SIZE>(row, n_blocks_done);
 *     int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 * #ifdef INDEX_BOUND_DEBUG
 *     if (row == 0 && col == 0)
 *     printf("%s:\n------------------------\n", __FUNCTION__);
 *     __syncthreads();
 * 
 *     printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d\n ",
 *            row, col, global_row, global_col, idx);
 * #endif
 * 
 *     __shared__ T L[TILE_SIZE][TILE_SIZE+1];
 * 
 *     L[row][col]= A[idx];
 *     __syncthreads();
 * 
 *     T fac;
 * 
 *     int k_max = TILE_SIZE;
 *     if (n_cols - global_pos<TILE_SIZE>(0, n_blocks_done) < TILE_SIZE)
 *         k_max = n_cols%TILE_SIZE;
 * 
 * #ifdef INDEX_BOUND_DEBUG
 *       printf("k_max : %d\n", k_max);
 * #endif
 * 
 *     for(int k=0; k < k_max; k++)
 *     {
 *         __syncthreads();
 *         fac = inv_sqrt(L[k][k]);
 *         __syncthreads();
 * 
 *         if ((row==k)&&(col>=k)) L[col][row]=(L[col][row])*fac;
 * 
 *         __syncthreads();
 * 
 * 
 *         if ((row>=col)&&(col>k)) L[row][col] -= L[col][k]*L[row][k];
 *     }
 * 
 *     __syncthreads();
 * 
 *     if (row>=col) A[idx] = L[row][col];
 * 
 * 
 * #ifdef INDEX_BOUND_DEBUG
 *     printf("A_%d : %f\n", idx,  L[row][col]);
 * #endif
 * }
 * 
 * 
@endcode
 <a name="plain-Kernelstrip_update"></a>
@code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __strip_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
 * {
 *     int boffy=n_blocks_done;
 * 
 *     int boffx = blockIdx.x + boffy + 1;
 * 
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 *     __shared__ T topleft[TILE_SIZE][TILE_SIZE+1];
 *     __shared__ T workingmat[TILE_SIZE][TILE_SIZE+1];
 * 
 *     int global_row = global_pos<TILE_SIZE>(row,n_blocks_done);
 *     int global_col = global_pos<TILE_SIZE>(col,n_blocks_done);
 * 
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     topleft[row][col]=A[idx];
 * 
 *     global_row = global_pos<TILE_SIZE>(row,boffx);
 *     int idx_w = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     workingmat[col][row] = A[idx_w];
 * 
 *     __syncthreads();
 * 
 * 
 *     int k_max = TILE_SIZE;
 * 
 *     if(row==0)
 *         for (int k=0; k < k_max; k++)
 *         {
 *             T sum=0.;
 *             for (int m = 0; m < k; m++)
 *                 sum += topleft[k][m]*workingmat[m][col];
 * 
 *             workingmat[k][col] = (workingmat[k][col] - sum)/topleft[k][k];
 *         }
 * 
 *     __syncthreads();
 * 
 *     A[idx_w] = workingmat[col][row];
 * }
 * 
@endcode
 <a name="plain-Kerneldiag_update"></a>
@code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __diag_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
 * {
 *     int boffx = blockIdx.x + n_blocks_done + 1;
 * 
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 *     int global_row = global_pos<TILE_SIZE>(row, boffx);
 *     int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 *     if ((global_row >= n_cols) || (global_col >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *     __shared__ T left[TILE_SIZE][TILE_SIZE+1];
 * 
 *     left[row][col]= A[idx];
 * 
 *     __syncthreads();
 * 
 *     T sum = 0.f;
 * 
 * 
 *     int k_max = TILE_SIZE;
 * 
 *     if(row>=col)
 *     {
 *         for(int kk=0; kk<k_max; kk++) sum += left[row][kk]*left[col][kk];
 * 
 *         global_col = global_pos<TILE_SIZE>(col, boffx);
 *         idx = lex_index_2D(global_row, global_col, leading_dim);
 * 
 *         A[idx] -= sum;
 *     }
 * }
 * 
 * 
@endcode
 <a name="plain-Kernello_update"></a>
@code
 * template<typename T, int TILE_SIZE>
 * __global__
 * void
 * __lo_update(T *A, int n_blocks_done, int n_blocks, int n_cols, int leading_dim)
 * {
 * 
 *     int col = threadIdx.x;
 *     int row = threadIdx.y;
 * 
 *     int boffy=blockIdx.y+n_blocks_done+1;
 *     int boffx=boffy+1;
 * 
 *     __shared__ T left[TILE_SIZE][TILE_SIZE];
 * 
 *     __shared__ T upt[TILE_SIZE][TILE_SIZE+1];
 * 
 * 
 *     int global_row_src = global_pos<TILE_SIZE>(row, boffy);
 *     int global_col_src = global_pos<TILE_SIZE>(col, n_blocks_done);
 * 
 *     if ((global_row_src >= n_cols) || (global_col_src >= n_cols))
 *         return;
 * 
 *     int idx = lex_index_2D(global_row_src, global_col_src, leading_dim);
 * 
 *     upt[row][col]=A[idx];
 *     __syncthreads();
 * 
 * 
 * #ifdef nSCHUR_DEBUG
 *     if (row == 0 && col == 0)
 *     printf("%s block (%d,%d):\n------------------------\n",
 *            __FUNCTION__, blockIdx.x, blockIdx.y);
 *     __syncthreads();
 * 
 * #endif
 * 
 *     for (;boffx<n_blocks;boffx++)
 *     {
 *         int global_row = global_pos<TILE_SIZE>(row, boffx);
 *         idx = lex_index_2D(global_row, global_col_src, leading_dim);
 * 
 *         left[row][col]= 0.;
 *         left[row][col]=A[idx];
 * 
 * #ifdef SCHUR_DEBUG
 *         printf("loading  left[%d][%d]=A[%d] == %f\n", row, col, idx, A[idx]);
 * #endif
 *         __syncthreads();
 * 
 *         if (global_row < n_cols)
 *         {
 *             T matrixprod=0.f;
 * 
 *             int k_max = TILE_SIZE;
 * 
 *             for (int kk=0;kk<k_max;kk++)
 *             {
 *                 matrixprod+=left[row][kk]*upt[col][kk];
 * 
 * #ifdef SCHUR_DEBUG
 *                 printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f, L_r%d : %f, U_1c : %f \n ",
 *                        row, col, global_row, global_col_src, idx, matrixprod, kk, left[row][kk], upt[col][kk]);
 * #endif
 *             }
 * 
 *             int global_col = global_pos<TILE_SIZE>(col, boffy);
 * 
 *             if (global_col < n_cols)
 *             {
 *                 idx = lex_index_2D(global_row, global_col, leading_dim);
 *                 A[idx] -= matrixprod;
 * 
 * #ifdef SCHUR_DEBUG
 *                 if (row == 0 && col == 0)
 *                     printf("%s:\n------------------------\n", __FUNCTION__);
 *                 __syncthreads();
 * 
 *                 printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f\n ",
 *                        row, col, global_row, global_col, idx, matrixprod);
 * #endif
 *             }
 *         }
 *     }
 * 
 * }
 * 
 * } // namespace Chol END
 * 
 * } // namespace step1 END
 * 
 * 
@endcode
 <a name="plain-Wrapperfunctions"></a>
@code
@endcode
 <a name="plain-Functionblackbox"></a>
@code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::blackbox(T * a_d, int n_cols, int leading_dim)
 * {
 *     cudaError_t error;
 * 
 *     int n_blocks = (n_cols+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);
 * 
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     dim3 logrid;
 * 
 *     for(int i=n_blocks; i>2; --i)
 *     {
 *         logrid.x=1;
 *         logrid.y=i-2;
 * 
 *         dim3 stripgrid(i-1);
 * 
 *         Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE>
 *                 <<<1, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__strip_update<T, DEFAULT_TILE_SIZE>
 *                 <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__diag_update<T, DEFAULT_TILE_SIZE>
 *                 <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 *         Chol::__lo_update<T, DEFAULT_TILE_SIZE>
 *                 <<< logrid, threads >>>(a_d, n_blocks-i, n_blocks, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 *     }
 * 
 *     if(n_blocks>1)
 *     {
 *         Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols,
 *                                                    leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *         Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
 *         cudaThreadSynchronize();
 * 
 *     }
 * 
 *     Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-1, n_cols, leading_dim);
 * 
 *     cudaThreadSynchronize();
 * 
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
@endcode
 <a name="plain-Functionsingle_thread"></a>
@code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::single_thread(T * a_d, int n_cols, int leading_dim)
 * {
 *     Chol::__single_thread<<<1,1>>>(a_d, n_cols, leading_dim);
 * }
 * 
 * 
 * 
@endcode
 <a name="plain-Functionfactorize_diag_block"></a>
@code
 * template<typename T>
 * cudaError_t step1::Kernels<T>::Cholesky::factorize_diag_block(T * a_d,
 *                                                               int n_blocks_done,
 *                                                               int n_cols,
 *                                                               int leading_dim)
 * {
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 *     Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks_done,  n_cols, leading_dim);
 *     cudaThreadSynchronize();
 * 
 *     return cudaGetLastError();
 * }
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Functionstrip_update"></a>
@code
 * template<typename T>
 * void
 * step1::Kernels<T>::Cholesky::strip_update(T *a_d,
 *                                           int n_blocks_done,
 *                                           int n_remaining_blocks, int n_cols,
 *                                           int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 stripgrid(n_remaining_blocks-1);
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
 *                                                  n_blocks_done, n_cols,
 *                                                  leading_dim);
 * 
 *     cudaThreadSynchronize();
 * 
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 * 
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Functiondiag_update"></a>
@code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::diag_update(T *a_d,
 *                                               int n_blocks_done,
 *                                               int n_remaining_blocks,
 *                                               int n_cols,
 *                                               int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 stripgrid(n_remaining_blocks-1);
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
 *                                                 n_blocks_done, n_cols,
 *                                                 leading_dim);
 * 
 *     cudaThreadSynchronize();
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Functionlo_update"></a>
@code
 * template<typename T>
 * void step1::Kernels<T>::Cholesky::lo_update(T *a_d,
 *                                             int n_blocks_done,
 *                                             int n_blocks,
 *                                             int n_remaining_blocks ,
 *                                             int n_cols,
 *                                             int leading_dim)
 * {
 *     cudaError_t error;
 *     dim3 logrid;
 *     logrid.x=1;
 *     logrid.y=n_remaining_blocks-2;
 *     dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
 * 
 *     Chol::__lo_update<T, DEFAULT_TILE_SIZE><<< logrid, threads >>>(a_d,
 *                                              n_blocks_done, n_blocks,  n_cols, leading_dim);
 *     cudaThreadSynchronize();
 *     error=cudaGetLastError();
 *     if (error != cudaSuccess)
 *     {
 *         printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
 *         exit(-1);
 *     }
 * }
 * 
 * 
 * 
 * template class step1::Kernels<float>;
 * template class step1::Kernels<double>;
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_1_H
 * #define CUDADriver_STEP_1_H
 * 
 * #include <step-1/SimParams.h>
 * 
 * #include <lac/FullMatrixAccessor.h>
 * #include <lac/blas++.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * 
 * 
 * namespace step1 {
 * 
@endcode
 <a name="plain-ClassCUDADriver"></a>
@code
 * template<typename T>
 * class CUDADriver : private Kernels<T> {
 * 
 * public:
 *     typedef typename blas_pp<T, cublas>::blas_wrapper_type BW;
 *     typedef typename blas_pp<T, cublas>::FullMatrixAccessor FullMatrixAccessor;
 *     typedef typename blas_pp<T, cublas>::Matrix Matrix;
 *     typedef typename blas_pp<T, cublas>::SubMatrix SubMatrix;
 *     typedef typename blas_pp<T, cublas>::MatrixSubCol MatrixSubCol;
 *     typedef typename blas_pp<T, cublas>::Vector Vector;
 *     typedef typename blas_pp<T, cublas>::SubColVector SubColVector;
 *     typedef typename blas_pp<T, cublas>::SubVectorBase SubVectorBase;
 * 
 *     typedef std::map<std::string,double> TimerName2Value;
 * 
 * 
 * 
 *     CUDADriver(const SimParams &p);
 * 
 * 
 * 
 *     double factorize(FullMatrixAccessor& A);
 * 
 *     double factorize(dealii::FullMatrix<T> &A);
 * 
 *     void chol_fac(FullMatrixAccessor& A, TimerName2Value& times);
 * 
 *     void lu_fac(FullMatrixAccessor& A, TimerName2Value& times);
 * 
 *     void single_thread_cholesky(FullMatrixAccessor& A);
 * 
 * 
 * private:
 *     Matrix A_d;
 * 
 *     const SimParams * params;
 * };
 * 
 * } // namespace step1 END
 * 
 * #endif // CUDADriver_STEP_1_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_1_HH
 * #define CUDA_DRIVER_STEP_1_HH
 * 
 * #include <base/CUDATimer.h>
 * 
 * #include <step-1/cuda_driver_step-1.h>
 * 
 * #include <step-1/cuda_kernel_wrapper_step-1.cu.h>
 * 
 * 
 * 
 * #include <QTime>
 * 
 * 
@endcode
 <a name="plain-ConstructorCUDADriver"></a>
@code
 * template<typename T>
 * step1::CUDADriver<T>::CUDADriver(const SimParams &p)
 *     :
 *       params(&p)
 * {}
 * 
 * 
@endcode
 <a name="plain-Functionfactorize"></a>
@code
 * template<typename T>
 * double
 * step1::CUDADriver<T>::factorize(FullMatrixAccessor& A)
 * {
 *     this->A_d = A;
 * 
 *     QTime t;
 *     t.start();
 *     this->cholesky.blackbox(this->A_d.array().val(), this->A_d.n_cols(), this->A_d.leading_dim );
 * 
 *     double kernel_time = t.elapsed()/1000.;
 * 
 *     A = this->A_d;
 * 
 *     return kernel_time;
 * }
 * 
@endcode
 <a name="plain-Functionchol_fac"></a>
@code
 * template<typename T>
 * void
 * step1::CUDADriver<T>::chol_fac(FullMatrixAccessor& A, TimerName2Value& times)
 * {
 *     this->A_d = A;
 * 
 *     QTime t;
 * 
 *     T* a_d = this->A_d.array().val();
 * 
 *     int n_rows = this->A_d.n_rows();
 * 
 *     int leading_dim = this->A_d.leading_dim;
 * 
 * 
 *     int n_blocks = (A.n_rows()+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);
 * 
 *     times["factorize_diagonal_block"] = 0.;
 *     times["strip_update"] = 0.;
 *     times["diag_update"] = 0.;
 *     times["lo_update"] = 0.;
 * 
 *     for(int i = n_blocks; i > 2; --i)
 *     {
 *         t.restart();
 * 
 * 
 *         cudaError_t error = this->cholesky.factorize_diag_block(a_d, n_blocks-i, n_rows, leading_dim);
 * 
 *         times["factorize_diagonal_block"]+=t.elapsed();
 * 
 *         AssertThrow(error == cudaSuccess, dealii::ExcMessage( cudaGetErrorString(error) ) );
 * 
 * 
 *         t.restart();
 *         this->cholesky.strip_update(a_d, n_blocks-i, i, n_rows, leading_dim);
 * 
 *         times["strip_update"]+=t.elapsed();
 * 
 * 
 *         t.restart();
 *         this->cholesky.diag_update(a_d, n_blocks-i, i, n_rows, leading_dim);
 * 
 *         times["diag_update"]+=t.elapsed();
 * 
 * 
 *         t.restart();
 *         this->cholesky.lo_update(a_d, n_blocks-i, n_blocks, i, n_rows, leading_dim);
 * 
 *         times["lo_update"]+=t.elapsed();
 *     }
 * 
 *     if(n_blocks>1)
 *     {
 *         this->cholesky.factorize_diag_block(a_d, n_blocks-2, n_rows, leading_dim);
 * 
 *         this->cholesky.strip_update(a_d, n_blocks-2, 2, n_rows, leading_dim);
 * 
 *         this->cholesky.diag_update(a_d, n_blocks-2, 2,  n_rows, leading_dim);
 *     }
 * 
 * 
 *     std::cout << "Cholesky decomposition..." << std::endl;
 *     this->cholesky.factorize_diag_block(a_d, n_blocks-1, n_rows, leading_dim);
 * 
 *     A = this->A_d;
 * }
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Functionsingle_thread_cholesky"></a>
@code
 * template<typename T>
 * void
 * step1::CUDADriver<T>::single_thread_cholesky(FullMatrixAccessor& A)
 * {
 *     this->A_d = A;
 * 
 *     this->cholesky.single_thread(this->A_d.array().val(),
 *                                  this->A_d.n_rows(), this->A_d.leading_dim );
 * 
 *     A = this->A_d;
 * }
 * 
 * #endif // CUDA_DRIVER_STEP_1_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef SIM_PARAMETER_H
 * #define SIM_PARAMETER_H
 * 
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <QDir>
 * 
 * namespace step1 {
 * 
@endcode
 <a name="plain-ClassSimParams"></a>
@code
 * struct SimParams {
 * 
 *     SimParams() {}
 * 
 *     static void declare(dealii::ParameterHandler & prm);
 * 
 *     void get(dealii::ParameterHandler & prm);
 * 
 *     int
 *     device,
 *     matrix_low,
 *     matrix_high,
 *     step_size,
 *     average_runs;
 * 
 *     bool use_double;
 * 
 *     QDir run_dir;
 * 
 * private:
 *     SimParams(const SimParams& / *other* /) {}
 * 
 *     SimParams& operator= (const SimParams& / *other* /)
 *     {
 *         return *this;
 *     }
 * 
 * };
 * }
 * 
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef SIM_PARAMETER_HH
 * #define SIM_PARAMETER_HH
 * 
 * #include <deal.II/base/parameter_handler.h>
 * #include <step-1/SimParams.h>
 * 
 * 
@endcode
 <a name="plain-Functiondeclare"></a>
@code
 * void
 * step1::SimParams::declare(dealii::ParameterHandler & prm)
 * {
 *     prm.enter_subsection("Simulation basics");
 * 
 *     prm.declare_entry("Run directory", "./test_me",
 *                       dealii::Patterns::Anything(),
 *                       "Specify a directory where results of "
 *                       "the test are to be stored. This can be either "
 *                       "an absolute path or path relative to the directory "
 *                       "where the program has been started. The default is "
 *                       "subdir called test_me-<date> where <date> will be replaced "
 *                       "by the date at which the program has been started. "
 *                       "this simplifies keeping the projects directory clean "
 *                       "");
 * 
 *     prm.leave_subsection();
 * 
 * 
 *     prm.enter_subsection("CUDA parameters");
 * 
 * 
 *     prm.declare_entry("Device", "0",
 *                       dealii::Patterns::Integer(),
 *                       "which CUDA-enabled GPU should be used");
 * 
 *     prm.declare_entry("Shared Memory", "true",
 *                       dealii::Patterns::Bool(),
 *                       "Whether shared (true) or L1 (false) memory should be used.");
 * 
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Testcase parameters");
 * 
 *     prm.declare_entry("Double-Precision","true",
 *                       dealii::Patterns::Bool(),
 *                       "Decide between double (true) or float (false) precision");
 * 
 *     prm.declare_entry("Matrix size - lower limit", "256",
 *                       dealii::Patterns::Integer(),
 *                       "Start value of the range of matrix sizes tested by the simulation");
 * 
 * 
 *     prm.declare_entry("Matrix size - upper limit", "513",
 *                       dealii::Patterns::Integer(),
 *                       "End value of the range of matrix sizes tested by the simulation");
 * 
 *     prm.declare_entry("Matrix size - step size", "1024",
 *                       dealii::Patterns::Integer(),
 *                       "Increment for the size of the test matrices");
 * 
 *     prm.declare_entry("Average - runs", "10",
 *                       dealii::Patterns::Integer(),
 *                       "Number of runs being used for averaging");
 * 
 *     prm.leave_subsection();
 * 
 * }
 * 
 * 
@endcode
 <a name="plain-Functionget"></a>
@code
 * void
 * step1::SimParams::get(dealii::ParameterHandler & prm)
 * {
 *     prm.enter_subsection("Simulation basics");
 * 
 *     run_dir.setPath(prm.get("Run directory").c_str());
 *     run_dir.makeAbsolute();
 * 
 *     if (!run_dir.exists())
 *         run_dir.mkpath(".");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("CUDA parameters");
 * 
 *     device        = prm.get_integer("Device");
 * 
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Testcase parameters");
 * 
 * 
 *     use_double   = prm.get_bool("Double-Precision");
 * 
 *     matrix_low   = prm.get_integer("Matrix size - lower limit");
 * 
 *     matrix_high  = prm.get_integer("Matrix size - upper limit");
 * 
 *     step_size    = prm.get_integer("Matrix size - step size");
 * 
 *     average_runs = prm.get_integer("Average - runs");
 * 
 * 
 *     prm.leave_subsection();
 * 
 * 
 * }
 * 
 * 
 * #endif
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #include <iostream>
 * #include <vector>
 * 
 * #include <QThread>
 * #include <QTime>
 * 
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/lac/full_matrix.h>
 * 
 * #include <QDir>
 * 
 * 
 * #include <step-1/SimParams.h>
 * 
 * #include <step-1/cuda_driver_step-1.h>
 * #include <step-1/cuda_driver_step-1.hh>
 * 
 * #include <lac/FullMatrixAccessor.h>
 * #include <lac/MatrixCreator.h>
 * 
 * 
 * namespace step1 {
 * 
@endcode
 <a name="plain-ClassCholeskyTest"></a>
@code
 * template<typename Number>
 * class CholeskyTest : public QThread
 * {
 *     static const unsigned int max_display_size = 20;
 * 
 * public:
 * 
 *     CholeskyTest(int n_r,
 *                  dealii::ConvergenceTable& s_table,
 *                  const SimParams &_params);
 * 
 * protected:
 * 
 *     void setup_and_assemble_test_matrix();
 * 
 *     void factorize();
 * 
 *     void check_results();
 * 
 *     void run();
 * 
 * 
 * 
 *     int n_rows;
 * 
 *     dealii::ConvergenceTable & speedup_table;
 * 
 *     dealii::FullMatrix<Number> A, L, L_T;
 * 
 * private:
 * 
 *     const SimParams * params;
 * };
 * 
 * 
@endcode
 <a name="plain-ClassCholesky"></a>
@code
 * class Cholesky {
 * 
 * public:
 *     template<typename T> static void cpu(std::vector<std::vector<T> > & A);
 * 
 *     template<typename T> static void cpu_tiled(T* A, int tile_size);
 * 
 *     template<typename T> static void LLtMult(T * A, const T * L, int n_rows);
 * };
 * 
 * 
@endcode
 <a name="plain-ConstructorCholeskyTest"></a>
@code
 * template<typename Number>
 * CholeskyTest<Number>::CholeskyTest(int n_r,
 *                                    dealii::ConvergenceTable& s_table,
 *                                    const SimParams &_params)
 *     :
 *       n_rows(n_r),
 *       speedup_table(s_table),
 *       params(&_params)
 * {}
 * 
 * 
 * 
@endcode
 <a name="plain-Functionsetup_and_assemble_test_matrix"></a>
@code
 * template<typename Number>
 * void CholeskyTest<Number>::setup_and_assemble_test_matrix()
 * {
 *     this->A.reinit(n_rows, n_rows);
 *     this->L.reinit(n_rows, n_rows);
 *     this->L_T.reinit(n_rows, n_rows);
 * 
 *     QTime t;
 * 
 *     this->speedup_table.add_value("n rows", n_rows);
 * 
 *     std::cout << "Initial Cholesky factor deal.II :" << std::endl;
 *     std::cout << "----------------------------------------------------"
 *               << std::endl;
 * 
 * #ifdef USE_HADAMARD
 *     MatrixCreator::extended_hadamard(n_rows, A);
 * #else
 *     t.start();
 * 
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *         for (unsigned int c = r; c <n_rows; ++c)
 *             L(r,c) = 1e-0*(r+2)*(c+2);
 * 
 * 
 *     if(false)
 *     qDebug("Time for CPU-based setup of Cholesky factor : %f s",
 *            t.elapsed()/1000.);
 * 
 *     if ( L.n_rows() < max_display_size)
 *         L.print(std::cout, 10, 5);
 *     if ( L.n_rows() < max_display_size)
 *         std::cout << std::endl;std::cout << std::endl;
 * 
 *     {
 *         t.restart();
 *         L_T.copy_from(L);
 *         if(false)
 *         qDebug("Time for CPU-based copying of Cholesky factor : %f s",
 *                t.elapsed()/1000.);
 * 
 *         t.restart();
 *         L_T.Tmmult(A, L, false);
 *         if(false)
 *         qDebug("Time for CPU-based multiplication of Cholesky factor"
 *                " : %f s",
 *                t.elapsed()/1000.);
 *         L_T.copy_from(A);
 *     }
 * 
 *     if ( A.n_rows() < max_display_size)
 *     {
 *         std::cout << "Matrix to factorize :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         A.print(std::cout, 10, 5);
 * 
 *         std::cout << std::endl;std::cout << std::endl;
 *     }
 * #endif
 * 
 * }
 * 
@endcode
 <a name="plain-Functionrun"></a>
@code
 * template<typename Number>
 * void CholeskyTest<Number>::run()
 * {
 *     this->setup_and_assemble_test_matrix();
 * 
 *     QTime t;
 *     double cpu_time, gpu_time;
 * 
 *     FullMatrixAccessor<Number> A_h_cpu(A, true);
 *     {
 *         t.restart();
 * 
 * 
 *         Cholesky::cpu_tiled<Number>(A_h_cpu.val(), A_h_cpu.n_rows() );
 * 
 * 
 *         cpu_time =  t.elapsed()/1000.;
 *         if (false)
 *         qDebug("Time for CPU-based Cholesky factorization : %f s",
 *                cpu_time);
 * 
 *         this->speedup_table.add_value("CPU factorization", cpu_time);
 *         this->speedup_table.set_precision("CPU factorization", 10);
 * 
 * 
 *         std::cout << "CPU-factorized Matrix (lower triangle contains "
 *                   << "transposed Cholesky factor) :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         if ( A.n_rows() < max_display_size)
 *             A_h_cpu.print(); //std::cout, 10, 5);
 * 
 *     }
 * 
 *     std::cout << std::endl;std::cout << std::endl;
 * 
 * 
 *     Assert( A.n_rows() == A.n_cols(),
 *             dealii::ExcMessage("Matrix not square! Cholesky impossible"));
 * 
 * 
 *     double kernel_time = 0;
 * 
 *     CUDADriver<Number> run(*params);
 *     t.restart();
 *     {
 *         FullMatrixAccessor<Number> A_h(A, true);
 *         kernel_time = run.factorize(A_h);
 * 
 *         gpu_time =  t.elapsed()/1000.;
 * 
 *         if (false)
 *         qDebug("Time for GPU-based Cholesky factorization %f"
 *                " including data transfer : %f s\n"
 *                "speed up factor factorization : %f netto : %f n_rows : %d\n",
 *                kernel_time,
 *                gpu_time,
 *                cpu_time/kernel_time,
 *                cpu_time/gpu_time,
 *                n_rows);
 * 
 * 
 *         std::cout << "GPU-factorized Matrix (lower triangle contains "
 *                   << "transposed Cholesky factor) :" << std::endl;
 *         std::cout << "----------------------------------------------------"
 *                   << std::endl;
 * 
 *         if ( A_h.n_rows() < max_display_size)
 *             A_h.print(); //std::cout, 10, 5);
 * 
 *         A_h_cpu -= A_h;
 * 
 *         std::cout << "difference of factorized matrices "
 *                         << " :" << std::endl;
 *               std::cout << "----------------------------------------------------"
 *                         << std::endl;
 * 
 *               if ( A_h_cpu.n_rows() < max_display_size)
 *                   A_h_cpu.print(); //std::cout, 10, 5);
 * 
 *               double F_norm =  A_h_cpu.frobenius_norm();
 * 
 *         std::cout << "||A_cpu - A_d ||_F = " << F_norm << "\n"
 *                   << "||A_cpu - A_d ||_F/n_el = " << F_norm/A_h_cpu.n_elements() << "\n"
 *                   << "||A_cpu - A_d ||_F/||A_d||_F = " << F_norm/A_h.frobenius_norm() << std::endl;
 * 
 *     }
 * 
 *     this->speedup_table.add_value("pure GPU fac", kernel_time);
 *     this->speedup_table.set_precision("pure GPU fac", 10);
 * 
 *     this->speedup_table.add_value("GPU fac incl data transfer", gpu_time);
 *     this->speedup_table.set_precision("GPU fac incl data transfer", 10);
 * 
 *     FullMatrixAccessor<Number> A_h(A, true);
 *     FullMatrixAccessor<Number> A_original = A_h;
 * 
 *     {
 *         typename CUDADriver<Number>::TimerName2Value times;
 * 
 *         run.chol_fac(A_h, times);
 * 
 * 
 *         typename CUDADriver<Number>::TimerName2Value::const_iterator
 *                 e=times.begin(),
 *                 end_t=times.end();
 * 
 *         for( ; e != end_t ; ++e)
 *         {
 *             this->speedup_table.add_value(e->first, e->second);
 *             this->speedup_table.set_precision(e->first,10);
 *         }
 *     }
 *     return;
 * }
 * 
 * 
@endcode
 <a name="plain-Functioncpu_tiled"></a>
@code
 * template<typename T>
 * void Cholesky::cpu_tiled(T* A,
 *                          int tile_size)
 * {
 * 
 *     for (int r = 0; r < tile_size; ++r)
 *     {
 *         T sum = 0.;
 *         int idx;
 *         int idx_c;
 * 
 *         for (int u = 0; u < r; ++u)
 *         {
 *             idx = r*tile_size + u;
 *             sum += A[idx] * A[idx];
 *         }
 *         idx = r*tile_size + r;
 *         A[idx] = sqrt(A[idx] - sum);
 * 
 *         for (int c = r+1; c < tile_size; ++c)
 *         {
 *             T tmp = 0.;
 * 
 *             for (int u = 0; u < r; ++u)
 *             {
 *                 idx_c = c*tile_size + u;
 *                 idx   = r*tile_size + u;
 *                 tmp += A[idx_c]*A[idx];
 *             }
 * 
 *             idx_c = c*tile_size + r;
 *             idx   = r*tile_size + c;
 *             A[idx_c]  = A[idx] - tmp;
 *             A[idx_c] /= A[r*tile_size + r];
 *         }
 *     }
 * }
 * 
 * 
@endcode
 <a name="plain-FunctionLLtMult"></a>
@code
 * template<typename T>
 * void Cholesky::LLtMult(T * A, const T * L, int n_rows)
 * {
 *     for (unsigned int r = 0; r < n_rows; ++r)
 *         for (unsigned int c = 0; c <=r; ++c)
 *         {
 *             unsigned int idx = c + (r*(r+1))/2;
 *             unsigned int k_max = std::min(r,c);
 * 
 *             A[idx] = 0.;
 * 
 *             for (unsigned int k = 0; k < k_max; ++k)
 *             {
 *                 unsigned int idx_k   = k + (r*(r+1))/2;
 *                 unsigned int idx_k_T = k + (c*(c+1))/2;
 * 
 *                 A[idx] += L[idx_k]*L[idx_k_T];
 *             }
 *         }
 * }
 * 
 * 
@endcode
 <a name="plain-ClassMyFancySimulation"></a>
@code
 * template<typename Number>
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     MyFancySimulation(SimParams &p);
 * 
 *     void run();
 * 
 *     static std::string precision_id();
 * 
 * private:
 *     const SimParams * params;
 * 
 * };
 * 
 * 
 * 
@endcode
 <a name="plain-Constructor"></a>
@code
 * template <typename Number>
 * step1::MyFancySimulation<Number>::MyFancySimulation(SimParams &p)
 *     :
 *       params(&p)
 * {
 *     cudaSetDevice(params->device); 
 * }
 * 
 * 
 * 
@endcode
 <a name="plain-Functionprecision_id"></a>
@code
 * template<>
 * std::string MyFancySimulation<float>::precision_id()
 * {
 *     return "float";
 * }
 * 
 * template<>
 * std::string MyFancySimulation<double>::precision_id()
 * {
 *     return "double";
 * }
 * 
 * }
 * 
 * 
@endcode
 <a name="plain-Functionrun"></a>
@code
 * template<typename Number>
 * void step1::MyFancySimulation<Number>::run()
 * {   
 * 
 *     std::ostringstream filename;
 * 
 *     filename << "chol_fac_times_" << params->matrix_low << "_" << precision_id().c_str() << ".dat";
 * 
 * 
 *     dealii::ConvergenceTable factorization_times;
 * 
 *     for (int n = params->matrix_low; n < params->matrix_high; n+=params->step_size)
 *     {
 *         CholeskyTest<Number> driver(n, factorization_times, *params);
 * 
 *         driver.start();
 *         / * driver.run();* /
 *         driver.wait();
 * 
 *         std::ofstream out(filename.str().c_str());
 *         factorization_times.write_text(out);
 *     }
 * 
 *     std::cout << "Done." << std::endl;
 * }
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Funktionmain"></a>
@code
 * int main(int argc, char *argv[])
 * {
 *     using namespace step1;
 * 
 *     SimParams params;
 * 
 *     dealii::ParameterHandler prm_handler;
 * 
 *     QDir cwd = QDir::current();
 * 
 *     const QDir launch_dir = cwd;
 *     cwd.setPath("../step-1");
 * 
 *     std::string prm_filename;
 *     if (argc == 1)
 *     {
 *         std::string tmp = argv[0];
 *         int found=tmp.find_last_of('/');
 *         prm_filename = tmp.substr(found+1);
 *         prm_filename += "-Decomp.prm";
 * 
 *         cwd.setPath("./prm");
 *     }
 *     else
 *     {
 *         QFileInfo tmp(argv[1]);
 * 
 *         QString prm_path = tmp.absolutePath();
 *         cwd.setPath(prm_path);
 *         cwd.makeAbsolute();
 *         prm_filename = tmp.fileName().toStdString();
 * 
 *         std::cout << "chosen prm file : " << tmp.absoluteFilePath().toStdString().c_str() << std::endl;
 *     }
 * 
 *     if (!cwd.exists() )
 *         launch_dir.mkpath( cwd.absolutePath() );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 *     SimParams::declare(prm_handler);
 *     prm_handler.read_input (prm_filename);
 * 
 *     QDir::setCurrent(launch_dir.absolutePath());
 * 
 *     params.get(prm_handler);
 * 
 *     cwd.setPath(params.run_dir.absolutePath());
 *     if (!cwd.exists())
 *         cwd.mkpath( "." );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     cwd.setPath("./log");
 *     cwd.makeAbsolute();
 *     if (!cwd.exists())
 *         cwd.mkpath(".");
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(prm_filename.c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 *     QDir::setCurrent(params.run_dir.absolutePath());
 * 
 * 
 *     if (!params.use_double) {
 * 
 *         MyFancySimulation<float> machma_float(params);
 * 
 *         machma_float.run();
 *     }
 *     else {
 * 
 *         MyFancySimulation<double> machma_double(params);
 * 
 *         machma_double.run();
 *     }
 * }
 * 
 * 
 * 
 @endcode
 */
