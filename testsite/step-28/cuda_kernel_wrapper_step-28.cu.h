// @sect3{File: cuda_kernel_wrapper_step-28.cu.h}
#ifndef CUDA_KERNEL_STEP_28_CU_H
#define CUDA_KERNEL_STEP_28_CU_H

#include <cuda.h>
#include <lac/release/ShapeData.h>

namespace step28 {


// @sect4{Struct: BEMKernels}
// A struct to hold the kernel wrappers. In this project, this is only a
// single wrapper, since only one part is parallelized using CUDA.
template<int dim>
struct BEMKernels {

    void assemble_bem_matrices(SciPAL::ShapeData<double> DL_matrix,
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
                               const double * W, const int W_n_rows, const int W_n_cols*/);

};


}

#endif // CUDA_KERNEL_STEP_28_CU_H
