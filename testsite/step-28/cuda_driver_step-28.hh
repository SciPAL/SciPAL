// @sect3{File: cuda_driver_step-28.hh}

#ifndef CUDA_DRIVER_STEP_28_HH
#define CUDA_DRIVER_STEP_28_HH

// The declaration of the interface to the CUDA-backend
// is contained in the following header.
#include <step-28/cuda_kernel_wrapper_step-28.cu.h>
#include <step-28/cuda_driver_step-28.h>

// We have to include for CUDA
#include <cuda_runtime_api.h>

// @sect4{Class: CUDADriver}
// @sect5{Constructor: CUDADriver}
//
// The constructor of the CUDADriver class allocates memory for the output
// on the device.
// The parameters are the required sizes.

template<int dim>
void step28::CUDADriver<dim>::reinit(const uint DL_SL_matrix_n_rows,
                                const uint DL_SL_matrix_n_cols)
{


    // Since we use the DeviceMatrix and DeviceVector types, the @p reinit
    // function takes care of the allocation
    this->d_DL_matrix.reinit(DL_SL_matrix_n_rows, DL_SL_matrix_n_cols);
    this->d_SL_matrix.reinit(DL_SL_matrix_n_rows, DL_SL_matrix_n_cols);

   // this->d_DL_matrix.print();
}


// @sect5{Destructor: CUDADriver}
//
// Again, the DeviceMatrix and DeviceVector types take care of deallocation in
// their respective destructors, nothing needs to be done here
template<int dim>
step28::CUDADriver<dim>::~CUDADriver()
{
}

// @sect5{Function: assemble_bem_matrices}
//
// This function takes care of data transfer between host and device, and calls
// the kernel wrapper function.
template<int dim>
double step28::CUDADriver<dim>::assemble_bem_matrices(FullMatrixAccessor &DL_matrix_k,
                             FullMatrixAccessor &SL_matrix_k,
                             const std::vector<double> &x_k,
                             const std::vector<double> &q_k,
                             const std::vector<double> &n_k,
                             const FullMatrixAccessor &W_k)
{

    // Copy input to device. DeviceMatrix and DeviceVector do this in their
    // assignment operators
    this->d_x                       = x_k;
    this->d_q                       = q_k;
    this->d_n                       = n_k;
    this->d_W                       = W_k;

    AssertThrow(q_k.size()%4 == 0,
                dealii::ExcMessage("The number of quadrature points is not a power of two. "
                                   "There exists an implementation for this case, but it "
                                   "is generally slower and therefore disabled."));



    // Run the kernel wrapper function:
    bem_kernels.assemble_bem_matrices(d_DL_matrix,
                                      d_SL_matrix,
                                      d_x,
                                      d_q, d_n, d_W
                /*this->d_DL_matrix.array().val(),
                                      this->d_SL_matrix.array().val(),
                                      this->d_SL_matrix.n_rows(),
                                      this->d_SL_matrix.n_cols(),
                                      this->d_x.array().val(), this->d_x.size()/dim,
                                      this->d_q.array().val(), this->d_q.size()/dim,
                                      this->d_n.array().val(), this->d_n.size()/dim,
                                      this->d_W.array().val(), this->d_W.n_rows(), this->d_W.n_cols()*/);

    // Copy output to host. This, again, this is eased by using DeviceMatrix and DeviceVector
    DL_matrix_k                 = this->d_DL_matrix;
    SL_matrix_k                 = this->d_SL_matrix;

    return 0;
}
#endif // CUDA_DRIVER_STEP_28_HH
