// @sect3{File: cuda_driver_step-28.h}
#ifndef CUDADriver_STEP_28_H
#define CUDADriver_STEP_28_H
#include <lac/FullMatrixAccessor.h>
#include <lac/blas++.h>
//deal.II includes
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <step-28/cuda_kernel_wrapper_step-28.cu.h>
// use project namespace again:
namespace step28 {
// Some typedefs to make life easier.
typedef typename blas_pp<double, cublas>::FullMatrixAccessor FullMatrixAccessor;
typedef typename blas_pp<double, cublas>::Matrix DeviceMatrix;
typedef typename blas_pp<double, cublas>::Vector DeviceVector;
// @sect4{Class: CUDADriver}
//
// This class manages the communication between host and device.
// It copies the input to the device, runs the required kernels, and copies
// output back to the host for single layer and double layer.
// Together with the BemForm classes it forms a variant of the pointer-to-implementation pattern.
template <int dim>
class CUDADriver {

public:
 void  reinit(const uint DL_SL_matrix_n_rows,
              const uint DL_SL_matrix_n_cols);

    ~CUDADriver ();

    double assemble_bem_matrices(FullMatrixAccessor &DL_matrix_k,
                                 FullMatrixAccessor &SL_matrix_k,
                                 const std::vector<double> &x_k,
                                 const std::vector<double> &q_k,
                                 const std::vector<double> &n_k,
                                 const FullMatrixAccessor &W_k);

    // @sect5{Function: point_to_double_vector}
    //
    // Function to write the contents of a std::vector<Point<dim> > into
    // a std::vector<double> of @p dim times the size.
    // This is done to read these vectors more efficiently on the GPU side.
    // We first copy all x components, then all y components and so on.
void point_to_double_vector(std::vector<double> &out,
                            const std::vector<dealii::Point<dim> > &in)
{
    const uint n_comps2copy = in.size()*dim;

    if (out.size() != n_comps2copy)
        out.resize(n_comps2copy);

    auto itr = in.begin(),
            enditr = in.end();
    uint stride = in.size();

    for(uint c = 0; itr != enditr; ++itr, ++c)
       for (int d = 0; d < dim; d++)
           out[c+d*stride]  = ((*itr)(d));
}



// The arrays for storing the components of the support points, quadrature points and normal vectors
// sorted by coordinate are non-public as they are only needed internally by the implementation
// of the BEM assembly.
protected:

std::vector<double> x_k;
std::vector<double> q_k;
std::vector<double> n_k;



private:
    // Structure to hold the kernels
    BEMKernels<dim> bem_kernels;

    // Input:
    DeviceMatrix d_W;
    DeviceVector d_x;
    DeviceVector d_q;
    DeviceVector d_n;

    // Output:
    DeviceMatrix d_DL_matrix;
    DeviceMatrix d_SL_matrix;


};

} // namespace step28 END

#endif // CUDADriver_STEP_28_H
