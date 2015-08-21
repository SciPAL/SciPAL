#ifndef CUDADriver_STEP_4_H
#define CUDADriver_STEP_4_H

// deal.II classes
#include <deal.II/lac/full_matrix.h>

// lab course related headers.
#include <lac/blas++.h>
#include <lac/cublas_algorithms.h>

// The declaration of the interface to the CUDA-backend
// are all contained in the following header.
#include <cuda_kernel_wrapper_step-4.cu.h>

namespace step4 {

    // @sect3{Class: CudaQRDecomposition}
    //
    // QR factorization using Householder transformations.
    // After the factorization @p Q_d and @p R_d contain the results.
    // @param T : the number type. Althoug it could be any of real
    // or complex floats or doubles we assume that only the real
    // versions are used.
    // @param blas : must be one of the BLAS wrapper types.
template <typename T, typename blas>
class CudaQRDecomposition {

public:
        // Some abbreviations
    typedef blas_pp<T, blas> BLAS;

    typedef typename BLAS::blas_wrapper_type  BW;
    typedef typename BLAS::FullMatrixAccessor FullMatrixAccessor;
    typedef typename BLAS::Matrix             Matrix;
    typedef typename BLAS::SubMatrix          SubMatrix;
    typedef typename BLAS::MatrixSubCol       MatrixSubCol;
    typedef typename BLAS::Vector             Vector;
    typedef typename BLAS::SubColVector       SubColVector;
    typedef typename BLAS::SubVectorBase      SubVectorBase;

    CudaQRDecomposition();

    ~CudaQRDecomposition ();


    T householder (const dealii::FullMatrix<T>& A);

    T lapack_based(const dealii::FullMatrix<T> &A);

    const Matrix & Q() const { return Q_d; }
    const Matrix & R() const { return R_d; }

private:
    Matrix Q_d, R_d;

    void check_matrix_dimensions(const dealii::FullMatrix<T>& A);
};

} // namespace step4 END

#endif // CUDADriver_STEP_4_H
