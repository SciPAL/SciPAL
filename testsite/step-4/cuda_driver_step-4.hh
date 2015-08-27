#ifndef CUDA_DRIVER_STEP_4_HH
#define CUDA_DRIVER_STEP_4_HH

#include <step-4/cuda_driver_step-4.h>

// Header containing the declarations of the wrapper
// functions - not the kernels.
// This header is the interface between the CUDA part compiled by nvcc and
// the host part compiled by gcc.
#include <step-4/cuda_kernel_wrapper_step-4.cu.h>

//#include <cutil_inline.h>


#include <helper_cuda.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include <iomanip>
#include <cmath>
#include <QtGlobal>

#include <base/sign.h>
#include <base/CUDATimer.h>

#include <lac/expression_templates_host.h>

#include /*<SciPal/*/  <lac/FullMatrixAccessor.h>

#include <lac/cublas_algorithms.h>


//! @sect4{Macro: print_intermediate }
//! Helper function for debugging, prints intermediate results.
//! since there is no nice way to conditionally define macros, we have to
//! comment and uncommet th definition by hand :(

#define DEBUG_RESULTS(var)\
//    std::cout << #var <<"("<<c<<"):\n"; var.print();


// @sect3{Implementation: CudaQRDecomposition Methods}
//
// @sect4{Constructor: CudaQRDecomposition}
//
// The CTor of the driver class initializes the BLAS library.
template <typename T, typename blas>
step4::CudaQRDecomposition<T, blas>::CudaQRDecomposition()
{
    BW::Init();
}



// @sect4{Destructor: CudaQRDecomposition}
//
// Shutdown the BLAS library.
template <typename T, typename blas>
step4::CudaQRDecomposition<T, blas>::~CudaQRDecomposition()
{
    BW::Shutdown();
}



// @sect4{Function: check_matrix_dimensions}
//
// Our QR factorization requires that the matrix which is to be
// factorized has more rows than columns.
// This is checked in this function and in case this requirement
// is not an exception is thrown. In debug mode the program aborts whereas
// in release mode you could catch the exception and try
// to somehow save the computation.
// @param A : Matrix to check.
template <typename T, typename blas>
void
step4::CudaQRDecomposition<T, blas>::check_matrix_dimensions(const dealii::FullMatrix<T>& A)
{
    const bool is_matrix_dim_correct = A.n_rows() >= A.n_cols();

    #ifdef QT_NO_DEBUG
    AssertThrow(is_matrix_dim_correct,
            dealii::ExcMessage("n_rows must be greater n_cols") );
    #else
    Assert(is_matrix_dim_correct,
       dealii::ExcMessage("n_rows must be greater n_cols") );
    #endif
}



// @sect4{Function: householder}
//
// As above we have to take the transpose copy of the input matrix.
// By having this function accepting a deal.II FullMatrix we enforce
// this conversion so that on both sides the storage format for the matrix
// which is to be factorized is unambiguous.
// Outside of this function it is row-major and inside column-major.
// The function returns the time needed for the factorization excluding the
// time needed for memory transfers.
// @param A : matrix to factorize.
template <typename T, typename blas>
T
step4::CudaQRDecomposition<T, blas>::householder(const dealii::FullMatrix<T> &A)
{
    // Check, whether there are more rows than columns.
    // Otherwise the factorization cannot be performed.
    check_matrix_dimensions(A);

        // Local variables for some frequently used properties
        // of @p A.
    int n_rows = A.n_rows();
    int n_cols = A.n_cols();
        // Unnecessary copies are ugly hence the name.
    ::FullMatrixAccessor<T> A_t_ugly(A,
                                   true /* transpose copy*/);

        // Reinitialize the device representation @p Q_d, @p R_d
        // of $Q$ and $R$.
        // During the factorization $A$ is overwritten with $R$ thus
        // we put the original $A$ into @p R_d right away.
        // Now copy from host to device.
    this->R_d = A_t_ugly;

        // Identity matrix for image space. This initializes $Q$.
    dealii::IdentityMatrix eye_image(n_rows);
    Q_d = eye_image;

    T alpha = 0, sigma = 0;

    Vector x(n_rows);
    Vector k(n_rows);
    Vector u(n_rows);

    SubMatrix A_rest(R_d, 0, 0);

    SubMatrix Q_0c(Q_d, 0, 0);

    CUDATimer timer;

    // Eliminate column-wise all entries below the diagonal by
    // Householder reflections.
    for(int c=0; c<n_cols-1; c++) {

        // Move the upper left corner of the view
        // to the next diagonal element.
        A_rest.reset(c, n_rows, c, n_cols);

        DEBUG_RESULTS(A_rest);

        // For $Q$ we have to mask all previously processed columns.
        Q_0c.reset(0, n_rows, c, n_rows);

        // Store the current diagonal element.
        T A_tt1 = R_d(c, c);

        // The non-zero part of the column Householder vector is initialized
        // with the subcolumn of @p A.
        MatrixSubCol sub_col_A(R_d, c, c);
        DEBUG_RESULTS(sub_col_A);

        SubColVector view_u(u, c, 0);
        SubColVector view_x(x, c, 0);
        SubColVector view_k(k, 0, 0);


        // Compute the i-th Householder vektor $\hat u^{(i)}$
        {
            // $ \alpha = sign (A_{tt})\|y\|_2$
            alpha    =  sign(A_tt1) * sub_col_A.l2_norm();
            T beta_t = A_tt1 + alpha;

            // If the following condition holds, the subcolumn is already
            // zero and there is nothing to do.
            if (std::fabs(alpha) < 1e-14)
                continue;

            sigma    = beta_t/alpha;

            view_u = sub_col_A;   // FIX ME

            u.add(c, alpha);

            Assert(std::fabs(beta_t) > 1e-14, dealii::ExcMessage("Div by 0!") );

            u *= 1./beta_t;
            DEBUG_RESULTS(u)

        }

        // The matrix-matrix product $(I - \sigma u u^T)A$
        // for updated the remaing part of @p A is replaced by
        // \f{eqnarray*}
        // A & = & A - \sigma u u^T A \\
        // & = & A -  u (\sigma A^Tu)^T \,. \\
        //\f}
        // To do this, we define $x := A^Tu$
        view_x = transpose(A_rest) * view_u; //FIX ME
        DEBUG_RESULTS(view_x);
        DEBUG_RESULTS(x)

        // and take into account $\sigma$ in the outer product
        // which due to the use of SubMatrixViews is executed only on
        // the truly affected part.
        // \f{eqnarray*} A(k:m,k:n) &=& A(k:m,k:n) - \sigma uu^T A(k:m,k:n) \\
        //                          &=& A(k:m,k:n) - \sigma ux^T\,,
        //\f}
        //FIX ME use expressions
        A_rest.add_scaled_outer_product(-sigma, view_u, view_x);
        DEBUG_RESULTS(A_rest);

        // Updating $Q$ is similar
        // \f{eqnarray*} Q(1:m,k:n) &=& Q(1:m,k:n) - \sigma Q(1:m,k:n)uu^T \\
        // k & := &   Q(1:m,k:n)u \\
        // \Rightarrow   Q(1:m,k:n) &=& Q(1:m,k:n) - \sigma ku^T
        //\f}


        view_k = Q_0c * view_u;
        DEBUG_RESULTS(k);

        Q_0c.add_scaled_outer_product(-sigma, view_k, view_u);
        DEBUG_RESULTS(Q_0c);
        DEBUG_RESULTS(Q_d);
        DEBUG_RESULTS(R_d);

    }

    timer.print_elapsed("Time spent on Householder-QR excluding host-device memory traffic : ");

    return timer.elapsed();
}





// @sect4{Function: lapack_based}
//
// To have some competition we also implement a wrapper fucntion for calling the
// QR-factorization provided by LAPACK.
template <typename T, typename blas>
T
step4::CudaQRDecomposition<T, blas>::lapack_based(const dealii::FullMatrix<T> &A)
{
#ifdef USE_LAPACK
    FullMatrixAccessor A_t_ugly(A,
                                true /*transpose while copying*/);

    __CLPK_integer m = A_t_ugly.n_rows();
    __CLPK_integer n = A_t_ugly.n_cols();

    __CLPK_integer lda = A_t_ugly.n_rows();

    __CLPK_integer lwork = -1;

    __CLPK_integer info;

    std::vector<float> tau(std::min(A_t_ugly.n_rows(),
                                    A_t_ugly.n_cols()));

    std::vector<float> work(1);

    float * A_val = const_cast<float*>(A_t_ugly.val());

    sgeqrf_(&m, &n, A_val, &lda, &tau[0], &work[0], &lwork, &info);


    if (info != 0) {
        printf("sgeqrf failed with error code %d\n", (int)info);
        return 0;
    }

    lwork = work[0];

    work.resize(lwork);

    sgeqrf_(&m, &n, A_val, &lda, &tau[0], &work[0], &lwork, &info);

    if (info != 0) {
        printf("sgeqrf failed with error code %d\n", (int)info);
        return 0;
    }
    this->R_d = A_t_ugly;
#endif
}



#endif // CUDA_DRIVER_STEP_4_HH
