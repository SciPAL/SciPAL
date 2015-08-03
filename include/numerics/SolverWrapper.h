#ifndef SOLVERWRAPPER_H
#define SOLVERWRAPPER_H

//including the headers in order, forces an overwrite of the standard complex format by the fftw format
#include <omp.h>
//STL
#include <complex>
#include <iostream>
//Solvers: Dense cuSolver
#include <cusolverDn.h>
//deal.II
#include <deal.II/base/parameter_handler.h>
//SciPAL
#include <base/VTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/Shape.h>
//Qt Includes
#include <QDir>
#include <QString>

namespace SciPAL {

// This is the template prototype which will be used later (for different types)
//template <typename T, typename BW> long
//SVD(const SciPAL::Matrix<T, BW>& A,
//    SciPAL::Matrix<T, BW>& U,
//    SciPAL::Vector<typename SciPAL::VTraits<T, BW::arch>::NumberType, BW>& S,
//    SciPAL::Matrix<T, BW>& Vt);

// This is the q&d implementation of a wrapper for the SVD in cusolver
long SVD(const SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas>& A,
         SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas>& U,
         SciPAL::Vector<double, cublas>& S,
         SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas>& Vt)
{
    unsigned int m = A.n_rows(), n = A.n_cols();

    if (m < n)
    {
        // TODO: Needs to be checked!
        SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> A_tmp(A.n_cols(), A.n_rows());      // create temporary A_tmp with transposed size
        A_tmp = SciPAL::adjoint<SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> >(A);       // copy content of A adjoint into A_tmp

        SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> U_tmp, Vt_tmp;
        long returnValue = SVD(A_tmp, U_tmp, S, Vt_tmp);                                        // performe SVD on A_tmp

        U.reinit(Vt_tmp.n_cols(), S.n_rows());
        Vt.reinit(S.n_rows(), U_tmp.n_rows());

        U = SciPAL::adjoint<SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> >(Vt_tmp);
        Vt = SciPAL::adjoint<SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> >(U_tmp);

        return returnValue;
    }
    U.reinit(m, n);
    SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> U_tmp(m, m);
    Vt.reinit(n, n);
    S.reinit(n);

    // --- CUDA solver initialization
    cusolverDnHandle_t solver_handle;
    cusolverDnCreate(&solver_handle);

    cusolverStatus_t stat;

    // --- device side SVD workspace and matrices
    int work_size = 0;
    int *devInfo;       cudaMalloc((void**)(&devInfo),          sizeof(int));

    stat = cusolverDnSgesvd_bufferSize(solver_handle, m, n, &work_size);
    if(stat != CUSOLVER_STATUS_SUCCESS ) std::cout << "Initialization of cuSolver failed. \n";

    SciPAL::Vector<SciPAL::CudaComplex<double>, cublas> work(work_size);

    // --- CUDA SVD execution
    stat = cusolverDnZgesvd(solver_handle, 'A', 'A', m, n, A.data_ptr, m, S.data_ptr, U_tmp.data_ptr, m, Vt.data_ptr, n, work.data_ptr, work_size, NULL, devInfo);

    cudaDeviceSynchronize();

    cusolverDnDestroy(solver_handle);

    /* Until the views are working...*/
    SciPAL::Matrix<SciPAL::CudaComplex<double>, cublas> ShrinkMatrix(U_tmp.n_cols(), S.n_rows());
    for (unsigned int i = 0; i < S.n_rows(); i++)
        ShrinkMatrix(i, i, 1.0);
    U = U_tmp * ShrinkMatrix;

//    std::cout << "Direct inside the SVD " << stat << std::endl;
//    std::cout << "U" << std::endl;
//    U.print();
//    std::cout << "S" << std::endl;
//    S.print();
//    std::cout << "Vt" << std::endl;
//    Vt.print();

    return stat;
}


} //namespace SciPAL end
#endif // SOLVERWRAPPER_H

