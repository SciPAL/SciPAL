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

Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
*/


#ifndef CUDA_DRIVER_STEP_42_HH
#define CUDA_DRIVER_STEP_42_HH


// The declaration of the interface to the CUDA-backend
// is contained in the following header.
#include <cuda_driver_step-42.h>
#include <cuda_kernel_wrapper_step-42.cu.h>



#include <lac/development/cublas_Matrix.h>


#include <lac/development/cublas_Vector.h>
#include <lac/blas++.h>
#include <base/CudaComplex.h>

// We have to include
#include <cuda_runtime_api.h>


        // @sect4{Constructor: CUDADriver}
        //
        // The constructor of the driver class allocates memory
        // on the GPU and copies data from host to device.
        // Furthermore, it keeps a pointer to the original host data for the
        // case that it has to be modified by the GPU results.
        // @param v_h : Pointer to a linear array in host-side memory that is to be copied to the GPU.
        // @param n : Number of entries of @p v_h.
step42::CUDADriver::CUDADriver() {

    BW::Init();
}




        // @sect4{Function: gemm_tests}
        //
        // This function tests the various special cases contained
        // in BLAS' gemm function.
void step42::CUDADriver::gemm_tests()
{

    // some dummy vectors
    const unsigned int n_rows = 4;
    const unsigned int n_cols = 3;
    const unsigned int n_elements = n_rows * n_cols;

    std::vector<Number>
            a(n_elements, 1.),
            b(n_elements, 2.),
            c(n_rows * n_rows, 1.23);

    for (unsigned int i = 0; i < b.size(); i++ )
        b[i] = i+1;


    SciPAL::Matrix<Number, BW>
            A(n_rows, n_cols, a),
            B(n_cols, n_rows, b),
            C(n_rows, n_rows, c);

     Number alpha = 1.1;
     Number beta = 2.;

     std::cout << "A : " << std::endl;
     A.print();

     std::cout << "B : " << std::endl;
     B.print();

     std::cout << "C : " << std::endl;
     C.print();



     std::cout << " ============ C = " << alpha << " * A ======" << std::endl;
     C = alpha * A; // * B + beta * C;
       std::cout << "C : " << std::endl;
     C.print();

     std::cout << " ============ C = A * B ======" << std::endl;
     C = A * B;   std::cout << "C : " << std::endl; C.print();

     std::cout << " ============ C = B * A ======" << std::endl;
     C = B * A;   std::cout << "C : " << std::endl; C.print();

     std::cout << " ============ C = alpha * A * B ======" << std::endl;
     C = alpha * A * B; std::cout << "C : " << std::endl; C.print();


     //sMMaM test
     std::cout << " ============ C = " << alpha << " * A * B + "<<"C======" << std::endl;
     c.clear();
     c.resize(n_rows * n_rows, 1.);
     SciPAL::Matrix<Number, BW>
             D(n_rows, n_rows, c);



     C = D;
     std::cout << "C : " << std::endl; C.print();
     C = alpha * A * B + // beta *
             C; std::cout << "C : " << std::endl; C.print();


    //gemm test
     std::cout << " ============ C = " << alpha << " * A * B + "  << beta  << " * C======" << std::endl;
     C = D;
     std::cout << "C : " << std::endl; C.print();

     C = alpha * A * B + beta * C;
     std::cout << "C : " << std::endl; C.print();
}


// @sect4{Function: gemv_tests}
//
// This function tests the various special cases contained
// in BLAS' gemv function.
//
// Currently, it is rather a test for the vector arithmetic.
// test vector expressions
void step42::CUDADriver::gemv_tests()
{
#ifndef nUSE_ARRAY_EXPRESSIONS
     const unsigned int n_rows = 4;
     const unsigned int n_cols = 4;
     const unsigned int n_elements = n_rows * n_cols;

     Number alpha = 1.1;
     Number beta = 2.;

     std::vector<Number>
                a(n_elements, 1.),
                b(n_elements, 2.);


     for (unsigned int i = 0; i < a.size(); i++ )
         a[i] = i+1;

     SciPAL::Vector<Number, BW> vA, vB(n_elements), vC;
       vA = a;
       // This sets all elements of vB to 2.3, note: vector needs to be initialized.
       vB = SciPAL::Literal<Number>(2.3);
       vC = a;


       std::cout << "vA : " << std::endl;
       vA.print();

       std::cout << "vB : " << std::endl;
       vB.print();

       std::cout << " ============ vC = " << alpha << " * vA ======" << std::endl;
       vC = alpha * vA;
         std::cout << "vC : " << std::endl;
       vC.print();

       std::cout << " ============ vC = " << alpha << " * vA + vB ======" << std::endl;
       vC = alpha * vA + vB;
         std::cout << "vC : " << std::endl;
       vC.print();

       std::cout << " ============ vA = sin(vC) ======" << std::endl;
       const unsigned int n_sin_elements = n_elements;
       std::vector<Number> d(n_sin_elements);
       for(uint i = 0; i < d.size(); i++)
           d[i] = i* 2.* M_PI / d.size();

       SciPAL::Vector<Number, BW> vD; vD = d; //(n_sin_elements, 1, d);
      vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
         std::cout << "sin(vD) : " << std::endl;
       vD.print();
       //vC = alpha * sin(vA) + vB;


       // After the element-wise sine of a vector we do the same for a matrix.
       SciPAL::Matrix<Number, BW>
               A(n_rows, n_cols, d);
       A = sin(A);
       A.print();


       std::cout << " ============ linear combination test ======" << std::endl;
       vC = 2.0 * vA;
       std::cout << "vC = 2.0 * vA" << std::endl;
       vC.print();


       vC = vA + vB;
       std::cout << "vC = vA + vB" << std::endl;
       vC.print();


       vC = vA + 2.0 * vB;
       std::cout << "vC = vA + 2.0 * vB" << std::endl;
       vC.print();

       vC = 2.0 * vA + 3.0 * vB;
       std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
       vC.print();

//       vC = 2.0 * vA + 3.0 * vB + 4.0 * vD;
//       std::cout << "vC = 2.0 * vA + 3.0 * vB + 4.0 * vD" << std::endl;
//       vC.print();

       //combined expr test
       vC =sin(2.0 * vA + 3.0 * vB);
       std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
       vC.print();

       //test pointwise sqrt
       vC = sqrt(vC);
       std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
       vC.print();

       //test pointwise *
       vC = vA && vB;
       std::cout << "vC = vA .* vB" << std::endl;
       vC.print();

       //test pointwise /
       vC = vA || vB;
       std::cout << "vC = vA ./ vB" << std::endl;
       vC.print();

       //combined expr test
       vC = (vA + vB) || (vA - vB);
       std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
       vC.print();

       //combined expr test
//       vC =abs(cos(2.0 * vA + 3.0 * vB));
//       std::cout << "vC = abs(cos(2.0 * vA + 3.0 * vB)" << std::endl;
//       vC.print();

#endif
}

void step42::CUDADriver::complex_tests()
{
    std::cout<<"Entering tests for complex number array exptessions."<<std::endl;
#ifndef nUSE_ARRAY_EXPRESSIONS
     const unsigned int n_rows = 4;
     const unsigned int n_cols = 4;
     const unsigned int n_elements = n_rows * n_cols;

     SciPAL::CudaComplex<Number> alpha(1.1);
     SciPAL::CudaComplex<Number> beta (2.);

     std::vector<std::complex<Number> >
                a(n_elements, 1.),
                b(n_elements, std::complex<Number>(2., 2.0));


     for (unsigned int i = 0; i < a.size(); i++ )
         a[i] = std::complex<Number>(i+1, (i+1)/2.); //generate some inputs

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vA, vB, vC;
       vA = a;
       vB = b;
       vC = a;
//               vA(n_elements, 1, a),
//               vB(n_elements, 1, b),
//               vC(n_elements, 1, a);

       std::cout << "vA : " << std::endl;
       vA.print();

       std::cout << "vB : " << std::endl;
       vB.print();

       std::cout << " ============ vC = " << alpha.real() << " * vA ======" << std::endl;
       vC = alpha * vA;
         std::cout << "vC : " << std::endl;
       vC.print();

       std::cout << " ============ vC = " << alpha.real() << " * vA + vB ======" << std::endl;
       vC = alpha * vA + vB;
         std::cout << "vC : " << std::endl;
       vC.print();

       std::cout << " ============ vA = sin(vC) ======" << std::endl;
       const unsigned int n_sin_elements = n_elements;
       std::vector<std::complex<Number> > d(n_sin_elements);
       for(uint i = 0; i < d.size(); i++)
           d[i] = std::complex<Number>(i* 2.* M_PI / d.size(), i* 4.* M_PI / d.size()) ;

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vD; vD = d; //(n_sin_elements, 1, d);
      vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
         std::cout << "sin(vD) : " << std::endl;
       vD.print();
       //vC = alpha * sin(vA) + vB;

       std::cout << " ============ Matrix A = sin(A) ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
               A(n_rows, n_cols, d);
       A = sin(A);
       A.print();

       std::cout << " ============ Matrix B = sqrt(A) ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
               B(A);
       B = sqrt(A);
       B.print();

       std::cout << " ============ Matrix A = .exp(B) ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
//       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
//               A(B);
       A = exp(B);
       A.print();

       std::cout << " ============ Matrix A = B*C ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
               C(A);
       A = B*C;
       A.print();


       std::cout << " ============ linear combination test ======" << std::endl;

       //This does not work : vC = 2.0 * vA; because of mismatching type.
       // Thus we hav to use numbers wrapped in CudaComplexes
       vC = beta * vA; ;
       std::cout << "vC = 2.0 * vA" << std::endl;
       vC.print();


       vC = vA + vB;
       std::cout << "vC = vA + vB" << std::endl;
       vC.print();


       vC = vA + beta * vB;
       std::cout << "vC = vA + 2.0 * vB" << std::endl;
       vC.print();

       vC = cplxNumber(2.0) * vA + cplxNumber(3.0) * vB;
       std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
       vC.print();

       //combined expr test
       vC =sin(cplxNumber(2.0) * vA + cplxNumber(3.0) * vB);
       std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
       vC.print();

       //test pointwise sqrt
       vC = sqrt(vC);
       std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
       vC.print();

       //test pointwise *
       vC = vA && vB;
       std::cout << "vC = vA .* vB" << std::endl;
       vC.print();

       //test pointwise /
       vC = vA || vB;
       std::cout << "vC = vA ./ vB" << std::endl;
       vC.print();

       //combined expr test
       vC = (vA + vB) || (vA - vB);
       std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
       vC.print();

       //combined expr test
//       vC =abs(cos(2.0 * vA + 3.0 * vB));
//       std::cout << "vC = abs(cos(2.0 * vA + 3.0 * vB)" << std::endl;
//       vC.print();

#endif
}

void step42::CUDADriver::cusolver_demonstration()
{
    // some dummy vectors
    const unsigned int n_rows = 3;
    const unsigned int n_cols = 3;
    const unsigned int n_elements = n_rows * n_cols;

    std::vector<cplxNumber> a(n_elements, 1.);

    for (unsigned int i = 0; i < a.size(); i++)
        a[i] = i+1;

    SciPAL::Matrix<cplxNumber, BW>
            A(n_rows, n_cols, a),
            U(n_rows, n_cols),
            Vt(n_cols, n_cols);

    SciPAL::Vector<Number, BW> S(n_cols);
    SciPAL::Vector<cplxNumber, BW> S_cplx(n_cols);

    A.print();

    int returnValue = SciPAL::SVD(A, U, S, Vt);

    std::cout << returnValue << std::endl << "Singular values: " << std::endl;

    S.print();

    std::cout << "Original matrix:" << std::endl;

    for (unsigned int i = 0; i < n_cols; i++)
        S_cplx.set(i, (cplxNumber)(S(i)));
    SciPAL::Matrix<cplxNumber, BW> A_tmp(n_rows,n_cols);
    A_tmp = U * SciPAL::diag<SciPAL::Vector<cplxNumber, BW> >(S_cplx);
    A = A_tmp * Vt;

    A.print();
}

void step42::CUDADriver::feature_demonstration()
{
    Number * bla = new Number[3];
    Number * bla2 = new Number[3];
    Number * bla3 = new Number[3];
    std::vector<Number*> h_testV(5);
    h_testV[0] = bla; std::cout<<"ptr1 " << bla << std::endl;
    h_testV[1] = bla2;std::cout<<"ptr2 " << bla2 << std::endl;
    h_testV[2] = bla3;std::cout<<"ptr3 " << bla3 << std::endl;

    SciPAL::Vector<Number*, cublas> d_testV(5);

 d_testV = h_testV;
 d_testV.print();




}

#endif // CUDA_DRIVER_STEP_42_HH
