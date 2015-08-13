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

#include<numeric>
// The declaration of the interface to the CUDA-backend
// is contained in the following header.
#include <cuda_driver_step-42.h>
#include <cuda_kernel_wrapper_step-42.cu.h>



#include <lac/cublas_Matrix.h>


#include <lac/cublas_Vector.h>
#include <lac/blas++.h>
#include <base/CudaComplex.h>

// We have to include
#include <cuda_runtime_api.h>

#include <lac/stack.h>

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

     //MMaM test
      std::cout << " ============ C = A * B + C======" << std::endl;
      C = D;
      std::cout << "C : " << std::endl; C.print();

      C =  A * B + C;
      std::cout << "C : " << std::endl; C.print();

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

       std::cout << " ============ vD = sin(vD) ======" << std::endl;
       const unsigned int n_sin_elements = n_elements;
       std::vector<Number> d(n_sin_elements);
       for(uint i = 0; i < d.size(); i++)
           d[i] = i* 2.* M_PI / d.size();

       SciPAL::Vector<Number, BW> vD; vD = d; //(n_sin_elements, 1, d);
      vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
         std::cout << "sin(vD) : " << std::endl;
       vD.print();
       //vC = alpha * sin(vA) + vB;


       std::cout << " ============ linear combination test ======" << std::endl;
       vC = 2.0 * vA;
       std::cout << "vC = 2.0 * vA" << std::endl;
       vC.print();

       std::cout << "vA : " << std::endl;
       vA.print();

       std::cout << "vB : " << std::endl;
       vB.print();
       vC = vA + vB;
       std::cout << "vC = vA + vB" << std::endl;
       vC.print();
       std::cout << "vA : " << std::endl;
       vA.print();

       std::cout << "vB : " << std::endl;
       vB.print();

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

       std::cout << " ============ pointwise arithmetic with matrices======" << std::endl;
{       // After the element-wise sine of a vector we do the same for a matrix.

       SciPAL::Matrix<Number, BW> A(n_rows, n_cols, d);
       A = sin(A);
       SciPAL::Matrix<Number, BW> B(n_rows, n_cols, b);
       SciPAL::Matrix<Number, BW> C(n_rows, n_cols);

       std::cout << " A =" << std::endl;

       A.print();
       std::cout << " B =" << std::endl;
       B.print();
       std::cout << " C =" << std::endl;
       C.print();

       C = abs(A);
       std::cout << "C = abs(A)" << std::endl;
       C.print();

       //test pointwise * on same matrix
//       C = C && C;
       std::cout << "C = C .* C" << std::endl;
       C.print();

       //test pointwise *
//       C = A && B;
       std::cout << "C = A .* B" << std::endl;
       C.print();

       //test pointwise /
       C = A || B;
       std::cout << "C = A ./ B" << std::endl;
       C.print();

       //concatuated expr test
       C = (A + B) || (A - B);
       std::cout << "C = (A + B) || (A - B)" << std::endl;
       C.print();

       //combined bin-un-expr test
       C = abs(A);
       C = (A + B) || C;
       //TODO not working:
//       C = (A + B) || abs(A);
       std::cout << "C =  (A + B) || abs(A)" << std::endl;
       C.print();

}
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

     std::vector<SciPAL::CudaComplex<Number> >
                a(n_elements, 1.),
                b(n_elements, SciPAL::CudaComplex<Number>(2., 2.0)),
                c(n_elements);
     std::iota(c.begin(), c.end(), 1);


     for (unsigned int i = 0; i < a.size(); i++ )
         a[i] = std::complex<Number>(i+1, (i+1)/2.); //generate some inputs

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> //vA, vB, vC;
//       vA = a;
//       vB = b;
//       vC = a;
               vA( a),
               vB( b),
               vC( a);

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
        std::vector<SciPAL::CudaComplex<Number> > d(n_sin_elements);
       for(uint i = 0; i < d.size(); i++)
           d[i] = SciPAL::CudaComplex<Number>(i* 2.* M_PI / d.size(), i* 4.* M_PI / d.size()) ;

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vD; vD = d; //(n_sin_elements, 1, d);
      vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
         std::cout << "sin(vD) : " << std::endl;
       vD.print();
       //vC = alpha * sin(vA) + vB;

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

//       //test pointwise *
//       vC = vA && vB;
//       std::cout << "vC = vA .* vB" << std::endl;
//       vC.print();

//       //test pointwise /
//       vC = vA || vB;
//       std::cout << "vC = vA ./ vB" << std::endl;
//       vC.print();

//       //combined expr test
//       vC = (vA + vB) || (vA - vB);
//       std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
//       vC.print();

       //combined expr test
//       vC =abs(cos(2.0 * vA + 3.0 * vB));
//       std::cout << "vC = abs(cos(2.0 * vA + 3.0 * vB)" << std::endl;
//       vC.print();

       std::cout << " ============ matrix test ======" << std::endl;
{       std::cout << " ============ Matrix A = sin(A) ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW> A(n_rows, n_cols, c);
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW> B(n_rows, n_cols, b);
       SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW> C(n_rows, n_cols);

       std::vector<Number> tmp(n_elements);
       std::iota(tmp.begin(), tmp.end(), 1);
       SciPAL::Matrix<Number, BW> D(n_rows, n_cols, tmp);
       SciPAL::Matrix<Number, BW> E(n_rows, n_cols, tmp);


       std::cout << " A =" << std::endl;

       A.print();
       std::cout << " B =" << std::endl;
       B.print();
       std::cout << " C =" << std::endl;
       C.print();
       std::cout << " D =" << std::endl;
       D.print();
       std::cout << " E =" << std::endl;
       E.print();

       C = abs(A);
       std::cout << "C = abs(A)" << std::endl;
       C.print();

       //test pointwise * on same matrix
       C = C && C;
       std::cout << "C = C .* C" << std::endl;
       C.print();

       //test pointwise *
       C = A && B;
       std::cout << "C = A .* B" << std::endl;
       C.print();

       //test && for mixed numbertypes
//       C = B && D;
       std::cout << "totally wrong: mixed numbertypes: C = A && D;" << std::endl;
       C.print();


       //test pointwise /
       C = A || B;
       std::cout << "C = A ./ B" << std::endl;
       C.print();

       //concatuated expr test
       C = (A + B) || (A - B);
       std::cout << "C = (A + B) || (A - B)" << std::endl;
       C.print();

       //combined bin-un-expr test
       D = abs(A);
//       C = (B || D) && D;
       std::cout << "C = (B || abs(A)) && abs(A);" << std::endl;
       C.print();

       std::cout << "D = abs(C)" << std::endl;
       D = abs(C);
       D.print();

       //combined bin-un-expr test
       D = abs(B);
//       B = (B || D) && E ;
       std::cout << "fails horribly: B = (B || abs(B)) && E;" << std::endl;
       B.print();

       std::cout << "D = abs(B)" << std::endl;
       D = abs(B);
       D.print();

       std::cout << " ============ Matrix B = sqrt(A) ======" << std::endl;
       // After the element-wise sine of a vector we do the same for a matrix.
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

       A = B*C;
       A.print();
}

#endif
}


void step42::CUDADriver::feature_demonstration()
{
    std::cout<<"=======Misceleaneous tests\n Test Vector of Pointers======="<<std::endl;
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


 std::cout<<"Test arbitrary copy directions"<<std::endl;
 SciPAL::Matrix<Number, blas> M1(2,2);
 M1(0,0,1.); M1(0,1,2.);
 M1(1,0,3.); M1(1,1,4.);

 std::cout<<"Matrix on host\n M1:"<<std::endl;
M1.print();

 SciPAL::Matrix<Number, blas> M2(2,2);
 std::cout<<"Other Matrix on host\n M2=M1:"<<std::endl;
 M2 = M1;
 M2.print();

 SciPAL::Matrix<Number, cublas> M3(2,2);
 std::cout<<" Matrix on GPU\n M3=M2:"<<std::endl;
 M3 = M2;
 M3.print();

 SciPAL::Matrix<Number, cublas> M4(2,2);
 std::cout<<"Other Matrix on GPU\n M4=M3:"<<std::endl;
 M4 = M3;
 M4.print();

 SciPAL::Matrix<Number, blas> M5(2,2);
 std::cout<<"back on host\n M5=M4:"<<std::endl;
 M5 = M4;
 M5.print();





}

void step42::CUDADriver::lin_combo(){
    std::cout<<"\n\nEntering tests for complex number array exptessions."<<std::endl;
     const unsigned int n_rows = 4;
     const unsigned int n_cols = 4;
     const unsigned int n_elements = n_rows * n_cols;

     SciPAL::CudaComplex<Number> alpha(1.1);
     SciPAL::CudaComplex<Number> beta (2.);

     std::vector<SciPAL::CudaComplex<Number> >
                a(n_elements, 1.),
                b(n_elements, std::complex<Number>(2., 2.0));


     for (unsigned int i = 0; i < a.size(); i++ )
         a[i] = SciPAL::CudaComplex<Number>(i+1, (i+1)/2.); //generate some inputs

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vA(n_elements), vB, vC;
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
        std::vector<SciPAL::CudaComplex<Number> > d(n_sin_elements);
       for(uint i = 0; i < d.size(); i++)
           d[i] = std::complex<Number>(i* 2.* M_PI / d.size(), i* 4.* M_PI / d.size()) ;

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vD; vD = d; //(n_sin_elements, 1, d);
      vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
         std::cout << "sin(vD) : " << std::endl;
       vD.print();
       //vC = alpha * sin(vA) + vB;


       std::cout << " ============ linear combination test ======" << std::endl;

       //This does not work : vC = 2.0 * vA; because of mismatching type.
       // Thus we hav to use numbers wrapped in CudaComplexes
       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vE(n_elements);
       vE = beta * vA; ;
       std::cout << "vE = 2.0 * vA" << std::endl;
       vE.print();

       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vF(n_elements);
       vF = vA + vB;
       std::cout << "vC = vA + vB" << std::endl;
       vC.print();


       SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vG(n_elements);
       vG = alpha*vA + beta*vB +alpha*vC + cplxNumber(3.0)*vD + alpha*vE + cplxNumber(4.0)*vF;
       vG.print();



}

void step42::CUDADriver::views(){
std::cout<<"entering tests for Views"<<std::endl;
    // some dummy vectors
    const unsigned int n_rows = 4;
    const unsigned int n_cols = 4;
    const unsigned int n_elements = n_rows * n_cols;

    std::vector<Number>
            a(n_elements, 1.),
            b(n_elements, 2.),
            c(n_rows * n_rows, 1.23);

    for (unsigned int i = 0; i < b.size(); i++ )
        b[i] = i+1;


    SciPAL::Matrix<Number, BW>
            A(n_rows, n_cols, b),
            B(n_cols, n_rows, b),
            C(2, 4, c);

    SciPAL::SubMatrixView<Number, BW> SV_A(A, 0,2,0,4);
    SciPAL::SubMatrixView<Number, BW> SV_B(B, 0,4,0,4);
    SciPAL::SubMatrixView<Number, BW> SV_C(C, 0,2,0,4);

    SV_C = SV_A * SV_B;

    std::cout << "Matrix A "<<std::endl;
    A.print();

    std::cout << "Matrix b "<<std::endl;
    B.print();

    std::cout << "C = A*B "<<std::endl;
    C.print();

}

void step42::CUDADriver::operator_precedence(){
const unsigned int n_rows = 2;
const unsigned int n_cols = 2;
const unsigned int n_elements = n_rows * n_cols;

Number alpha = 1.1;
Number beta = 2.;

std::vector<Number>
           a(n_elements, 1.),
           b(n_elements, 2.);


for (unsigned int i = 0; i < a.size(); i++ )
    a[i] = i+1;

SciPAL::Vector<Number, BW> vA, vB(n_elements), vC, vD(n_elements);
  vA = a;
  // This sets all elements of vB to 2.3, note: vector needs to be initialized.
  vB = SciPAL::Literal<Number>(2.3);
  vC = a;


  std::cout << "vA : " << std::endl;
  vA.print();

  std::cout << "vB : " << std::endl;
  vB.print();

  std::cout << "vC : " << std::endl;
  vC.print();

    std::cout << "vD = (vA && vB ) + vC" << std::endl;
    vD = (vA && vB ) + vC;
    vD.print();

    std::cout << "vA : " << std::endl;
    vA.print();

    std::cout << "vB : " << std::endl;
    vB.print();

    std::cout << "vC : " << std::endl;
    vC.print();
    std::cout << std::endl;

    std::cout << "vD = vC + (vA && vB)" << std::endl;
   vD = vC + vA * vB;
    vD.print();

    SciPAL::Literal<Number> test(0);
    test = (SciPAL::transpose<SciPAL::Vector<Number, BW> >(vA)) * vB;
    std::cout<< test << std::endl;
}
void step42::CUDADriver::stacks_of_LAOs(){

    // some dummy vectors
    const unsigned int n_rows = 4;
    const unsigned int n_cols = 4;
    const unsigned int n_elements = n_rows * n_cols;
//  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
    std::cout <<" Stack of vector test with real numbers\n";
    typedef Number test_nmbr;
    typedef SciPAL::Vector<test_nmbr, BW> test_LAO;
    std::vector<std::vector<test_nmbr>> h_test(3, std::vector<test_nmbr>(n_elements));
    for(auto &i : h_test)
        std::iota(i.begin(), i.end(), 1);

    std::vector<SciPAL::Vector<test_nmbr, BW> > d_test(3, SciPAL::Vector<test_nmbr, BW>(n_elements));

    //element wise copy works
    for(uint ii=0; ii<h_test.size(); ii++)
        d_test[ii] = h_test[ii];

    //or copy whole stacks. note: this converts the std::vector in a SciPAL::Vector
//    d_test = h_test;

    std::cout<<"initialization \n";
    d_test[0].print();
    d_test[2] = 2.0 * d_test[0] + 3.0*d_test[1];
    std::cout<<"result of d_test[2] = 2.0 * d_test[0] + 3.0*d_test[1] \n";
    d_test[2].print();
}
////  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     {
//        std::cout <<" Stack of vector test with cplx numbers\n";
//        typedef cplxNumber test_nmbr;
//        typedef SciPAL::Vector<test_nmbr, BW> test_LAO;
//        std::vector<std::vector<test_nmbr>> h_test(3, std::vector<test_nmbr>(n_elements));
//        for(auto &i : h_test)
//            std::iota(i.begin(), i.end(), 1);

//        SciPAL::Stack<test_LAO> d_test(3, 1);

//        d_test = h_test;

//        std::cout<<"initialization \n";
//        d_test[0].print();
//        d_test[2] = SciPAL::Literal<test_nmbr>(2.0) * d_test[0] + SciPAL::Literal<test_nmbr>(3.0)*d_test[1];
//        std::cout<<"result of d_test[2] = 2.0 * d_test[0] + 3.0*d_test[1] \n";
//        d_test[2].print();
//    }

////  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     {
//        std::cout <<" Stack of mtx test with real numbers\n";
//        typedef Number test_nmbr;
//        typedef SciPAL::Matrix<test_nmbr, BW> test_LAO;
//        typedef SciPAL::Matrix<test_nmbr, blas> test_LAO_h;

//        std::vector<test_nmbr> a(n_elements);

//        std::iota(a.begin(), a.end(), 2);

//        test_LAO_h host_init_mtx(n_rows, n_cols, a);
//        std::cout<<"host init mtx \n";
//        host_init_mtx.print();

//        //this is some how ugly how we have to create the stack of mtx on the host
//        std::vector<test_LAO_h> h_test(3);

//        for(uint ii = 0; ii< h_test.size(); ii++)
//            h_test[ii] = host_init_mtx;

//        SciPAL::Stack<test_LAO> d_test(3,1);
//        d_test = h_test;

//        std::cout<<"host initialization \n";
//        h_test[0].print();

//        std::cout<<"initialization \n";
//        d_test[0].print();
//        d_test[2] = SciPAL::Literal<test_nmbr>(2.0) * d_test[0] + SciPAL::Literal<test_nmbr>(3.0)*d_test[1];
//        std::cout<<"result of d_test[2] = 2.0 * d_test[0] + 3.0*d_test[1] \n";
//        d_test[2].print();
//    }

////  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//         {
//            std::cout <<" Stack of mtx test with cplx numbers\n";
//            typedef cplxNumber test_nmbr;
//            typedef SciPAL::Matrix<test_nmbr, BW> test_LAO;
//            typedef SciPAL::Matrix<test_nmbr, blas> test_LAO_h;

//            std::vector<test_nmbr> a(n_elements);

//            std::iota(a.begin(), a.end(), 2);

//            test_LAO_h host_init_mtx(n_rows, n_cols, a);
//            std::cout<<"host init mtx \n";
//            host_init_mtx.print();

//            //this is some how ugly how we have to create the stack of mtx on the host
//            std::vector<test_LAO_h> h_test(3);

//            for(uint ii = 0; ii< h_test.size(); ii++)
//                h_test[ii] = host_init_mtx;

//            //in this case we CAN NOT write SciPAL::Stack<test_LAO> d_test(3,1);
//            //as for the real case, because internally this maps to a construction
//            //of a dealii::Identity matrix, which is sadly not defined for SciPAL's
//            //complex numbers. The usage of this constructor becomes only a problem
//            //if we try to copy uninitialized Stack elements.
//            SciPAL::Stack<test_LAO> d_test;
//            d_test = h_test;

//            std::cout<<"host initialization \n";
//            h_test[0].print();

//            std::cout<<"initialization \n";
//            d_test[0].print();
//            d_test[2] = SciPAL::Literal<test_nmbr>(2.0) * d_test[0] + SciPAL::Literal<test_nmbr>(3.0)*d_test[1];
//            std::cout<<"result of d_test[2] = 2.0 * d_test[0] + 3.0*d_test[1] \n";
//            d_test[2].print();
//        }
}

#endif //CUDA_DRIVER
