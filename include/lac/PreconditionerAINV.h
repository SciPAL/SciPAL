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


#ifndef PreconditionerAINV_H
#define PreconditionerAINV_H

#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include "../../testsite/step-8/cuda_driver_step-8.h"
#include <base/subscriptor.h>
#include "../../testsite/step-8/csrmatrixview.h"
#include <cutil.h>
#include <lac/cublas_BlancMatrix.h>



template<typename Number>
class PreconditionerAINV : public virtual dealii::Subscriptor
{
        public:

                typedef dealii::Vector<Number> VECTOR;
                typedef typename step8::CUDAVectorView<Number>::Type CUDAVector;
                typedef typename step8::CudaCSRMatrixView<Number>::Type CUDASparseMatrix;
                typedef BlancMatrix<Number> MATRIX;
                typedef dealii::SparseMatrix<Number> SparseMatrix;


                PreconditionerAINV();


                void initialize (const MATRIX &A);
                                 //!typename BaseClass::AdditionalData parameters = typename BaseClass::AdditionalData());
                void initialize (const SparseMatrix &A);

                void Tvmult (CUDAVector &dst, const CUDAVector &src) const;

                void vmult (CUDAVector &dst, const CUDAVector &src) const;

                void postprocessor (CUDAVector &dst, const CUDAVector &src) const;

                void Mmult(CUDASparseMatrix &dst, CUDASparseMatrix &src);

                void transpose();

                void AINVPrecond( CUDASparseMatrix &B, CUDASparseMatrix &d);



        protected:
                Number * _nonZeroValues;
                int * _columnIndices;
                int * _rowStart;

                int _nnz;
                int _nrows;

                MATRIX *inputMatrix;
                CUDASparseMatrix M_r;
                CUDASparseMatrix M_l;
                CUDASparseMatrix A;

};



#endif //! PreconditionerAINV_H
