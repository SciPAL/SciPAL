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


#ifndef PreconditionerAINV_HH
#define PreconditionerAINV_HH

#include <lac/PreconditionerAINV.h>
#include <cuda_kernel_wrapper_step-8.cu.h>
#include <lac/SparseMatrixAccessor.hh>

//Constructors and initializers...

template<typename Number>
PreconditionerAINV<Number>::PreconditionerAINV(){

}

template<typename Number>
void PreconditionerAINV<Number>::initialize (const MATRIX &_A)//!,
               //!  typename BaseClass::AdditionalData parameters = typename BaseClass::AdditionalData())
{
    SparseMatrixAccessor<Number> accessable_A(_A);

    //! TODO: Transpose sparsity pattern of M_r / M_l



    std::vector<Number> vectorOfZeros(_A.n_nonzero_elements());
    std::vector<Number> dVector(_A.m());

    for (int i=0; i<_A.m(); i++)
    {
        dVector[i]=0;
    }

    for (int i=0; i<_A.n_nonzero_elements(); i++)
    {
        vectorOfZeros[i]=0;
    }


    this->A.reinit(accessable_A.get_values(),_A.get_sparsity_pattern());

    CUDASparseMatrix AT(this->A);
    AT.transpose();

    this->M_l.reinit(&vectorOfZeros[0],AT.get_sparsity_pattern());
    this->M_r.reinit(&vectorOfZeros[0],AT.get_sparsity_pattern());

    CUDASparseMatrix Z(&vectorOfZeros[0],AT.get_sparsity_pattern());

    dealii::SparsityPattern DPattern(_A.m(),_A.n(),1,true);

    DPattern.compress();
    CUDASparseMatrix D(&dVector[0], DPattern);

    //!Calculating M_r

    this->A.AINVPrecond(Z,D);

    D.diag_invert();
    Z.Mmult(this->M_r,D);

    AT.AINVPrecond(this->M_l,D);
    M_r.transpose();

}


template<typename Number>
void PreconditionerAINV<Number>::initialize (const SparseMatrix &_SA)//!,
               //!  typename BaseClass::AdditionalData parameters = typename BaseClass::AdditionalData())
{

    MATRIX _A;
#warning eltdim sollte hier noch gesetzt werden
    _A.reinit(_SA.get_sparsity_pattern(),1);  //eltdim = 1 atm; has to be changed in future versions
    _A.copy_from(_SA);


    SparseMatrixAccessor<Number> accessable_A(_A);

    //! TODO: Transpose sparsity pattern of M_r / M_l


    //dummy-vector
    std::vector<Number> vectorOfZeros(_A.n_nonzero_elements());
    std::vector<Number> dVector(_A.m());

    //fill dummy vectors with zeros
    for (int i=0; i<_A.m(); i++)
    {
        dVector[i]=0;
    }

    for (int i=0; i<_A.n_nonzero_elements(); i++)
    {
        vectorOfZeros[i]=0;
    }


    this->A.reinit(accessable_A.get_values(),_A.get_sparsity_pattern());

    CUDASparseMatrix AT(this->A);
    AT.transpose();

    this->M_l.reinit(&vectorOfZeros[0],AT.get_sparsity_pattern());
    this->M_r.reinit(&vectorOfZeros[0],AT.get_sparsity_pattern());

    CUDASparseMatrix Z(&vectorOfZeros[0],AT.get_sparsity_pattern());

    dealii::SparsityPattern DPattern(_A.m(),_A.n(),1,true);

    DPattern.compress();
    CUDASparseMatrix D(&dVector[0], DPattern);

    //!Calculating M_r

    this->A.AINVPrecond(Z,D);

    D.diag_invert();
    Z.Mmult(this->M_r,D);

    AT.AINVPrecond(this->M_l,D);
    M_r.transpose();

}


//----------------------------------------------------------


template<typename Number>
void PreconditionerAINV<Number>::AINVPrecond( CUDASparseMatrix &B, CUDASparseMatrix &d) {

    step8::SpMVKernels<Number> spmv;
    const int BLOCK_SIZE=256;

    spmv.AINV(this->_nrows,BLOCK_SIZE, d._nonZeroValues,
              this->_nonZeroValues, this->_columnIndices, this->_rowStart,
              B._nonZeroValues,B._columnIndices,B._rowStart);


printf("NNZ B=%i\n",B._nnz);

}


// @sect4{Funktion: vmult}
//!
//! Diese Funktion steuert den Aufruf des Kernels,
//! der die Matrix-Vektor Multiplikation ausf&uuml;hrt.
template<typename Number>
void PreconditionerAINV<Number>::vmult(CUDAVector &dst, const CUDAVector &src) const
{

    step8::SpMVKernels<Number> SpMV;
        //! Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
        //! Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
        CUDAVector y = dst; //!(&dst(0), dst.size());
        CUDAVector x = src; //!(&src(0), src.size());


        const unsigned int NUM_BLOCKS = this->_nrows;
        const unsigned int BLOCK_SIZE = 32;
        //!const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

        //! TO DO: Werden genuegend Bloecke gestartet?

        //! Berechne das Ergebnis der Matrix-Vektor Multiplikation
        //! mit Hilfe der Kernel-Funktion auf der GPU.

        //! prec(A)*p = p_
        SpMV.csr(
                this->M_l.getNrows(),           //M_l oder M_r zuerst?
                this->M_l.getRowStart(),
                this->M_l.getColumnIndices(),
                this->M_l.getNonZeroValues(),
                x.getValues_d(),
                y.getValues_d());

        //! p__ = A*p_

        A.vmult(x,y);

        //! y = prec(A)*p__
        SpMV.csr(
                this->M_r.getNrows(),
                this->M_r.getRowStart(),
                this->M_r.getColumnIndices(),
                this->M_r.getNonZeroValues(),
                x.getValues_d(),
                y.getValues_d());

        //! Speichere Ergebnis auf dem Speicher der CPU.
        y.copyToHost();
}




template<typename Number>
void PreconditionerAINV<Number>::postprocessor (CUDAVector &dst, const CUDAVector &src) const
{

    step8::SpMVKernels<Number> SpMV;
    //! Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
    //! Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
    CUDAVector y=dst;//!(&dst(0), dst.size());
    CUDAVector x=src;//!(&src(0), src.size());


    const unsigned int NUM_BLOCKS = this->_nrows;
    const unsigned int BLOCK_SIZE = 32;
    //!const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

    //! TO DO: Werden genuegend Bloecke gestartet?

    //! Berechne das Ergebnis der Matrix-Vektor Multiplikation
    //! mit Hilfe der Kernel-Funktion auf der GPU.


    SpMV.csr(
            this->M_r.getNrows(),
            this->M_r.getRowStart(),
            this->M_r.getColumnIndices(),
            this->M_r.getNonZeroValues(),
            x.getValues_d(),
            y.getValues_d());

    //! Speichere Ergebnis auf dem Speicher der CPU.
    y.copyToHost();
}


// @sect4{Funktion: mmult}
//!
//! Diese Funktion steuert den Aufruf des Kernels,
//! der die Matrix-Matrix Multiplikation ausf&uuml;hrt
//! The source Matrix has to be transposed!!!
template<typename Number>
void PreconditionerAINV<Number>::Mmult(CUDASparseMatrix &dst, CUDASparseMatrix &src)
{
    step8::SpMVKernels<Number> SpMV;

    //!Values is going to be divided into Blocks each of which has a block size (which equals the threads per block) defined in BLOCK_SIZE
    const unsigned int BLOCK_SIZE = 2;
    //!then the number of blocks can be calculated by dividing the number of nonzero elements by the block size
    unsigned int NUM_BLOCKS = this->_nnz / BLOCK_SIZE;
    if (!NUM_BLOCKS) NUM_BLOCKS = 1;
    //!These definitions above have no exception catching procedures, this might be a TODO

    //!Call of the corresponding Matrix-Multiplication Kernel wrapper in spmv
    SpMV.CSRApproxMMult(NUM_BLOCKS,BLOCK_SIZE, this->_nonZeroValues, this->_columnIndices, this->_rowStart,
                        src._nonZeroValues, src._columnIndices, src._rowStart,
                        dst._nonZeroValues, dst._columnIndices, dst._rowStart);

}

template<typename Number>
void PreconditionerAINV<Number>::transpose()
{
    unsigned int n_rows = this->_nrows;

    //! Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Nicht-null-Elemente ...
    size_t sizeValVector = sizeof(Number) * this->_nnz;
    size_t sizeRowStart = sizeof(int) * (n_rows+1);
    size_t sizeColInd = sizeof(int) * this->_nnz;
    Number v_h[this->_nnz];
    int rowStart[n_rows+1];
    int columnIndices[this->_nnz];

    //! ... und weise den Speicherbereich zu.
   //! cudaMalloc((void**)&this->_nonZeroValues, sizeValVector);
    cudaMemcpy(v_h, this->_nonZeroValues , sizeValVector, cudaMemcpyDeviceToHost);
    cudaMemcpy(rowStart, this->_rowStart, sizeRowStart, cudaMemcpyDeviceToHost);
    cudaMemcpy(columnIndices, this->_columnIndices, sizeColInd, cudaMemcpyDeviceToHost);





    //! !!!ONLY DEFINED FOR N x N MARICES!!!
    dealii::SparsityPattern newPat(n_rows,n_rows,n_rows,false);
    Number tmp_values[this->_nnz];
    std::vector<int> counter(n_rows);


    //!init
    for (std::vector<int>::iterator k=counter.begin(); k<=counter.end(); k++)
        *k = 0;

    //!building new pattern
    int row=0;
    int col;

    for (int i=0; i<this->_nnz;i++)
    {

        while (rowStart[row]<=i) row++; row--;
        col = columnIndices[i];
        newPat.add(col,row);
    }
    newPat.compress();
    const size_t *startrow =newPat.get_rowstart_indices();
    int posInVal;

    for (int i=0; i<this->_nnz;i++)
    {
        col = columnIndices[i];
        posInVal = startrow[col];
        posInVal += counter[col];
        counter[col]++;
        tmp_values[posInVal] = v_h[i];
    }
    //!for (int i=0; i<_nnz; i++) printf(" %.3f \n",tmp_values[i]);
    //!Rekonstrukt the Matrix
    //! Speichere die Anzahl der Zeilen
    //! und der Nicht-null-Elemente
    this->_nrows = newPat.n_rows();

    //! ... und weise den Speicherbereich zu.
    cudaMalloc((void**)&this->_nonZeroValues, sizeValVector);
    cudaMemcpy(this->_nonZeroValues, tmp_values, sizeValVector, cudaMemcpyHostToDevice);

    //! Allociere den ben&ouml;tigten Speicher auf der CPU Seite ...
    std::vector<int> rowSt(this->_nrows+1);

    //! Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Zeilenanf&auml;nge.
    size_t sizeRowVector = sizeof(int) * rowSt.size();

    //! ... und weise den Speicherbereich auf der GPU Seite zu.
    cudaMalloc((void**)&this->_rowStart, sizeRowVector);

    //! Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Spaltenindizes.
    size_t sizeColVector = sizeof(int) * this->_nnz;

    //! Allociere den ben&ouml;tigten Speicher auf der CPU Seite ...
    std::vector<int> colIn(this->_nnz);

    //! ... und teile den Speicherbereich auf der GPU zu.
    cudaMalloc((void**)&this->_columnIndices, sizeColVector);

    //! Initialisiere tempor&auml;re Variablen.
    int rowLength = 0;
    int rowOffset = 0;

    //! Iteriere &uuml;ber alle Zeilen.
    for(int i = 0; i < this->_nrows; i++)
    {
            //! Speichere die L&auml;nge der aktuellen Zeile.
            rowLength = newPat.row_length(i);

            //! Speichere den Zeiger auf das erste Element.
            rowSt[i] = rowOffset;

            //! Iteriere &uuml;ber die Anzahl der Zeileneintr&auml;ge.
            for(int j = 0; j < rowLength; j++)
            {
                    //! Speichere den Spaltenindex des aktuellen Elements.
                    colIn[rowOffset+j] = newPat.column_number(i, j);
            }

            //! Aktualisiere die momentane Position.
            rowOffset += rowLength;
    }

    //! F&uuml:llen des letzten Elements mit der Anzahl der Nicht-null-Elemente.
    rowSt[this->_nrows] = this->_nnz;

    //! Gebe die berechneten Zeilenanf&auml;nge und die Spaltenindizes an die GPU weiter.
    cudaMemcpy(this->_rowStart, &rowSt[0], sizeRowVector, cudaMemcpyHostToDevice);
    cudaMemcpy(this->_columnIndices, &colIn[0], sizeColVector, cudaMemcpyHostToDevice);

}




// @sect4{Transponierte Matrix-Vektor Multiplikation: CudaCSRMatrixView}
//!
//! Diese Funktion steuert den Aufruf des Kernels,
//! der die transponierte Matrix-Vektor Multiplikation ausf&uuml;hrt.

template<typename Number>
void PreconditionerAINV<Number>::Tvmult (CUDAVector &dst, const CUDAVector &src) const
{
        //! do nothing
}



#endif //! PreconditionerAINV_HH
