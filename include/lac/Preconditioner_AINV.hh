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


#ifndef PRECONDITIONER_AINV_HH
#define PRECONDITIONER_AINV_HH

#include "csrmatrixview.h"

// Einbinden der ben&ouml<tigten Kernel.

#include "cuda_kernel_wrapper_step-8.cu.h"
#include <lac/Preconditioner_AINV.h>
#include "lac/SparseMatrixAccessor.h"

//const int MAX_THREADS = (30 * 1024);
const int WARP_SIZE = 32;

// @sect4{Konstruktor: CudaCSRMatrixView}
//
// Der Konstruktor der CudaCSRMatrixView Klasse allokiert den Speicher
// auf der GPU und kopiert die Daten vom host auf die GPU.
// @param v_h : Werte-Vektor
// @param sp : SparsityPattern Objekt


template<typename Number>
Preconditioner_AINV<Number>::Preconditioner_AINV(){

}


template<typename Number>
void Preconditioner_AINV<Number>::AINVInvert(CudaCSRMatrixView<Number> &Z, CudaCSRMatrixView &Wtransposed, CudaCSRMatrixView<Number> &D)
{

    CudaCSRMatrixView<Number> AT(*this);
    AT.transpose();


    AT.AINVPrecond(Wtransposed,D);
    this->AINVPrecond(Z,D);


    Wtransposed.transpose();
    Z.transpose();

}

template<typename Number>
void Preconditioner_AINV<Number>::AINVPrecond( CudaCSRMatrixView<Number> &B, CudaCSRMatrixView<Number> &d) {

    SpMVKernels<Number> spmv;
    const int BLOCK_SIZE=256;

    spmv.AINV(this->_nrows,BLOCK_SIZE, d._nonZeroValues,
              this->_nonZeroValues, this->_columnIndices, this->_rowStart,
              B._nonZeroValues,B._columnIndices,B._rowStart);


printf("NNZ B=%i\n",B._nnz);

}



// @sect4{Funktion: vmult}
//
// Diese Funktion steuert den Aufruf des Kernels,
// der die Matrix-Vektor Multiplikation ausf&uuml;hrt.
template<typename Number>
void Preconditioner_AINV<Number>::vmult_right(CUDAVector &dst, const CUDAVector &src) const
{

    SpMVKernels<Number> SpMV;
        // Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
        // Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
        CUDAVectorView<double> y(&dst(0), dst.size());
        CUDAVectorView<double> x(&src(0), src.size());


        const unsigned int NUM_BLOCKS = this->_nrows;
        const unsigned int BLOCK_SIZE = 32;
        //const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

        // TO DO: Werden genuegend Bloecke gestartet?

        // Berechne das Ergebnis der Matrix-Vektor Multiplikation
        // mit Hilfe der Kernel-Funktion auf der GPU.


        SpMV.csr(
                this->M_r->_nrows,
                this->M_r->_rowStart,
                this->M_r->_columnIndices,
                this->M_r->_nonZeroValues,
                x.__values_d,
                y.__values_d);

        // Speichere Ergebnis auf dem Speicher der CPU.
        y.copyToHost();
}

// @sect4{Funktion: vmult}
//
// Diese Funktion steuert den Aufruf des Kernels,
// der die Matrix-Vektor Multiplikation ausf&uuml;hrt.
template<typename Number>
void Preconditioner_AINV<Number>::vmult_left(CUDAVector &dst, const CUDAVector &src) const
{

    SpMVKernels<Number> SpMV;
        // Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
        // Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
        CUDAVectorView<double> y(&dst(0), dst.size());
        CUDAVectorView<double> x(&src(0), src.size());


        const unsigned int NUM_BLOCKS = this->_nrows;
        const unsigned int BLOCK_SIZE = 32;
        //const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

        // TO DO: Werden genuegend Bloecke gestartet?

        // Berechne das Ergebnis der Matrix-Vektor Multiplikation
        // mit Hilfe der Kernel-Funktion auf der GPU.
        SpMV.csr(
                this->M_l->_nrows,
                this->M_l->_rowStart,
                this->M_l->_columnIndices,
                this->M_l->_nonZeroValues,
                x.__values_d,
                y.__values_d);

        // Speichere Ergebnis auf dem Speicher der CPU.
        y.copyToHost();
}

// @sect4{Funktion: vmult}
//
// Diese Funktion steuert den Aufruf des Kernels,
// der die Matrix-Vektor Multiplikation ausf&uuml;hrt.
template<typename Number>
void Preconditioner_AINV<Number>::vmult(CUDAVector &dst, const CUDAVector &src, const struct additional_data) const
{

    SpMVKernels<Number> SpMV;
        // Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
        // Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
        CUDAVectorView<double> y(&dst(0), dst.size());
        CUDAVectorView<double> x(&src(0), src.size());


        const unsigned int NUM_BLOCKS = this->_nrows;
        const unsigned int BLOCK_SIZE = 32;
        //const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

        // TO DO: Werden genuegend Bloecke gestartet?

        // Berechne das Ergebnis der Matrix-Vektor Multiplikation
        // mit Hilfe der Kernel-Funktion auf der GPU.

        // prec(A)*p = p_
        SpMV.csr(
                this->M_l->_nrows,
                this->M_l->_rowStart,
                this->M_l->_columnIndices,
                this->M_l->_nonZeroValues,
                x.__values_d,
                y.__values_d);

        // p__ = A*p_
        A->vmult(x,y);


        // y = prec(A)*p__
        SpMV.csr(
                this->M_r->_nrows,
                this->M_r->_rowStart,
                this->M_r->_columnIndices,
                this->M_r->_nonZeroValues,
                x.__values_d,
                y.__values_d);

        // Speichere Ergebnis auf dem Speicher der CPU.
        y.copyToHost();
}


template<typename Number>
void Preconditioner_AINV<Number>::postprocessor() (CUDAVector &dst, const CUDAVector &src, const struct additional_data) const
{

    SpMVKernels<Number> SpMV;
    // Die &uuml;bergebenen Vektoren sollen deal.II Vektoren sein.
    // Wir m&uuml;ssen also zun&auml;chst CUDAVector daraus machen.
    CUDAVectorView<double> y(&dst(0), dst.size());
    CUDAVectorView<double> x(&src(0), src.size());


    const unsigned int NUM_BLOCKS = this->_nrows;
    const unsigned int BLOCK_SIZE = 32;
    //const unsigned int NUM_BLOCKS = std::min(MAX_BLOCKS, (this->_nrows + BLOCK_SIZE - 1)/BLOCK_SIZE);

    // TO DO: Werden genuegend Bloecke gestartet?

    // Berechne das Ergebnis der Matrix-Vektor Multiplikation
    // mit Hilfe der Kernel-Funktion auf der GPU.


    SpMV.csr(
            this->M_r->_nrows,
            this->M_r->_rowStart,
            this->M_r->_columnIndices,
            this->M_r->_nonZeroValues,
            x.__values_d,
            y.__values_d);

    // Speichere Ergebnis auf dem Speicher der CPU.
    y.copyToHost();
}


// @sect4{Funktion: mmult}
//
// Diese Funktion steuert den Aufruf des Kernels,
// der die Matrix-Matrix Multiplikation ausf&uuml;hrt
// The source Matrix has to be transposed!!!
template<typename Number>
void Preconditioner_AINV<Number>::Mmult(CudaCSRMatrixView<Number> &dst, CudaCSRMatrixView<Number> &src)
{
    SpMVKernels<Number> SpMV;

    //Values is going to be divided into Blocks each of which has a block size (which equals the threads per block) defined in BLOCK_SIZE
    const unsigned int BLOCK_SIZE = 2;
    //then the number of blocks can be calculated by dividing the number of nonzero elements by the block size
    unsigned int NUM_BLOCKS = this->_nnz / BLOCK_SIZE;
    if (!NUM_BLOCKS) NUM_BLOCKS = 1;
    //These definitions above have no exception catching procedures, this might be a TODO

    //Call of the corresponding Matrix-Multiplication Kernel wrapper in spmv
    SpMV.CSRApproxMMult(NUM_BLOCKS,BLOCK_SIZE, this->_nonZeroValues, this->_columnIndices, this->_rowStart,
                        src._nonZeroValues, src._columnIndices, src._rowStart,
                        dst._nonZeroValues, dst._columnIndices, dst._rowStart);

}

template<typename Number>
void Preconditioner_AINV<Number>::transpose()
{
    unsigned int n_rows = this->_nrows;

    // Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Nicht-null-Elemente ...
    size_t sizeValVector = sizeof(Number) * this->_nnz;
    size_t sizeRowStart = sizeof(int) * (n_rows+1);
    size_t sizeColInd = sizeof(int) * this->_nnz;
    Number v_h[this->_nnz];
    int rowStart[n_rows+1];
    int columnIndices[this->_nnz];

    // ... und weise den Speicherbereich zu.
   // cudaMalloc((void**)&this->_nonZeroValues, sizeValVector);
    cudaMemcpy(v_h, this->_nonZeroValues , sizeValVector, cudaMemcpyDeviceToHost);
    cudaMemcpy(rowStart, this->_rowStart, sizeRowStart, cudaMemcpyDeviceToHost);
    cudaMemcpy(columnIndices, this->_columnIndices, sizeColInd, cudaMemcpyDeviceToHost);





    // !!!ONLY DEFINED FOR N x N MARICES!!!
    dealii::SparsityPattern newPat(n_rows,n_rows,n_rows,false);
    Number tmp_values[this->_nnz];
    std::vector<int> counter(n_rows);


    //init
    for (std::vector<int>::iterator k=counter.begin(); k<=counter.end(); k++)
        *k = 0;

    //building new pattern
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
    //for (int i=0; i<_nnz; i++) printf(" %.3f \n",tmp_values[i]);
    //Rekonstrukt the Matrix
    // Speichere die Anzahl der Zeilen
    // und der Nicht-null-Elemente
    this->_nrows = newPat.n_rows();

    // ... und weise den Speicherbereich zu.
    cudaMalloc((void**)&this->_nonZeroValues, sizeValVector);
    cudaMemcpy(this->_nonZeroValues, tmp_values, sizeValVector, cudaMemcpyHostToDevice);

    // Allociere den ben&ouml;tigten Speicher auf der CPU Seite ...
    std::vector<int> rowSt(this->_nrows+1);

    // Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Zeilenanf&auml;nge.
    size_t sizeRowVector = sizeof(int) * rowSt.size();

    // ... und weise den Speicherbereich auf der GPU Seite zu.
    cudaMalloc((void**)&this->_rowStart, sizeRowVector);

    // Berechne die Gr&ouml;&szlig;e des ben&ouml;tigten Speichers f&uuml;r die Spaltenindizes.
    size_t sizeColVector = sizeof(int) * this->_nnz;

    // Allociere den ben&ouml;tigten Speicher auf der CPU Seite ...
    std::vector<int> colIn(this->_nnz);

    // ... und teile den Speicherbereich auf der GPU zu.
    cudaMalloc((void**)&this->_columnIndices, sizeColVector);

    // Initialisiere tempor&auml;re Variablen.
    int rowLength = 0;
    int rowOffset = 0;

    // Iteriere &uuml;ber alle Zeilen.
    for(int i = 0; i < this->_nrows; i++)
    {
            // Speichere die L&auml;nge der aktuellen Zeile.
            rowLength = newPat.row_length(i);

            // Speichere den Zeiger auf das erste Element.
            rowSt[i] = rowOffset;

            // Iteriere &uuml;ber die Anzahl der Zeileneintr&auml;ge.
            for(int j = 0; j < rowLength; j++)
            {
                    // Speichere den Spaltenindex des aktuellen Elements.
                    colIn[rowOffset+j] = newPat.column_number(i, j);
            }

            // Aktualisiere die momentane Position.
            rowOffset += rowLength;
    }

    // F&uuml:llen des letzten Elements mit der Anzahl der Nicht-null-Elemente.
    rowSt[this->_nrows] = this->_nnz;

    // Gebe die berechneten Zeilenanf&auml;nge und die Spaltenindizes an die GPU weiter.
    cudaMemcpy(this->_rowStart, &rowSt[0], sizeRowVector, cudaMemcpyHostToDevice);
    cudaMemcpy(this->_columnIndices, &colIn[0], sizeColVector, cudaMemcpyHostToDevice);

}


template<typename Number>
void Preconditioner_AINV<Number>::initialize (const MATRIX &A)//,
               //  typename BaseClass::AdditionalData parameters = typename BaseClass::AdditionalData())
{
    SparseMatrixAccessor<Number> accessable_A(A);
    this->A->reinit(accessable_A.get_values(),accessable_A.get_sparsity_pattern());
    // TODO: Transpose sparsity pattern of M_r / M_l

    std::vector<Number> vectorOfZeros(A.n_nonzero_elements());
    std::vector<Number> dVector(accessable_A.get_cols());

    for (int i=0; i<accessable_A.get_cols(); i++)
    {
        dVector[i]=0;
    }

    for (int i=0; i<A.n_nonzero_elements(); i++)
    {
        vectorOfZeros[i]=0;
    }

    this->M_l->reinit(vectorOfZeros,A.get_sparsity_pattern());
    this->M_r->reinit(vectorOfZeros,A.get_sparsity_pattern());

    CudaCSRMatrixView<Number> AT(this->A);
    AT.transpose();

    CudaCSRMatrixView<Number> Z(vectorOfZeros,A.get_sparsity_pattern());

    dealii::SparsityPattern DPattern(accessable_A.get_cols(),accessable_A.get_cols(),1,true);

    DPattern.compress();
    CudaCSRMatrixView<Number> D(dVector, DPattern);

    //Calculating M_r

    this->A->AINVPrecond(Z,D);

    D.diagInvert;
    Z.Mmult(M_r,D);

    AT.AINVPrecond(M_l,D);


}
// @sect4{Transponierte Matrix-Vektor Multiplikation: CudaCSRMatrixView}
//
// Diese Funktion steuert den Aufruf des Kernels,
// der die transponierte Matrix-Vektor Multiplikation ausf&uuml;hrt.
template<typename Number>
void Preconditioner_AINV<Number>::Tvmult (CUDAVector &dst, const CUDAVector &src) const
{
        // do nothing
}



#endif // PRECONDITIONER_AINV_HH
