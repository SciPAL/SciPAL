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

Copyright  S. C. Kramer , J. Hagemann  2010 - 2016, T. Köhler 2014 - 2016
*/
namespace SciPAL {
namespace LAOOperations {


void permute_inv_col(const Vector<int,BW>& PVT)
{
    Assert(this->n_cols() == PVT.size(),
           dealii::ExcMessage("Pivot length does not match matrix dimensions."));
    // Swap the cols inversely according to PVT:
    Matrix<T,BW> tmp(this->n_rows(),this->n_cols());
    ColVectorView<T,BW,const Matrix<T,BW>> tmp_col(tmp,0);
    ColVectorView<T,BW,const Matrix<T,BW>> this_col(*this,0);
    for(unsigned int i = 0; i < this->n_cols(); i++)
    {
        int j = PVT(i) - 1;
        this_col.reset(i);
        tmp_col.reset(j);
        // copy col i into col temp j
        tmp_col = this_col;
    }
    *this = tmp;
}


// Permutation of the matrix' rows according to the respective lapack definition in the LU decomposition
// "IPIV is INTEGER array, dimension (min(M,N))
// The pivot indices; for 1 <= i <= min(M,N), row i of the
// matrix was interchanged with row IPIV(i)."

void permute_row(const Vector<int,BW>& PVT)
{
    Assert(this->n_rows() == PVT.size(),
           dealii::ExcMessage("Pivot length does not match matrix dimensions."));
    // Swap the rows according to PVT:
    for(int i = 0; i < this->n_rows(); i++)
    {
        int j = PVT(i) - 1;
        if(i != j)
            BW::swap(this->n_cols(), this->data() + i, this->stride * this->leading_dim, this->data() + j, this->stride * this->leading_dim);
    }
}

void permute_inv_row(const Vector<int,BW>& PVT)
{
    Assert(this->n_rows() == PVT.size(),
           dealii::ExcMessage("Pivot length does not match matrix dimensions."));
    // Swap the rows inversely according to PVT:
    for(int i = (this->n_rows() - 1); i >= 0; i--)
    {
        int j = PVT(i) - 1;
        if(i != j)
            BW::swap(this->n_cols(), this->data() + i, this->stride * this->leading_dim, this->data() + j, this->stride * this->leading_dim);
    }
}

static size_t count_n_rows(const Matrix<T, BW>& src)
{
    return src.n_rows();
}
template<typename... Mtx>
static size_t count_n_rows(const Matrix<T, BW> first, const Mtx&... Matrices)
{
    return first.n_rows() + count_n_rows(Matrices...);
}
static size_t count_n_cols(const Matrix<T, BW>& src)
{
    return src.n_cols();
}
template<typename... Mtx>
static size_t count_n_cols(const Matrix<T, BW> first, const Mtx&... Matrices)
{
    return first.n_cols() + Matrix<T, BW>::count_n_cols(Matrices...);
}

void insert(int current_row, int current_col, const Matrix<T, BW>& src)
{
    //todo assertion for range
    size_t m = src.n_rows();
    size_t n = src.n_cols();
    int next_row = current_row + m;
    int next_col = current_col + n;

    SubMatrixView<T,BW> insertion(*this, current_row, next_row, current_col, next_col);

    insertion = src;
}

template<typename... Mtx>
void insert(int current_row, int current_col,const Matrix<T, BW>& first, const Mtx&... Matrices)
{
    //todo assertion for range
    size_t m = first.n_rows();
    size_t n = first.n_cols();
    int next_row = current_row + m;
    int next_col = current_col + n;

    SubMatrixView<T,BW> insertion(*this, current_row, next_row, current_col, next_col);

    insertion = first;
    this->insert(next_row, next_col, Matrices...);
}

template<typename... Mtx>
void direct_sum(const Mtx&... Matrices)
{

    size_t m = Matrix<T, BW>::count_n_rows(Matrices...);
    size_t n = Matrix<T, BW>::count_n_cols(Matrices...);
    this->reinit(m, n);
    this->insert(0, 0, Matrices...);
}

// Permutation of the matrix' colums according to the respective lapack definition in the QP3 decomposition
// "On exit, if JPVT(J)=K, then the J-th column of A*P was the the K-th column of A."
void permute_col(const Vector<int,BW>& PVT)
{
    Assert(this->n_cols() == PVT.size(),
           dealii::ExcMessage("Pivot length does not match matrix dimensions."));
    // Swap the cols  according to PVT:

    Matrix<T,BW> tmp(this->n_rows(),this->n_cols());
    ColVectorView<T,BW,Matrix<T,BW>> tmp_col(tmp,0);
    ColVectorView<T,BW,Matrix<T,BW>> this_col(this,0);
    for(unsigned int i = 0; i < this->n_cols(); i++)
    {
        int j = PVT(i) - 1;
        this_col.reset(j);
        tmp_col.reset(i);
        // copy col j into col temp i
        tmp_col = this_col;
    }
    *this = tmp;
}

// @sect4{Funktion: get_amax}
//!
//! returns the index of the absolutely largest value
template <typename T, //! numbertype
          typename BW, //! blas type
          typename T_src,
          template <typename, typename, typename> class LAO //! template for result type
          >
unsigned int index_of_largest_element(const LAO<T, BW, T_src>& src) const
{
    return BW::amax(src.size(), src.data(), src.stride);
}

// @sect4{Funktion: get_amax}
//!
//! returns the index of the absolutely largest value
template <typename T, //! numbertype
          typename BW, //! blas type
          template <typename, typename> class LAO //! template for result type
          >
unsigned int index_of_largest_element(const LAO<T, BW>& src) const
{
    return BW::amax(src.size(), src.data(), src.stride);
}

}//LAOOperations

}//SciPAL



//ToDO write Blas cublas wrapper for accum sum?
#pragma omp declare reduction(+ : SciPAL::CudaComplex<float>, SciPAL::CudaComplex<float> : omp_out += omp_in) initializer(omp_priv={0.0, 0.0})

#pragma omp declare reduction(+ : SciPAL::CudaComplex<double>, SciPAL::CudaComplex<double> : omp_out += omp_in) initializer(omp_priv={0.0, 0.0})
