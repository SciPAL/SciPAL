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


#ifndef HADAMARD_H
#define HADAMARD_H

#include <iostream>
#include <vector>
#include <deal.II/lac/full_matrix.h>
#include <lac/FullMatrixAccessor.h>
#include <cmath>

struct MatrixCreator{
    template <typename T>
    static void hadamard(int log2N, dealii::FullMatrix<T>& H);

    template <typename T>
    static void normalized_hadamard(int log2N, dealii::FullMatrix<T>& H);

    template <typename T>
            static void extended_hadamard(size_t n_rows, dealii::FullMatrix<T>& H);

    template <typename T>
            static void extended_normalized_hadamard(size_t n_rows, dealii::FullMatrix<T>& H);

    template <typename T>
    static void vandermonde_row_wise(dealii::FullMatrix<T>& V,const dealii::Vector<T>& x){
        vandermonde(V,x,true);
    }

    template <typename T>
    static void vandermonde_col_wise(dealii::FullMatrix<T>& V,const dealii::Vector<T>& x){
        vandermonde(V,x,false);
    }

    template <typename T>
    static void toepliz(dealii::FullMatrix<T>& M, const dealii::Vector<T>& r);

private:
    template <typename T>
    static void vandermonde(dealii::FullMatrix<T>& V,const dealii::Vector<T>& x, bool rowwise = true);
};

template <typename T>
void MatrixCreator::hadamard(int log2N, dealii::FullMatrix<T>& H){
    int N = std::pow(2.0,log2N);
    H.reinit(N,N);
    H(0,0) = 1;

    for (int n = 1; n < N; n += n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                H(i+n,j)   =  H(i,j);
                H(i,j+n)   =  H(i,j);
                H(i+n,j+n) = -1 * H(i,j);
            }
        }
    }
}

template <typename T>
        void MatrixCreator::extended_normalized_hadamard(size_t n_rows, dealii::FullMatrix<T>& H){
        H.reinit(n_rows,n_rows);
        size_t block_size;
        size_t n_rows_rest=n_rows;
        size_t block_begin=0;
        while(true){
            block_size=std::pow(2.,std::floor(log2(n_rows_rest)));
            n_rows_rest-=block_size;

            dealii::FullMatrix<T> H_tmp(block_size);
            normalized_hadamard(log2(block_size),H_tmp);
            H.fill(H_tmp,block_begin,block_begin,0,0);
            block_begin+=block_size;

            if(n_rows_rest==0)
                break;
        }
    }

template <typename T>
        void MatrixCreator::extended_hadamard(size_t n_rows, dealii::FullMatrix<T>& H){
        H.reinit(n_rows,n_rows);
        size_t block_size;
        size_t n_rows_rest=n_rows;
        size_t block_begin=0;
        while(true){
            block_size=std::pow(2.,std::floor(log2(n_rows_rest)));
            n_rows_rest-=block_size;

            dealii::FullMatrix<T> H_tmp(block_size);
            hadamard(log2(block_size),H_tmp);
            H.fill(H_tmp,block_begin,block_begin,0,0);
            block_begin+=block_size;

            if(n_rows_rest==0)
                break;
        }

    }

template <typename T>
void MatrixCreator::normalized_hadamard(int log2N, dealii::FullMatrix<T>& H)
{
    MatrixCreator::hadamard(log2N, H);

    H /= std::sqrt(H.n_rows() );
}



template <typename T>
void MatrixCreator::vandermonde(dealii::FullMatrix<T>& V,const dealii::Vector<T>& x, bool rowwise){
    int N = x.size();
    V.reinit(N,N);
    if(rowwise == false){
        for(int i = 0; i<N ; i++){
            for(int j = 0; j<N; j++){
                V(i,j) = pow(x(i),j);
            }
        }
    }
    if(rowwise == true){
        for(int i = 0; i<N ; i++){
            for(int j = 0; j<N; j++){
                V(j,i) = pow(x(i),j);
            }
        }
    }
}

template <typename T>
void MatrixCreator::toepliz(dealii::FullMatrix<T>& M, const dealii::Vector<T>& r){
    if (r.size()%2 == 0){
        std::cout<<"Toepliz: Nur ungerade Anzahl an Elemente sind Erlaubt"<<std::endl;
    }else{
        int N = (r.size()/2)+1;
        M.reinit(N,N);
        for(int i = 0; i<N; i++){
            for (int j = 0; j<N; j++){
                M(i,j) = r(((N-1)-j)+i);
            }
        }
    }


}

#endif //! HADAMARD_H
