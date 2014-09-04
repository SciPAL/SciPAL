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


#ifndef cublas_BlancMatrix_H
#define cublas_BlancMatrix_H


#include <iomanip>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>

#include <QtGlobal>
#include <QFile>
#include <QTextStream>

#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>



using namespace dealii;

template <typename number>
        class BlancMatrix  : public virtual Subscriptor
{
public:

    typedef number value_type;
    typedef typename numbers::NumberTraits<number>::real_type real_type;

    /**
     * @name Constructors and initalization.
     */


    BlancMatrix();
    BlancMatrix( const BlancMatrix& );
    explicit BlancMatrix( const SparsityPattern &sparsity, const int eltdim );

    virtual ~BlancMatrix();

    template<typename somenumber> BlancMatrix<number> & copy_from (const BlancMatrix<somenumber> &matrix);
    template<typename somenumber> BlancMatrix<number> & copy_from (const SparseMatrix<somenumber> &matrix);


    BlancMatrix<number>& operator = (const BlancMatrix<number> &);
    BlancMatrix<number> & operator= (const IdentityMatrix  &id);
    BlancMatrix & operator = (const double d);

    virtual void reinit (const SparsityPattern &sparsity, const int eltdim);
    virtual void clear ();


    /**
     * @name Information on the matrix
     */

    bool empty () const;
    unsigned int m () const;
    unsigned int n () const;
    unsigned int n_nonzero_elements () const;

    unsigned int dimension() const;

    unsigned int n_actually_nonzero_elements (const double threshold = 0.) const;
    const SparsityPattern & get_sparsity_pattern () const;
    unsigned int memory_consumption () const;

    /**
     * @name Modifying entries
     */

    void set (const unsigned int i,
              const unsigned int j,
              const unsigned int k,
              const unsigned int l,
              const number value);

    void add (const unsigned int i,
              const unsigned int j,
              const unsigned int k,
              const unsigned int l,
              const number value);

    BlancMatrix & operator *= (const number factor);
    BlancMatrix & operator /= (const number factor);


    template <typename somenumber>
            void add (const number factor,
                      const BlancMatrix<somenumber> &matrix);


    /**
     * @name Entry Access
     */

    number global_entry (const unsigned int j) const;

    number operator () (const unsigned int i,
                        const unsigned int j,
                        const unsigned int k,
                        const unsigned int l ) const;

    number el (const unsigned int i,
               const unsigned int j,
               const unsigned int k,
               const unsigned int l
               ) const;

    number diag_element (const unsigned int i) const;
    number & diag_element (const unsigned int i);

    /**
     * @name Matrix vector multiplications
     */

    template <class OutVector, class InVector>
            void vmult (OutVector& dst,
                        const InVector& src) const;



    template <typename somenumber>
            void residual (Vector<somenumber>       &dst,
                           const Vector<somenumber> &x,
                           const Vector<somenumber> &b) const;


    /**
     * @name Matrix norms
     */

    real_type l1_norm () const;
    real_type linfty_norm () const;
    real_type frobenius_norm () const;

    /**
     * @name Preconditioning methods
     */



    template <typename somenumber>
            void precondition_Jacobi (Vector<somenumber>       &dst,
                                      const Vector<somenumber> &src,
                                      const number              omega = 1.) const;






    /**
     * @name Input/Output
     */

    void print (std::ostream &out) const;
    void print_formatted (std::ostream       &out,
                          const unsigned int  precision   = 3,
                          const bool          scientific  = true,
                          const unsigned int  width       = 0,
                          const char         *zero_string = " ",
                          const double        denominator = 1.) const;

    void print_pattern(std::ostream& out,
                       const double threshold = 0.) const;


    /** @addtogroup Exceptions
     * @{ */

    /**
     * Exception
     */


    DeclException2 (ExcInvalidIndex,
                    int, int,
                    << "The entry with index <" << arg1 << ',' << arg2
                    << "> does not exist.");
    /**
                                        * Exception
                                        */
    DeclException1 (ExcInvalidIndex1,
                    int,
                    << "The index " << arg1 << " is not in the allowed range.");
    /**
                                        * Exception
                                        */
    DeclException0 (ExcDifferentSparsityPatterns);
    /**
                                        * Exception
                                        */
    DeclException2 (ExcIteratorRange,
                    int, int,
                    << "The iterators denote a range of " << arg1
                    << " elements, but the given number of rows was " << arg2);
    /**
                                        * Exception
                                        */
    DeclException0 (ExcSourceEqualsDestination);


    void calc_inverse() const;
    void clear_inverse();

protected:
    number *val;

private:

    SmartPointer<const SparsityPattern > cols;

    unsigned int max_len;
    int eltdim;

    number* inv_val;

    mutable bool calced_inverse;

    template <class OutVector, class InVector>
            void threaded_vmult (OutVector       &dst,
                                 const InVector &src,
                                 const int        range_begin,
                                 const int        range_end
                                 ) const;

    template <typename somenumber>
            void threaded_residual (Vector<somenumber>       &dst,
                                    const Vector<somenumber> &u,
                                    const Vector<somenumber> &b,
                                    const int        range_begin,
                                    const int        range_end
                                    ) const ;


};



/*---------------------- Inline functions -----------------------------------*/
template <typename number>
        inline
        unsigned int BlancMatrix<number>::m () const
{
    Assert (cols != 0, ExcNotInitialized());
    return cols->n_rows();
}


template <typename number>
        inline
        unsigned int BlancMatrix<number>::n () const
{
    Assert (cols != 0, ExcNotInitialized());
    return cols->n_cols();
}


template <typename number>
        inline
        unsigned int BlancMatrix<number>::dimension () const
{

    return this->eltdim;
}

//! Inline the set() and add()
//! functions, since they will be
//! called frequently.

template <typename number>
        inline
        void
        BlancMatrix<number>::set (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k,
                                  const unsigned int l,
                                  const number       value)
{
    Assert (numbers::is_finite(value),
            ExcMessage("The given value is not finite but either "
                       "infinite or Not A Number (NaN)"));


    const unsigned int index = cols->operator()(i, j)*eltdim*eltdim + k*eltdim+l;

    //    std::cout << "i: " << i << "\tj: "<< j << "\tk: " <<k<<"\tl: " << l << "\t\t" << index << std::endl;


    //! it is allowed to set elements of
    //! the matrix that are not part of
    //! the sparsity pattern, if the
    //! value to which we set it is zero
    if (index == SparsityPattern::invalid_entry)
    {
        Assert ((cols->operator()(i, j) != SparsityPattern::invalid_entry) ||
                (value == 0.),
                ExcInvalidIndex(i, j));
        return;
    }

    this->val[index] = value;

}







template <typename number>
        inline
        void
        BlancMatrix<number>::add (const unsigned int i,
                                  const unsigned int j,
                                  const unsigned int k,
                                  const unsigned int l,
                                  const number       value)
{
    Assert (numbers::is_finite(value),
            ExcMessage("The given value is not finite but either "
                       "infinite or Not A Number (NaN)"));

    if (value == 0)
        return;

    const unsigned int index = cols->operator()(i, j)*eltdim*eltdim + k*eltdim+l;

    //! it is allowed to add elements to
    //! the matrix that are not part of
    //! the sparsity pattern, if the
    //! value to which we set it is zero
    if ( cols->operator()(i, j) == SparsityPattern::invalid_entry)
    {
        Assert ((index != SparsityPattern::invalid_entry) ||
                (value == 0.),
                ExcInvalidIndex(i, j));
        return;
    }



    val[index] += value;
}





//!funktioniert
template <typename number>
        inline
        BlancMatrix<number> &
        BlancMatrix<number>::operator *= (const number factor)
                                         {
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());

    number             *val_ptr    = &val[0];
    const number *const end_ptr    = &val[cols->n_nonzero_elements()*eltdim*eltdim];

    while (val_ptr != end_ptr)
        *val_ptr++ *= factor;

    return *this;
}


//!funktioniert
template <typename number>
        inline
        BlancMatrix<number> &
        BlancMatrix<number>::operator /= (const number factor)
                                         {
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert (factor !=0, ExcDivideByZero());

    const number factor_inv = 1. / factor;

    number             *val_ptr    = &val[0];
    const number *const end_ptr    = &val[cols->n_nonzero_elements()*eltdim*eltdim];

    while (val_ptr != end_ptr)
        *val_ptr++ *= factor_inv;

    return *this;
}

template <typename number>
        template <typename somenumber>
        void
        BlancMatrix<number>::add (const number factor,
                                  const BlancMatrix<somenumber> &matrix)
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

    number             *val_ptr    = &val[0];
    const somenumber   *matrix_ptr = &matrix.val[0];
    const number *const end_ptr    = &val[cols->n_nonzero_elements()*(eltdim*eltdim)];

    while (val_ptr != end_ptr)
        *val_ptr++ += factor * *matrix_ptr++;
}

//!funktioniert
template <typename number>
        inline
        number BlancMatrix<number>::operator () (const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const unsigned int l
                                                 ) const
{
    Assert (cols != 0, ExcNotInitialized());
    AssertThrow (cols->operator()(i,j) != SparsityPattern::invalid_entry,
                 ExcInvalidIndex(i,j));


    const unsigned int index = cols->operator()(i, j)*eltdim*eltdim + k*eltdim+l;


    return val[index];
}

template <typename number>
        inline
        number BlancMatrix<number>::global_entry (const unsigned int j) const
{

    //  Assert( eltdim != 1, ExcIndexRange(this->dimension(),1,1))
    Assert (cols != 0, ExcNotInitialized());
    Assert (j < cols->n_nonzero_elements()*(eltdim*eltdim),
            ExcIndexRange (j, 0, cols->n_nonzero_elements()*eltdim*eltdim));


    return val[j];
}


//!funktioniert
template <typename number>
        inline
        number BlancMatrix<number>::el (const unsigned int i,
                                        const unsigned int j,
                                        const unsigned int k,
                                        const unsigned int l ) const
{
    Assert (cols != 0, ExcNotInitialized());
    const unsigned int index = cols->operator()(i, j);

    if (index != SparsityPattern::invalid_entry){

        const unsigned int index = cols->operator()(i, j)*eltdim*eltdim + k*eltdim+l;



        return val[index];


    }else
        return 0;
}



//!---------------------------------------------------------------------------


template <typename number>
        BlancMatrix<number>::BlancMatrix ()
            :
            cols(0, "BlancMatrix"),
            val(0),
            inv_val(0),
            max_len(0),
            eltdim(0),
            calced_inverse(false)
{}



template <typename number>
        BlancMatrix<number>::BlancMatrix (const BlancMatrix &m)
            :
            Subscriptor (m),
            cols(0, "BlancMatrix"),
            val(0),
            inv_val(0),
            max_len(0),
            eltdim(0),
            calced_inverse(false)
{
    Assert (m.cols==0, ExcInvalidConstructorCall());
    Assert (m.val==0, ExcInvalidConstructorCall());
    Assert (m.inv_val==0, ExcInvalidConstructorCall());
    Assert (m.max_len==0, ExcInvalidConstructorCall());
    Assert (m.eltdim==0, ExcInvalidConstructorCall());
}



template <typename number>
        BlancMatrix<number>&
        BlancMatrix<number>::operator = (const BlancMatrix<number> &m)
                                        {
    Assert (m.cols==0, ExcInvalidConstructorCall());
    Assert (m.val==0, ExcInvalidConstructorCall());
    Assert (m.inv_val==0, ExcInvalidConstructorCall());
    Assert (m.max_len==0, ExcInvalidConstructorCall());
    Assert (m.eltdim==0, ExcInvalidConstructorCall());


    return *this;
}



template <typename number>
        BlancMatrix<number>::BlancMatrix (const SparsityPattern &c, const int eltdim=1)
            :
            cols(0, "BlancMatrix"),
            val(0),
            inv_val(0),
            max_len(0),
            eltdim(0),
            calced_inverse(false)
{
    reinit (c,eltdim);
}

template <typename number>
        BlancMatrix<number>::~BlancMatrix ()
{
    cols = 0;
    eltdim = 0;
    if (val != 0)
        delete[] val;

    if (calced_inverse)
        delete[] inv_val;
}

template <typename number>
        template <typename somenumber>
        BlancMatrix<number> &
        BlancMatrix<number>::copy_from (const BlancMatrix<somenumber> &matrix)
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

    std::copy (&matrix.val[0], &matrix.val[cols->n_nonzero_elements()*(eltdim*eltdim)],
               &val[0]);

    return *this;
}

template <typename number>
        template <typename somenumber>
        BlancMatrix<number> &
        BlancMatrix<number>::copy_from (const SparseMatrix<somenumber> &matrix)
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());


    for(int i = 0; i < m(); i++)
    {
        for(int j = 0; j < cols->row_length(i); j++)
        {
            number value = matrix.raw_entry(i,j);
            val[cols->get_rowstart_indices()[i]+j] = value;
        }
    }
    this->eltdim = 1;

    return *this;
}



template <typename number>
        BlancMatrix<number> &
        BlancMatrix<number>::operator = (const double d)
                                        {
    Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

    Assert (cols != 0, ExcNotInitialized());
    Assert (cols->compressed || cols->empty(), SparsityPattern::ExcNotCompressed());

    if (val)
        std::fill_n (&val[0], cols->n_nonzero_elements()*(eltdim*eltdim), d);

    return *this;
}

template <typename number>
        void
        BlancMatrix<number>::reinit (const SparsityPattern &sparsity, const int eltdim)
{
    cols = &sparsity;
    this->eltdim = eltdim;

    if (cols->empty())
    {
        if (val != 0)
            delete[] val;
        if (calced_inverse)
            delete[] inv_val;

        val = 0;
        inv_val = 0;
        max_len = 0;
        return;
    }

    const unsigned int N = cols->n_nonzero_elements() * (eltdim*eltdim);
    if (N > max_len || max_len == 0)
    {
        if (val != 0)
            delete[] val;
        val = new number[N];

        if (calced_inverse)
            delete[] inv_val;
        inv_val = new number[N];

        max_len = N;
    }

    if (val != 0)
        std::fill_n (&val[0], N, 0);

    if (calced_inverse)
        std::fill_n (&inv_val[0], N, 0);
}

template <typename number>
        void
        BlancMatrix<number>::clear ()
{
    cols = 0;
    eltdim = 0;
    if (val) delete[] val;
    if (inv_val) delete[] inv_val;
    val = 0;
    inv_val = 0;
    max_len = 0;
    calced_inverse=false;
}

template <typename number>
        bool
        BlancMatrix<number>::empty () const
{
    if (cols == 0)
        return true;
    else
        return cols->empty();
}



template <typename number>
        unsigned int
        BlancMatrix<number>::n_nonzero_elements () const
{
    Assert (cols != 0, ExcNotInitialized());
    return cols->n_nonzero_elements ();
}

template <typename number>
        const SparsityPattern &
        BlancMatrix<number>::get_sparsity_pattern () const
{
    Assert (cols != 0, ExcNotInitialized());
    return *cols;
}

template <typename number>
         unsigned int
         BlancMatrix<number>::memory_consumption () const
{
            return sizeof(*this) + max_len*sizeof(number);
}

template <typename number>
        void BlancMatrix<number>::print (std::ostream &out) const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());

    for (unsigned int i=0; i<cols->n_rows(); ++i)
        for (unsigned int j=cols->get_rowstart_indices()[i]; j<cols->get_rowstart_indices()[i+1]; ++j)
            for(unsigned int k=0; k<eltdim; ++k)
                for(unsigned int l=0; l<eltdim; ++l){
                    if(eltdim==1)
                        out << "(" << i << "," << cols->get_column_numbers()[j] << ") " << val[j] << std::endl;
                    else
                        out << "(" << i << "," << cols->get_column_numbers()[j] << ") (" << k<< "," << l<< ") "<< val[j*eltdim*eltdim+k*eltdim+l] << std::endl;
                }



//    for(int i = 0; i < this->n_nonzero_elements()*eltdim*eltdim; i++)
//        std::cout << i << "\t" << val[i] << std::endl;

    std::cout << std::endl;

    AssertThrow (out, ExcIO());
}



// man sollte ein Möglichkeit finden, dass auch Blockmatrizen ausgegeben werden können
template <typename number>
        void BlancMatrix<number>::print_formatted (std::ostream &out,
                                                   const unsigned int precision,
                                                   const bool scientific,
                                                   const unsigned int width_,
                                                   const char* zero_string,
                                                   const double denominator) const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());

//    if(eltdim != 1)
//    {
//        std::cout << __FUNCTION__ << std::endl;
//        std::cout << "this function works only useful with eltdim = 1" << std::endl;
//        return;
//    }
    unsigned int width = width_;

    std::ios::fmtflags old_flags = out.flags();
    unsigned int old_precision = out.precision (precision);

    if (scientific)
    {
        out.setf (std::ios::scientific, std::ios::floatfield);
        if (!width)
            width = precision+7;
    } else {
        out.setf (std::ios::fixed, std::ios::floatfield);
        if (!width)
            width = precision+2;
    }

    // oberer Deckel
    for (unsigned int i=0; i<eltdim*m()*(width+1)+n(); ++i)
        out << "-";
    out << std::endl;

    for (unsigned int i=0; i<m(); ++i)  {
        out << "|";

        for(unsigned int k=0; k<eltdim; ++k) {

            for (unsigned int j=0; j<n(); ++j) {

                for (unsigned int l=0; l<eltdim; ++l) {

                  if ((*cols)(i,j) != SparsityPattern::invalid_entry){
                    out << std::setw(width)
                        << val[cols->operator()(i,j)*(eltdim*eltdim)+k*eltdim+l] * denominator << ' ';
                  }else{
                    out << std::setw(width) << zero_string << ' ';
                  }
                }
                out << "|" ;
            }
            out << std::endl << "|";

        }


        for (unsigned int i=0; i<eltdim*m()*(width+1)+(n()-1); ++i)
            out << "-";
        out << "|" << std::endl;


      }

    // unterer Deckel
//    for (unsigned int i=0; i<m()*(width+1)+2; ++i)
//        out << "-";
//    out << std::endl;


    AssertThrow (out, ExcIO());

    //! reset output format
    out.precision(old_precision);
    out.flags (old_flags);
}


//!funktioniert
template <typename number>
        typename BlancMatrix<number>::real_type
        BlancMatrix<number>::l1_norm () const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());

    Vector<real_type> column_sums( n()*eltdim );
    const unsigned int n_rows = m();

    for (unsigned int i=0; i<n_rows; ++i){
        for (unsigned int j=cols->get_rowstart_indices()[i]; j<cols->get_rowstart_indices()[i+1]; ++j){

            for(unsigned int k=0; k < eltdim; ++k){
                for(unsigned int l=0; l < eltdim; ++l){
                    column_sums( cols->get_column_numbers()[j] ) += numbers::NumberTraits<number>::abs(val[j*eltdim*eltdim+k*eltdim+l]);
                }
            }
        }
    }



    return column_sums.linfty_norm();
}


//!funktioniert
template <typename number>
        typename BlancMatrix<number>::real_type
        BlancMatrix<number>::linfty_norm () const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());

    Vector<real_type> column_sums( n()*eltdim );
    const unsigned int n_rows = m();

    for (unsigned int i=0; i<n_rows; ++i){
         for (unsigned int j=cols->get_rowstart_indices()[i]; j<cols->get_rowstart_indices()[i+1]; ++j){
            for(unsigned int k=0; k < eltdim; ++k){
                for(unsigned int l=0; l < eltdim; ++l){
                    column_sums(i) += numbers::NumberTraits<number>::abs(val[j*eltdim*eltdim+k*eltdim+l]);
                }
            }
        }
    }

    return column_sums.linfty_norm();
}



template <typename number>
        typename BlancMatrix<number>::real_type
        BlancMatrix<number>::frobenius_norm () const
{
    //! simply add up all entries in the
    //! sparsity pattern, without taking any
    //! reference to rows or columns
    real_type norm_sqr = 0;
    const unsigned int n_rows = m();

    for (unsigned int i=0; i<n_rows; ++i){
         for (unsigned int j=cols->get_rowstart_indices()[i]; j<cols->get_rowstart_indices()[i+1]; ++j){

            for(unsigned int k=0; k < eltdim; ++k){
                for(unsigned int l=0; l < eltdim; ++l){
                    norm_sqr += numbers::NumberTraits<number>::abs_square(val[j*eltdim*eltdim+k*eltdim+l]);
                }
            }
        }
    }


    return std::sqrt (norm_sqr);
}
#warning "dies muss wieder weg"
//#define DEAL_II_USE_MT 1

template <typename number>
        template <class OutVector, class InVector>
        void
        BlancMatrix<number>::vmult (OutVector& dst,
                                    const InVector& src) const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert(m() == (dst.size() /eltdim), ExcDimensionMismatch(m(),dst.size()/eltdim));
    Assert(n() == (src.size() /eltdim), ExcDimensionMismatch(n(),src.size()/eltdim));


    const unsigned int n_rows = m();
    dst = 0.;

//    std::cout << "deal makro: " << DEAL_II_USE_MT << std::endl;
//    std::cout << "default threads: " << multithread_info.n_default_threads << std::endl;
//    std::cout << "ratio: " << n_rows/multithread_info.n_default_threads << std::endl;

    if ( DEAL_II_USE_MT &&
         (multithread_info.n_default_threads > 1) &&
         (n_rows/multithread_info.n_default_threads > 2000))
    {

        const unsigned int n_threads = multithread_info.n_default_threads;
//        std::cout << "MT on" << std::endl;

        //! then spawn threads. since
        //! some compilers have trouble
        //! finding out which
        //! 'encapsulate' function to
        //! take of all those possible
        //! ones if we simply drop in
        //! the address of an overloaded
        //! template member function,
        //! make it simpler for the
        //! compiler by giving it the
        //! correct type right away:
        typedef
                void (BlancMatrix<number>::*mem_fun_p)
                (OutVector &,
                 const InVector &,
                 const int,
                 const int
                 ) const;
        const mem_fun_p comp
                = (&BlancMatrix<number>::
                   template threaded_vmult<OutVector,InVector>);

        Threads::ThreadGroup<> threads;

        for (unsigned int i=0; i<n_threads; ++i)
        {

            int range_begin = i*n_rows/n_threads ;
            int range_end = (i+1)*n_rows/n_threads ;

            threads += Threads::spawn<void, BlancMatrix<number> > (*this, comp) (dst, src,
                                                                                 range_begin, range_end
                                                                                 );
        }

        threads.join_all();
        return;
    }
    else
    {

        if(eltdim == 1)
        {


            for(int i = 0; i < m(); i++)
            {
                number buffer = 0.0;
                int rowstart = cols->get_rowstart_indices()[i];

                for(int j = 0; j < cols->row_length(i); j++)
                {

                    buffer += this->val[rowstart+j] * src( cols->column_number(i,j) );
                }
                dst(i) = buffer;
            }




        }
        else
        {
            int k = 0;
            int l = 0;
            int p = 0;
            int q = 0;



            for (int i=0,s=0,t=eltdim ; i<m() ; i++,s=t,t+=eltdim )
            {

                int rowstart = cols->get_rowstart_indices()[i]*eltdim*eltdim;



                for (int j=0,r=0 ; j<cols->row_length(i) ; j++ ) {

                    p = cols->column_number(i,j) * eltdim;
                    q = p + eltdim;

                    for ( k=s ; k<t ; k++ ) {

                        number buffer = 0.0;

                        for ( l=p ; l<q ; l++,r++ ) {

                            buffer += this->val[rowstart+r] * src(l);
                        } /*l*/

                        dst(k) += buffer;

                    } /*k*/
                } /*j*/
            } /*i*/

        }
    }
}


template <typename number>
        template <class OutVector, class InVector>
        void
        BlancMatrix<number>::threaded_vmult (OutVector       &dst,
                                             const InVector &src,
                                             const int        range_begin,
                                             const int        range_end
                                             ) const
{
    //! this function should not be called
    //! when not in parallel mode.
    Assert (DEAL_II_USE_MT, ExcInternalError());


    if(eltdim == 1)
    {
        for(int i = range_begin; i < range_end; i++)
        {
            number buffer = 0.0;

            for(int j = 0; j < cols->row_length(i); j++)
            {
                int rowstart = cols->get_rowstart_indices()[i];
                buffer += this->val[rowstart+j] * src( cols->column_number(i,j) );
            }
            dst(i) = buffer;
        }

    }
    else
    {

        int k = 0;
        int l = 0;
        int p = 0;
        int q = 0;



        for (int i=range_begin,s=i*eltdim,t=(i+1)*eltdim ; i<range_end ; i++,s=t,t+=eltdim )
        {

            int rowstart = cols->get_rowstart_indices()[i]*eltdim*eltdim;

            for (int j=0,r=0 ; j<cols->row_length(i) ; j++ ) {

                p = cols->column_number(i,j) * eltdim;


                q = p + eltdim;


                for ( k=s ; k<t ; k++ ) {

                    number buffer = 0.0;

                    for ( l=p ; l<q ; l++,r++ ) {

                        buffer += this->val[rowstart+r] * src(l);
                    } /*l*/

                    dst(k) += buffer;

                } /*k*/
            } /*j*/
        } /*i*/

    }

}




template <typename number>
        template <typename somenumber>
        void
        BlancMatrix<number>::residual (Vector<somenumber>       &dst,
                                       const Vector<somenumber> &u,
                                       const Vector<somenumber> &b) const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert(m() == dst.size()/eltdim, ExcDimensionMismatch(m(),dst.size()/eltdim));
    Assert(m() == b.size()/eltdim, ExcDimensionMismatch(m(),b.size()/eltdim));
    Assert(n() == u.size()/eltdim, ExcDimensionMismatch(n(),u.size()/eltdim));

    Assert (&u != &dst, ExcSourceEqualsDestination());


    const unsigned int n_rows = m();

    if(eltdim==1)
        dst=0;
    else
        dst=b;


    if ( DEAL_II_USE_MT &&
         (multithread_info.n_default_threads > 1) &&
         (n_rows/multithread_info.n_default_threads > 2000))
    {

        const unsigned int n_threads = multithread_info.n_default_threads;

        //!std::cout << "MT on" << std::endl;

        //! then spawn threads. since
        //! some compilers have trouble
        //! finding out which
        //! 'encapsulate' function to
        //! take of all those possible
        //! ones if we simply drop in
        //! the address of an overloaded
        //! template member function,
        //! make it simpler for the
        //! compiler by giving it the
        //! correct type right away:


        typedef
                void (BlancMatrix<number>::*mem_fun_p)
                (Vector<somenumber>       &,
                 const Vector<somenumber> &,
                 const Vector<somenumber> &,
                 const int,
                 const int ) const;

        const mem_fun_p comp_residual = &BlancMatrix<number>::
                                        template threaded_residual<somenumber>;





        Threads::ThreadGroup<> threads;

        for (unsigned int i=0; i<n_threads; ++i)
        {

            int range_begin = i*n_rows/n_threads ;
            int range_end = (i+1)*n_rows/n_threads ;

            threads += Threads::spawn<void, BlancMatrix<number> > (*this, comp_residual) (dst, u,b,
                                                                                          range_begin, range_end
                                                                                          );
        }

        threads.join_all();
        return;
    }
    else
    {


        if(eltdim == 1)
        {

            for(int i = 0; i < m(); i++)
            {
                number buffer = b(i);
                for(int j = 0; j < cols->row_length(i); j++)
                {

                    int rowstart = cols->get_rowstart_indices()[i];
                    buffer -= this->val[rowstart+j] * u( cols->column_number(i,j) );
                }
                dst(i) = buffer;
            }
        }
        else
        {
            int k = 0;
            int l = 0;
            int p = 0;
            int q = 0;



            for (int i=0,s=0,t=eltdim ; i<m() ; i++,s=t,t+=eltdim )
            {

                int rowstart = cols->get_rowstart_indices()[i]*eltdim*eltdim;

                const number* m_values_row = &this->val[rowstart];
                //!            const unsigned int* colidxs_row     = &cols->get_column_numbers()[i];


                for (int j=0,r=0 ; j<cols->row_length(i) ; j++ ) {

                    p = cols->column_number(i,j) * eltdim;
                    q = p + eltdim;

                    for ( k=s ; k<t ; k++ ) {

                        number buffer = 0;

                        for ( l=p ; l<q ; l++,r++ ) {

                            buffer += m_values_row[r] * u(l);
                        } /*l*/

                        dst(k) -= buffer;

                    } /*k*/
                } /*j*/
            } /*i*/
        }
    }
}



template <typename number>
        template <typename somenumber>
        void
        BlancMatrix<number>::threaded_residual (Vector<somenumber>       &dst,
                                                const Vector<somenumber> &u,
                                                const Vector<somenumber> &b,
                                                const int        range_begin,
                                                const int        range_end
                                                ) const
{
    //! this function should not be called
    //! when not in parallel mode.
    Assert (DEAL_II_USE_MT, ExcInternalError());


    if(eltdim == 1)
    {

        for(int i = range_begin; i < range_end; i++)
        {
            number buffer = b(i);
            for(int j = 0; j < cols->row_length(i); j++)
            {

                int rowstart = cols->get_rowstart_indices()[i];
                buffer -= this->val[rowstart+j] * u( cols->column_number(i,j) );
            }
            dst(i) = buffer;
        }

    }
    else
    {

        int k = 0;
        int l = 0;
        int p = 0;
        int q = 0;



        for (int i=range_begin,s=(i)*eltdim,t=(i+1)*eltdim ; i<range_end ; i++,s=t,t+=eltdim )
        {

            int rowstart = cols->get_rowstart_indices()[i]*eltdim*eltdim;

            const number* m_values_row = &this->val[rowstart];
            //!            const unsigned int* colidxs_row     = &cols->get_column_numbers()[i];


            for (int j=0,r=0 ; j<cols->row_length(i) ; j++ ) {

                p = cols->column_number(i,j) * eltdim;
                q = p + eltdim;

                for ( k=s ; k<t ; k++ ) {

                    number buffer = 0;

                    for ( l=p ; l<q ; l++,r++ ) {

                        buffer += m_values_row[r] * u(l);
                    } /*l*/

                    dst(k) -= buffer;

                } /*k*/
            } /*j*/
        } /*i*/

    }

}

template <typename number>
        void BlancMatrix<number>::calc_inverse() const
{

    if(calced_inverse == true) return;

    number* diagelt_values = NULL;
    number* inv_diag_val = NULL;

    std::cout << "calculate the diag inv" << std::endl;

    dealii::FullMatrix<number> diag_element(eltdim);
    dealii::FullMatrix<number> inv_element(eltdim);

#warning multithreading wie bei vmult auch

    for ( int i = 0; i < m(); i++ ) {
        diagelt_values = & this->val    [ cols->get_rowstart_indices()[i] * (eltdim*eltdim)];
        inv_diag_val   = & this->inv_val[ cols->get_rowstart_indices()[i] * (eltdim*eltdim)];

        //! Diagonaleintraege umkopieren
        for(int x=0; x<eltdim; x++){
            for(int y=0; y<eltdim; y++){
                diag_element(x,y) = diagelt_values[x*eltdim+y];
            }
        }


        //! Das Diagonalelement invertieren
        inv_element.invert(diag_element);


        //! Das Inverse an die richtige Stelle in der inversen Matrix setzen
        for(int x=0; x<eltdim; x++){
            for(int y=0; y<eltdim; y++){
                inv_diag_val[x*eltdim+y] = inv_element(x,y);
            }
        }






    } /*i*/

    calced_inverse = true;

}

template <typename number>
        void BlancMatrix<number>::clear_inverse()
{
    if(calced_inverse == false) return;


    this->calced_inverse = false;

}


template <typename number>
        template <typename somenumber>
        void
        BlancMatrix<number>::precondition_Jacobi (Vector<somenumber>       &dst,
                                                  const Vector<somenumber> &src,
                                                  const number              om) const
{
    Assert (cols != 0, ExcNotInitialized());
    Assert (val != 0, ExcNotInitialized());
    Assert (cols->optimize_diagonal(),
            typename SparsityPattern::ExcDiagonalNotOptimized());

    Assert (dst.size()/eltdim == n(), ExcDimensionMismatch (dst.size()/eltdim, n()));
    Assert (src.size()/eltdim == n(), ExcDimensionMismatch (src.size()/eltdim, n()));

    const unsigned int n = src.size()/eltdim;


    if(eltdim == 1)
    {
        if(om != 1. )
            for(int i = 0; i < n; i++)
                dst(i) = om* src(i)/ this->val[ cols->get_rowstart_indices()[i] ];

        else
            for(int i = 0; i < n; i++)
                dst(i) = src(i)/ this->val[ cols->get_rowstart_indices()[i] ];
    } else {



        if(calced_inverse!=true)    calc_inverse();


        number* diagelt_values = NULL;


#warning multithreading wie bei vmult auch

        int i,s,t,k,r,l;

        for ( i=0,s=0,t=eltdim ; i< n ; i++,s=t,t+=eltdim ) {
            diagelt_values = & this->inv_val[ cols->get_rowstart_indices()[i] * (eltdim*eltdim)];


            for ( k=s,r=0 ; k<t ; k++ ) {
                number buffer = 0.0;
                for ( l=s ; l<t ; l++,r++ ) {
                    buffer +=  om * diagelt_values[r] * src(l);

                } /*l*/
                dst(k) =  buffer;

            } /*k*/
        } /*i*/

    }


}

#endif
