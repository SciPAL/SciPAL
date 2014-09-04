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


//----------------------------  solver_gmres.h  ---------------------------
//    $Id: solver_gmres.h 23876 2011-06-28 18:21:51Z kanschat $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_gmres.h  ---------------------------
#ifndef SOLVER_GMRES_H
#define SOLVER_GMRES_H


#include<QDebug>
#include <float.h>

#include <lac/solver_cuda.h>

//#include <deal.II/base/config.h>
//#include <deal.II/base/subscriptor.h>
//#include <deal.II/base/logstream.h>
#include <deal.II/lac/householder.h>
//#include <deal.II/lac/solver.h>
//#include <deal.II/lac/solver_control.h>
//#include <deal.II/lac/full_matrix.h>
//#include <deal.II/lac/vector.h>

#include <vector>
#include <cmath>


/*!@addtogroup Solvers */
/*@{*/

namespace internal
{
				   /**
				    * A namespace for a helper class
				    * to the GMRES solver.
				    */
  namespace SolverGMRES
  {
				     /**
				      * Class to hold temporary
				      * vectors.  This class
				      * automatically allocates a new
				      * vector, once it is needed.
				      *
				      * A future version should also
				      * be able to shift through
				      * vectors automatically,
				      * avoiding restart.
				      */

    template <class VECTOR>
    class TmpVectors
    {
      public:
					 /**
					  * Constructor. Prepares an
					  * array of @p VECTOR of
					  * length @p max_size.
					  */
	TmpVectors(const unsigned int    max_size,
           dealii::VectorMemory<VECTOR> &vmem);

                                         /**
					  * Delete all allocated vectors.
					  */
	~TmpVectors();

					 /**
					  * Get vector number
					  * @p i. If this vector was
					  * unused before, an error
					  * occurs.
					  */
	VECTOR& operator[] (const unsigned int i) const;

					 /**
					  * Get vector number
					  * @p i. Allocate it if
					  * necessary.
					  *
					  * If a vector must be
					  * allocated, @p temp is
					  * used to reinit it to the
					  * proper dimensions.
					  */
	VECTOR& operator() (const unsigned int i,
			    const VECTOR      &temp);

      private:
					 /**
					  * Pool were vectors are
					  * obtained from.
					  */
    dealii::VectorMemory<VECTOR> &mem;

					 /**
					  * Field for storing the
					  * vectors.
					  */
	std::vector<VECTOR*> data;

					 /**
					  * Offset of the first
					  * vector. This is for later
					  * when vector rotation will
					  * be implemented.
					  */
	unsigned int offset;
    };
  }
}

/**
 * Implementation of the Restarted Preconditioned Direct Generalized
 * Minimal Residual Method. The stopping criterion is the norm of the
 * residual.
 *
 * The AdditionalData structure contains the number of temporary
 * vectors used. The size of the Arnoldi basis is this number minus
 * three. Additionally, it allows you to choose between right or left
 * preconditioning. The default is left preconditioning. Finally it
 * includes a flag indicating whether or not the default residual is
 * used as stopping criterion.

 * <h3>Left versus right preconditioning</h3>
 *
 * @p AdditionalData allows you to choose between left and right
 * preconditioning. As expected, this switches between solving for the
 * systems <i>P<sup>-1</sup>A</i> and <i>AP<sup>-1</sup></i>,
 * respectively.
 *
 * A second consequence is the type of residual which is used to
 * measure convergence. With left preconditioning, this is the
 * <b>preconditioned</b> residual, while with right preconditioning,
 * it is the residual of the unpreconditioned system.
 *
 * Optionally, this behavior can be overridden by using the flag
 * AdditionalData::use_default_residual. A <tt>true</tt> value refers
 * to the behavior described in the previous paragraph, while
 * <tt>false</tt> reverts it. Be aware though that additional
 * residuals have to be computed in this case, impeding the overall
 * performance of the solver.
 *
 * <h3>The size of the Arnoldi basis</h3>
 *
 * The maximal basis size is controlled by
 * AdditionalData::max_n_tmp_vectors, and it is this number minus 2.
 * If the number of iteration steps exceeds this number, all basis
 * vectors are discarded and the iteration starts anew from the
 * approximation obtained so far.
 *
 * Note that the minimizing property of GMRes only pertains to the
 * Krylov space spanned by the Arnoldi basis. Therefore, restarted
 * GMRes is <b>not</b> minimizing anymore. The choice of the basis
 * length is a trade-off between memory consumption and convergence
 * speed, since a longer basis means minimization over a larger
 * space.
 *
 * For the requirements on matrices and vectors in order to work with
 * this class, see the documentation of the Solver base class.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann.
 */
template <class VECTOR = dealii::Vector<double> >
class SolverGMRES : public Solver<VECTOR>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver.
				      */
    struct AdditionalData
    {
					 /**
					  * Constructor. By default, set the
					  * number of temporary vectors to 30,
					  * i.e. do a restart every
					  * 28 iterations. Also
					  * set preconditioning from left and
					  * the residual of the stopping
					  * criterion to the default residual.
					  */
	AdditionalData (const unsigned int max_n_tmp_vectors = 30,
                        const bool right_preconditioning = false,
                        const bool use_default_residual = true);

					 /**
					  * Maximum number of
					  * temporary vectors. This
					  * parameter controls the
					  * size of the Arnoldi basis,
					  * which for historical
					  * reasons is
					  * #max_n_tmp_vectors-2.
					  */
	unsigned int    max_n_tmp_vectors;

					 /**
					  * Flag for right
					  * preconditioning.
					  *
					  * @note Change between left
					  * and right preconditioning
					  * will also change the way
					  * residuals are
					  * evaluated. See the
					  * corresponding section in
					  * the SolverGMRES.
					  */
	bool right_preconditioning;

					 /**
					  * Flag for the default
					  * residual that is used to
					  * measure convergence.
					  */
	bool use_default_residual;
    };

				     /**
				      * Constructor.
				      */
    SolverGMRES (dealii::SolverControl        &cn,
         dealii::VectorMemory<VECTOR> &mem,
		 const AdditionalData &data=AdditionalData());

				     /**
				      * Constructor. Use an object of
				      * type GrowingVectorMemory as
				      * a default to allocate memory.
				      */
    SolverGMRES (dealii::SolverControl        &cn,
		 const AdditionalData &data=AdditionalData());

				     /**
				      * Solve the linear system $Ax=b$
				      * for x.
				      */
    template<class MATRIX, class PRECONDITIONER>
    void
    solve (const MATRIX         &A,
	   VECTOR               &x,
	   const VECTOR         &b,
	   const PRECONDITIONER &precondition);

    DeclException1 (ExcTooFewTmpVectors,
		    int,
		    << "The number of temporary vectors you gave ("
		    << arg1 << ") is too small. It should be at least 10 for "
		    << "any results, and much more for reasonable ones.");

  protected:
				     /**
				      * Includes the maximum number of
				      * tmp vectors.
				      */
    AdditionalData additional_data;

				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();

				     /**
				      * Transformation of an upper
				      * Hessenberg matrix into
				      * tridiagonal structure by givens
				      * rotation of the last column
				      */
    void givens_rotation (dealii::Vector<double> &h,  dealii::Vector<double> &b,
              dealii::Vector<double> &ci, dealii::Vector<double> &si,
			  int col) const;
				     /**
				      * Projected system matrix
				      */
    dealii::FullMatrix<double> H;
				     /**
				      * Auxiliary matrix for inverting @p H
				      */
    dealii::FullMatrix<double> H1;

  private:
    				     /**
				      * No copy constructor.
				      */
    SolverGMRES (const SolverGMRES<VECTOR>&);
};

/**
 * Implementation of the Generalized minimal residual method with flexible
 * preconditioning method.
 *
 * This version of the GMRES method allows for the use of a different
 * preconditioner in each iteration step. Therefore, it is also more
 * robust with respect to inaccurate evaluation of the
 * preconditioner. An important application is also the use of a
 * Krylov space method inside the preconditioner.
 *
 * FGMRES needs two vectors in each iteration steps yielding a total
 * of <tt>2 * SolverFGMRESAdditionalData::max_basis_size+1</tt>
 * auxiliary vectors.
 *
 * Caveat: documentation of this class is not up to date. There are
 * also a few parameters of GMRES we would like to introduce here.
 *
 * @author Guido Kanschat, 2003
 */
template <class VECTOR = dealii::Vector<double> >
class SolverFGMRES : public Solver<VECTOR>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver.
				      */
    struct AdditionalData
    {
					 /**
					  * Constructor. By default,
					  * set the number of
					  * temporary vectors to 30,
					  * preconditioning from left
					  * and the residual of the
					  * stopping criterion to the
					  * default residual
					  * (cf. class documentation).
					  */
	AdditionalData(const unsigned int max_basis_size = 30,
		       const bool /*use_default_residual*/ = true)
			:
			max_basis_size(max_basis_size)
	  {}

					 /**
					  * Maximum number of
					  * tmp vectors.
					  */
	unsigned int    max_basis_size;
    };

				     /**
				      * Constructor.
				      */
    SolverFGMRES (dealii::SolverControl        &cn,
          dealii::VectorMemory<VECTOR> &mem,
		  const AdditionalData &data=AdditionalData());

				     /**
				      * Constructor. Use an object of
				      * type GrowingVectorMemory as
				      * a default to allocate memory.
				      */
    SolverFGMRES (dealii::SolverControl        &cn,
		  const AdditionalData &data=AdditionalData());

				     /**
				      * Solve the linear system $Ax=b$
				      * for x.
				      */
    template<class MATRIX, class PRECONDITIONER>
    void
    solve (const MATRIX         &A,
	   VECTOR               &x,
	   const VECTOR         &b,
	   const PRECONDITIONER &precondition);

  private:
                                     /**
				      * Additional flags.
				      */
    AdditionalData additional_data;
				     /**
				      * Projected system matrix
				      */
    dealii::FullMatrix<double> H;
				     /**
				      * Auxiliary matrix for inverting @p H
				      */
    dealii::FullMatrix<double> H1;
};

/*@}*/
/* --------------------- Inline and template functions ------------------- */


#ifndef DOXYGEN
namespace internal
{
  namespace SolverGMRES
  {
    template <class VECTOR>
    inline
    TmpVectors<VECTOR>::
    TmpVectors (const unsigned int    max_size,
        dealii::VectorMemory<VECTOR> &vmem)
		    :
		    mem(vmem),
		    data (max_size, 0),
		    offset(0)
    {}


    template <class VECTOR>
    inline
    TmpVectors<VECTOR>::~TmpVectors ()
    {
      for (typename std::vector<VECTOR*>::iterator v = data.begin();
	   v != data.end(); ++v)
	if (*v != 0)
	  mem.free(*v);
    }


    template <class VECTOR>
    inline VECTOR&
    TmpVectors<VECTOR>::operator[] (const unsigned int i) const
    {
      Assert (i+offset<data.size(),
	      ExcIndexRange(i, -offset, data.size()-offset));

      Assert (data[i-offset] != 0, ExcNotInitialized());
      return *data[i-offset];
    }


    template <class VECTOR>
    inline VECTOR&
    TmpVectors<VECTOR>::operator() (const unsigned int i,
				    const VECTOR      &temp)
    {
      Assert (i+offset<data.size(),
	      ExcIndexRange(i,-offset, data.size()-offset));
      if (data[i-offset] == 0)
	{
	  data[i-offset] = mem.alloc();
	  data[i-offset]->reinit(temp);
	}
      return *data[i-offset];
    }
  }
}



template <class VECTOR>
inline
SolverGMRES<VECTOR>::AdditionalData::
AdditionalData (const unsigned int max_n_tmp_vectors,
                const bool         right_preconditioning,
                const bool         use_default_residual)
                :
                max_n_tmp_vectors(max_n_tmp_vectors),
                right_preconditioning(right_preconditioning),
                use_default_residual(use_default_residual)
{}


template <class VECTOR>
SolverGMRES<VECTOR>::SolverGMRES (dealii::SolverControl        &cn,
                                  dealii::VectorMemory<VECTOR> &mem,
                                  const AdditionalData &data)
		:
		Solver<VECTOR> (cn,mem),
		additional_data(data)
{}



template <class VECTOR>
SolverGMRES<VECTOR>::SolverGMRES (dealii::SolverControl        &cn,
                                  const AdditionalData &data) :
		Solver<VECTOR> (cn),
		additional_data(data)
{}



template <class VECTOR>
inline
void
SolverGMRES<VECTOR>::givens_rotation (dealii::Vector<double> &h,
                      dealii::Vector<double> &b,
                      dealii::Vector<double> &ci,
                      dealii::Vector<double> &si,
				      int     col) const
{
  for (int i=0 ; i<col ; i++)
    {
      const double s = si(i);
      const double c = ci(i);
      const double dummy = h(i);
      h(i)   =  c*dummy + s*h(i+1);
      h(i+1) = -s*dummy + c*h(i+1);
    };

  const double r = 1./std::sqrt(h(col)*h(col) + h(col+1)*h(col+1));
  si(col) = h(col+1) *r;
  ci(col) = h(col)   *r;
  h(col)  =  ci(col)*h(col) + si(col)*h(col+1);
  b(col+1)= -si(col)*b(col);
  b(col) *=  ci(col);
}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverGMRES<VECTOR>::solve (const MATRIX         &A,
			    VECTOR               &x,
			    const VECTOR         &b,
			    const PRECONDITIONER &precondition)
{
				   // this code was written a very
				   // long time ago by people not
				   // associated with deal.II. we
				   // don't make any guarantees to its
				   // optimality or that it even works
				   // as expected...

//TODO:[?] Check, why there are two different start residuals.
//TODO:[GK] Make sure the parameter in the constructor means maximum basis size

  dealii::deallog.push("GMRES");
  const unsigned int n_tmp_vectors = additional_data.max_n_tmp_vectors;

				   // Generate an object where basis
				   // vectors are stored.
  ::internal::SolverGMRES::TmpVectors<VECTOR> tmp_vectors (n_tmp_vectors, this->memory);

				   // number of the present iteration; this
				   // number is not reset to zero upon a
				   // restart
  unsigned int accumulated_iterations = 0;

				   // matrix used for the orthogonalization
				   // process later
  H.reinit(n_tmp_vectors, n_tmp_vectors-1);

				   // some additional vectors, also used
				   // in the orthogonalization
  dealii::Vector<double>
    gamma(n_tmp_vectors),
    ci   (n_tmp_vectors-1),
    si   (n_tmp_vectors-1),
    h    (n_tmp_vectors-1);


  unsigned int dim = 0;

  dealii::SolverControl::State iteration_state = dealii::SolverControl::iterate;

				   // switch to determine whether we want a
				   // left or a right preconditioner. at
				   // present, left is default, but both
				   // ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;
				   // Per default the left
				   // preconditioned GMRes uses the
				   // preconditioned residual and the
				   // right preconditioned GMRes uses
				   // the unpreconditioned residual as
				   // stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

				   // define two aliases
  VECTOR &v = tmp_vectors(0, x);
  VECTOR &p = tmp_vectors(n_tmp_vectors-1, x);

				   // Following vectors are needed
				   // when not the default residuals
				   // are used as stopping criterion
  VECTOR *r=0;
  VECTOR *x_=0;
  dealii::Vector<double> *gamma_=0;
  if (!use_default_residual)
    {
      r=this->memory.alloc();
      x_=this->memory.alloc();
      r->reinit(x);
      x_->reinit(x);

      gamma_ = new dealii::Vector<double> (gamma.size());
    }

                                   ///////////////////////////////////
				   // outer iteration: loop until we
				   // either reach convergence or the
				   // maximum number of iterations is
				   // exceeded. each cycle of this
				   // loop amounts to one restart
  do
    {
				       // reset this vector to the
				       // right size
      h.reinit (n_tmp_vectors-1);

      if (left_precondition)
	{
	  A.vmult(p,x);
	  p.sadd(-1.,1.,b);
	  precondition.vmult(v,p);
	}
      else
	{
	  A.vmult(v,x);
	  v.sadd(-1.,1.,b);
	};

      double rho = v.l2_norm();

				       // check the residual here as
				       // well since it may be that we
				       // got the exact (or an almost
				       // exact) solution vector at
				       // the outset. if we wouldn't
				       // check here, the next scaling
				       // operation would produce
				       // garbage
      if (use_default_residual)
	{
	  iteration_state = this->control().check (
	    accumulated_iterations, rho);

      if (iteration_state != dealii::SolverControl::iterate)
	    break;
	}
      else
	{
      dealii::deallog << "default_res=" << rho << std::endl;

	  if (left_precondition)
	    {
	      A.vmult(*r,x);
	      r->sadd(-1.,1.,b);
	    }
	  else
	    precondition.vmult(*r,v);

	  double res = r->l2_norm();
	  iteration_state = this->control().check (
	    accumulated_iterations, res);

      if (iteration_state != dealii::SolverControl::iterate)
	    {
	      this->memory.free(r);
	      this->memory.free(x_);

	      delete gamma_;
	      break;
	    }
	}

      gamma(0) = rho;

      v *= 1./rho;

				       // inner iteration doing at
				       // most as many steps as there
				       // are temporary vectors. the
				       // number of steps actually
				       // been done is propagated
				       // outside through the @p dim
				       // variable
      for (unsigned int inner_iteration=0;
	   ((inner_iteration < n_tmp_vectors-2)
	    &&
        (iteration_state==dealii::SolverControl::iterate));
	   ++inner_iteration)
	{
	  ++accumulated_iterations;
	  // yet another alias
	  VECTOR& vv = tmp_vectors(inner_iteration+1, x);

	  if (left_precondition)
	    {
	      A.vmult(p, tmp_vectors[inner_iteration]);
	      precondition.vmult(vv,p);
	    } else {
	      precondition.vmult(p, tmp_vectors[inner_iteration]);
	      A.vmult(vv,p);
	    };

	  dim = inner_iteration+1;

					   /* Orthogonalization */
	  for (unsigned int i=0 ; i<dim ; ++i)
	    {
	      h(i) = vv * tmp_vectors[i];
	      vv.add(-h(i), tmp_vectors[i]);
	    };

	  				   /* Re-orthogonalization */
	  for (unsigned int i=0 ; i<dim ; ++i)
	    {
	      double htmp = vv * tmp_vectors[i];
	      h(i) += htmp;
	      vv.add(-htmp, tmp_vectors[i]);
	    }

	  const double s = vv.l2_norm();
	  h(inner_iteration+1) = s;
//TODO: s=0 is a lucky breakdown. Handle this somehow decently

	  vv *= 1./s;

					   /*  Transformation into
					       triagonal structure  */
	  givens_rotation(h,gamma,ci,si,inner_iteration);

					   /*  append vector on matrix  */
	  for (unsigned int i=0; i<dim; ++i)
	    H(i,inner_iteration) = h(i);

					   /*  default residual  */
	  rho = std::fabs(gamma(dim));

	  if (use_default_residual)
	    iteration_state = this->control().check (
	      accumulated_iterations, rho);
  	  else
  	    {
          dealii::deallog << "default_res=" << rho << std::endl;

          dealii::Vector<double> h_(dim);
	      *x_=x;
	      *gamma_=gamma;
	      H1.reinit(dim+1,dim);

	      for (unsigned int i=0; i<dim+1; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  H1(i,j) = H(i,j);

	      H1.backward(h_,*gamma_);

	      if (left_precondition)
		for (unsigned int i=0 ; i<dim; ++i)
		  x_->add(h_(i), tmp_vectors[i]);
	      else
		{
		  p = 0.;
		  for (unsigned int i=0; i<dim; ++i)
		    p.add(h_(i), tmp_vectors[i]);
		  precondition.vmult(*r,p);
		  x_->add(1.,*r);
		};
	      A.vmult(*r,*x_);
	      r->sadd(-1.,1.,b);
					       // Now *r contains the
					       // unpreconditioned
					       // residual!!
	      if (left_precondition)
		{
		  const double res=r->l2_norm();

		  iteration_state = this->control().check (
		    accumulated_iterations, res);
		}
	      else
		{
		  precondition.vmult(*x_, *r);
		  const double preconditioned_res=x_->l2_norm();

		  iteration_state = this->control().check (
		    accumulated_iterations, preconditioned_res);
		}
  	    }
	};
				       // end of inner iteration. now
				       // calculate the solution from
				       // the temporary vectors
      h.reinit(dim);
      H1.reinit(dim+1,dim);

      for (unsigned int i=0; i<dim+1; ++i)
	for (unsigned int j=0; j<dim; ++j)
	  H1(i,j) = H(i,j);

      H1.backward(h,gamma);

      if (left_precondition)
	for (unsigned int i=0 ; i<dim; ++i)
	  x.add(h(i), tmp_vectors[i]);
      else
	{
	  p = 0.;
	  for (unsigned int i=0; i<dim; ++i)
	    p.add(h(i), tmp_vectors[i]);
	  precondition.vmult(v,p);
	  x.add(1.,v);
	};
				       // end of outer iteration. restart if
				       // no convergence and the number of
				       // iterations is not exceeded
    }
  while (iteration_state == dealii::SolverControl::iterate);

  if (!use_default_residual)
    {
      this->memory.free(r);
      this->memory.free(x_);

      delete gamma_;
    }

  dealii::deallog.pop();
				   // in case of failure: throw
				   // exception
//  if (this->control().last_check() != SolverControl::success)
//    throw dealii::SolverControl::NoConvergence (this->control().last_step(),
//					this->control().last_value());
				   // otherwise exit as normal
}



template<class VECTOR>
double
SolverGMRES<VECTOR>::criterion ()
{
				   // dummy implementation. this function is
				   // not needed for the present implementation
				   // of gmres
  Assert (false, ExcInternalError());
  return 0;
}


//----------------------------------------------------------------------//

template <class VECTOR>
SolverFGMRES<VECTOR>::SolverFGMRES (dealii::SolverControl        &cn,
                    dealii::VectorMemory<VECTOR> &mem,
				    const AdditionalData &data)
		:
		Solver<VECTOR> (cn, mem),
		additional_data(data)
{}



template <class VECTOR>
SolverFGMRES<VECTOR>::SolverFGMRES (dealii::SolverControl        &cn,
				    const AdditionalData &data)
		:
		Solver<VECTOR> (cn),
		additional_data(data)
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverFGMRES<VECTOR>::solve (
  const MATRIX& A,
  VECTOR& x,
  const VECTOR& b,
  const PRECONDITIONER& precondition)
{
  dealii::deallog.push("FGMRES");

  dealii::SolverControl::State iteration_state = dealii::SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

				   // Generate an object where basis
				   // vectors are stored.
  typename ::internal::SolverGMRES::TmpVectors<VECTOR> v (basis_size, this->memory);
  typename ::internal::SolverGMRES::TmpVectors<VECTOR> z (basis_size, this->memory);

				   // number of the present iteration; this
				   // number is not reset to zero upon a
				   // restart
  unsigned int accumulated_iterations = 0;

				   // matrix used for the orthogonalization
				   // process later
  H.reinit(basis_size+1, basis_size);

				   // Vectors for projected system
  dealii::Vector<double> projected_rhs;
  dealii::Vector<double> y;

  // Iteration starts here

  VECTOR* aux = this->memory.alloc();
  aux->reinit(x);
  do
    {
      A.vmult(*aux, x);
      aux->sadd(-1., 1., b);

      double beta = aux->l2_norm();
      if (this->control().check(accumulated_iterations,beta)
      == dealii::SolverControl::success)
	break;

      H.reinit(basis_size+1, basis_size);
      double a = beta;

      for (unsigned int j=0;j<basis_size;++j)
	{
	  v(j,x).equ(1./a, *aux);

	  precondition.vmult(z(j,x), v[j]);
	  A.vmult(*aux, z[j]);

					   // Gram-Schmidt
	  for (unsigned int i=0;i<=j;++i)
	    {
	      H(i,j) = *aux * v[i];
	      aux->add(-H(i,j), v[i]);
	    }
	  H(j+1,j) = a = aux->l2_norm();

					   // Compute projected solution

	  if (j>0)
	    {
	      H1.reinit(j+1,j);
	      projected_rhs.reinit(j+1);
	      y.reinit(j);
	      projected_rhs(0) = beta;
	      H1.fill(H);
          dealii::Householder<double> house(H1);
	      double res = house.least_squares(y, projected_rhs);
	      iteration_state = this->control().check(++accumulated_iterations, res);
          if (iteration_state != dealii::SolverControl::iterate)
		break;
	    }
	}
				       // Update solution vector
      for (unsigned int j=0;j<y.size();++j)
	x.add(y(j), z[j]);

    } while (iteration_state == dealii::SolverControl::iterate);

  this->memory.free(aux);

  dealii::deallog.pop();
				   // in case of failure: throw
				   // exception
  if (this->control().last_check() != dealii::SolverControl::success)
    throw dealii::SolverControl::NoConvergence (this->control().last_step(),
					this->control().last_value());
}

#endif // DOXYGEN


#endif
