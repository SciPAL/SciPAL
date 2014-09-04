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


#ifndef SOLVER_BICGSTAB_H
#define SOLVER_BICGSTAB_H



#include <lac/solver_cuda.h>


template <class VECTOR = dealii::Vector<double> >
class SolverBicgstab : public Solver<VECTOR>
{
  public:
    				     /**
				      * There are two possibilities to
				      * compute the residual: one is an
				      * estimate using the computed value @p
				      * tau. The other is exact computation
				      * using another matrix vector
				      * multiplication. This increases the
				      * costs of the algorithm, so it is
				      * should be set to false whenever the
				      * problem allows it.
				      *
				      * Bicgstab is susceptible to breakdowns, so
				      * we need a parameter telling us, which
				      * numbers are considered zero.
				      */
    struct AdditionalData
    {
					 /**
					  * Constructor.
					  *
					  * The default is to perform an
					  * exact residual computation and
					  * breakdown parameter 1e-10.
					  */
	AdditionalData(const bool   exact_residual = true,
		       const double breakdown      = 1.e-10) :
			exact_residual(exact_residual),
			breakdown(breakdown)
	  {}
					 /**
					  * Flag for exact computation of residual.
					  */
	bool exact_residual;
					 /**
					  * Breakdown threshold.
					  */
	double breakdown;
    };

				     /**
				      * Constructor.
				      */
    SolverBicgstab (dealii::SolverControl        &cn,
            dealii::VectorMemory<VECTOR> &mem,
		    const AdditionalData &data=AdditionalData());

				     /**
				      * Constructor. Use an object of
				      * type GrowingVectorMemory as
				      * a default to allocate memory.
                      */
    SolverBicgstab (dealii::SolverControl        &cn,
		    const AdditionalData &data=AdditionalData());
    
				     /**
				      * Virtual destructor.
				      */
    virtual ~SolverBicgstab ();

				     /**
				      * Solve primal problem only.
				      */
    template<class MATRIX, class PRECONDITIONER>
    void
    solve (const MATRIX &A,
	   VECTOR       &x,
	   const VECTOR &b,
	   const PRECONDITIONER& precondition);

  protected:
				     /**
				      * Computation of the stopping criterion.
				      */
    template <class MATRIX>
    double criterion (const MATRIX& A, const VECTOR& x, const VECTOR& b);

				     /**
				      * Interface for derived class.
				      * This function gets the current
				      * iteration vector, the residual
				      * and the update vector in each
				      * step. It can be used for a
				      * graphical output of the
				      * convergence history.
				      */
    virtual void print_vectors(const unsigned int step,
			       const VECTOR& x,
			       const VECTOR& r,
			       const VECTOR& d) const;

				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vx;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vr;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vrbar;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vp;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vy;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vz;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vs;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vt;
				     /**
				      * Auxiliary vector.
				      */
    VECTOR *Vv;
				     /**
				      * Right hand side vector.
				      */
    const VECTOR *Vb;
  
				     /**
				      * Auxiliary value.
				      */
    double alpha;
				     /**
				      * Auxiliary value.
				      */
    double beta;
				     /**
				      * Auxiliary value.
				      */
    double omega;
				     /**
				      * Auxiliary value.
				      */
    double rho;
				     /**
				      * Auxiliary value.
				      */
    double rhobar;
  
				     /**
				      * Current iteration step.
				      */
    unsigned int step;
  
				     /**
				      * Residual.
				      */
    double res;

				     /**
				      * Additional parameters.
				      */
    AdditionalData additional_data;
  
  private:
				     /**
				      * Everything before the iteration loop.
				      */
    template <class MATRIX>
    dealii::SolverControl::State start(const MATRIX& A);

				     /**
				      * The iteration loop itself.
				      */
    template<class MATRIX, class PRECONDITIONER>
    bool
    iterate(const MATRIX& A, const PRECONDITIONER& precondition);

    void cleanup();
  
};

/*@}*/
/*-------------------------Inline functions -------------------------------*/



template<class VECTOR>
SolverBicgstab<VECTOR>::SolverBicgstab (dealii::SolverControl &cn,
                    dealii::VectorMemory<VECTOR> &mem,
					const AdditionalData &data)
		:
		Solver<VECTOR>(cn,mem),
		additional_data(data)
{}



template<class VECTOR>
SolverBicgstab<VECTOR>::SolverBicgstab (dealii::SolverControl &cn,
					const AdditionalData &data)
		:
		Solver<VECTOR>(cn),
		additional_data(data)
{}



template<class VECTOR>
SolverBicgstab<VECTOR>::~SolverBicgstab ()
{}



template <class VECTOR>
template <class MATRIX>
double
SolverBicgstab<VECTOR>::criterion (const MATRIX& A, const VECTOR& x, const VECTOR& b)
{
  A.vmult(*Vt, x);
  Vt->sadd(-1.,1.,b);
  res = Vt->l2_norm();
  
  return res;
}



template <class VECTOR >
template <class MATRIX>
dealii::SolverControl::State
SolverBicgstab<VECTOR>::start(const MATRIX& A)
{
  A.vmult(*Vr, *Vx);
  Vr->sadd(-1.,1.,*Vb);
  res = Vr->l2_norm();
  
  Vp->reinit(*Vx);
  Vv->reinit(*Vx);
  *Vrbar = *Vr;
  return this->control().check(step, res);
}



template<class VECTOR>
void
SolverBicgstab<VECTOR>::print_vectors(const unsigned int,
				      const VECTOR&,
				      const VECTOR&,
				      const VECTOR&) const
{}

template <class VECTOR>
void
SolverBicgstab<VECTOR>::cleanup()
{
  this->memory.free(Vr);
  this->memory.free(Vrbar);
  this->memory.free(Vp);
    this->memory.free(Vy);
    this->memory.free(Vz);
    this->memory.free(Vs);
    this->memory.free(Vt);
    this->memory.free(Vv);
  dealii::deallog.pop();
}


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
bool
SolverBicgstab<VECTOR>::iterate(const MATRIX& A,
				const PRECONDITIONER& precondition)
{
//TODO:[GK] Implement "use the length of the computed orthogonal residual" in the BiCGStab method.
  dealii::SolverControl::State state = dealii::SolverControl::iterate;
  alpha = omega = rho = 1.;

  VECTOR& r = *Vr;
  VECTOR& rbar = *Vrbar;
  VECTOR& p = *Vp;
  VECTOR& y = *Vy;
  VECTOR& z = *Vz;
  VECTOR& s = *Vs;
  VECTOR& t = *Vt;
  VECTOR& v = *Vv;
  
  do
    {
      ++step;
      
      rhobar = r*rbar;
      beta   = rhobar * alpha / (rho * omega);
      rho    = rhobar;
      p.sadd(beta, 1., r, -beta*omega, v);
      precondition.vmult(y,p);
      A.vmult(v,y);
      rhobar = rbar * v;

      alpha = rho/rhobar;

//TODO:[?] Find better breakdown criterion

      if (std::fabs(alpha) > 1.e10)
	return true;
      
      s.equ(1., r, -alpha, v);

				       // check for early success, see
				       // the lac/bicgstab_early
				       // testcase as to why this is
				       // necessary
      if (this->control().check(step, s.l2_norm()/Vb->l2_norm())
      == dealii::SolverControl::success)
	{
	  Vx->add(alpha, y);
	  print_vectors(step, *Vx, r, y);
	  return false;
	}

      precondition.vmult(z,s);
      A.vmult(t,z);
      rhobar = t*s;
      omega = rhobar/(t*t);
      Vx->add(alpha, y, omega, z);
      r.equ(1., s, -omega, t);

      if (additional_data.exact_residual)
	res = criterion(A, *Vx, *Vb);
      else
	res = r.l2_norm();
      
      state = this->control().check(step, res);
      print_vectors(step, *Vx, r, y);
    }
  while (state == dealii::SolverControl::iterate);
  return false;
}


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverBicgstab<VECTOR>::solve(const MATRIX &A,
			      VECTOR       &x,
			      const VECTOR &b,
			      const PRECONDITIONER& precondition)
{
  dealii::deallog.push("Bicgstab");
  Vr    = this->memory.alloc(); Vr->reinit(x);
  Vrbar = this->memory.alloc(); Vrbar->reinit(x);
  Vp    = this->memory.alloc();
  Vy    = this->memory.alloc(); Vy->reinit(x);
  Vz    = this->memory.alloc(); Vz->reinit(x);
  Vs    = this->memory.alloc(); Vs->reinit(x);
  Vt    = this->memory.alloc(); Vt->reinit(x);
  Vv    = this->memory.alloc();

  Vx = &x;
  Vb = &b;

  step = 0;

  bool state;
  
  do 
    {
      if (step != 0)
    dealii::deallog << "Restart step " << step << std::endl;
      if (start(A) == dealii::SolverControl::success)
	break;
      state = iterate(A, precondition);
    }
  while (state);

  this->memory.free(Vr);
  this->memory.free(Vrbar);
  this->memory.free(Vp);
  this->memory.free(Vy);
  this->memory.free(Vz);
  this->memory.free(Vs);
  this->memory.free(Vt);
  this->memory.free(Vv);
  
  dealii::deallog.pop();
  
				   // in case of failure: throw
				   // exception
  if (this->control().last_check() != dealii::SolverControl::success)
    throw dealii::SolverControl::NoConvergence (this->control().last_step(),
					this->control().last_value());
				   // otherwise exit as normal
}

#endif //! SOLVER_BICGSTAB_H
