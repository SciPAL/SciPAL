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


#ifndef SOLVER_CG_H
#define SOLVER_CG_H


#include <lac/solver_cuda.h>
#include <deal.II/lac/tridiagonal_matrix.h>

template <class VECTOR = dealii::Vector<double> >
                         class SolverCG : public ::Solver<VECTOR>
{
  public:
                                     /**
                                      * Standardized data struct to pipe
                                      * additional data to the solver.
                                      */
    struct AdditionalData
    {
                                         /**
                                          * Write coefficients alpha and beta
                                          * to the log file for later use in
                                          * eigenvalue estimates.
                                          */
        bool log_coefficients;

                                         /**
                                          * Compute the condition
                                          * number of the projected
                                          * matrix.
                                          *
                                          * @note Requires LAPACK support.
                                          */
        bool compute_condition_number;

                                         /**
                                          * Compute the condition
                                          * number of the projected
                                          * matrix in each step.
                                          *
                                          * @note Requires LAPACK support.
                                          */
        bool compute_all_condition_numbers;

                                         /**
                                          * Compute all eigenvalues of
                                          * the projected matrix.
                                          *
                                          * @note Requires LAPACK support.
                                          */
        bool compute_eigenvalues;

                                         /**
                                          * Constructor. Initialize data
                                          * fields.  Confer the description of
                                          * those.
                                          */

        AdditionalData (const bool log_coefficients = false,
                        const bool compute_condition_number = false,
                        const bool compute_all_condition_numbers = false,
                        const bool compute_eigenvalues = false);
    };

                                     /**
                                      * Constructor.
                                      */
    SolverCG (dealii::SolverControl        &cn,
              dealii::VectorMemory<VECTOR> &mem,
              const AdditionalData &data = AdditionalData());

                                     /**
                                      * Constructor. Use an object of
                                      * type GrowingVectorMemory as
                                      * a default to allocate memory.
                                      */
    SolverCG (dealii::SolverControl        &cn,
              const AdditionalData &data=AdditionalData());

                                     /**
                                      * Virtual destructor.
                                      */
    virtual ~SolverCG ();

                                     /**
                                      * Solve the linear system $Ax=b$
                                      * for x.
                                      */
    template <class MATRIX, class PRECONDITIONER>
    void
    solve (const MATRIX         &A,
           VECTOR               &x,
           const VECTOR         &b,
           const PRECONDITIONER &precondition);

  protected:
                                     /**
                                      * Implementation of the computation of
                                      * the norm of the residual. This can be
                                      * replaced by a more problem oriented
                                      * functional in a derived class.
                                      */
    virtual double criterion();

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
                                      * Temporary vectors, allocated through
                                      * the @p VectorMemory object at the start
                                      * of the actual solution process and
                                      * deallocated at the end.
                                      */
    VECTOR *Vr;
    VECTOR *Vp;
    VECTOR *Vz;
    VECTOR *VAp;

                                     /**
                                      * Within the iteration loop, the
                                      * square of the residual vector is
                                      * stored in this variable. The
                                      * function @p criterion uses this
                                      * variable to compute the convergence
                                      * value, which in this class is the
                                      * norm of the residual vector and thus
                                      * the square root of the @p res2 value.
                                      */
    double res2;

                                     /**
                                      * Additional parameters.
                                      */
    AdditionalData additional_data;

  private:
    void cleanup();
};

/*@}*/

/*------------------------- Implementation ----------------------------*/


template <class VECTOR>
inline
SolverCG<VECTOR>::AdditionalData::
AdditionalData (const bool log_coefficients,
                const bool compute_condition_number,
                const bool compute_all_condition_numbers,
                const bool compute_eigenvalues)
                :
                log_coefficients (log_coefficients),
                compute_condition_number(compute_condition_number),
                compute_all_condition_numbers(compute_all_condition_numbers),
                compute_eigenvalues(compute_eigenvalues)
{}


template <class VECTOR>
SolverCG<VECTOR>::SolverCG (dealii::SolverControl        &cn,
                            dealii::VectorMemory<VECTOR> &mem,
                            const AdditionalData &data)
                :
                Solver<VECTOR>(cn,mem),
                additional_data(data)
{}



template <class VECTOR>
SolverCG<VECTOR>::SolverCG (dealii::SolverControl        &cn,
                            const AdditionalData &data)
                :
                Solver<VECTOR>(cn),
                additional_data(data)
{}



template <class VECTOR>
SolverCG<VECTOR>::~SolverCG ()
{}



template <class VECTOR>
double
SolverCG<VECTOR>::criterion()
{
  return std::sqrt(res2);
}



template <class VECTOR>
void
SolverCG<VECTOR>::cleanup()
{
  this->memory.free(Vr);
  this->memory.free(Vp);
  this->memory.free(Vz);
  this->memory.free(VAp);
  dealii::deallog.pop();
}



template <class VECTOR>
void
SolverCG<VECTOR>::print_vectors(const unsigned int,
                                const VECTOR&,
                                const VECTOR&,
                                const VECTOR&) const
{}



template <class VECTOR>
template <class MATRIX, class PRECONDITIONER>
void
SolverCG<VECTOR>::solve (const MATRIX         &A,
                         VECTOR               &x,
                         const VECTOR         &b,
                         const PRECONDITIONER &precondition)
{
  dealii::SolverControl::State conv=dealii::SolverControl::iterate;

  dealii::deallog.push("cg");

                                   //! Memory allocation
  Vr  = this->memory.alloc();
  Vp  = this->memory.alloc();
  Vz  = this->memory.alloc();
  VAp = this->memory.alloc();
                                   //! Should we build the matrix for
                                   //! eigenvalue computations?
  bool do_eigenvalues = additional_data.compute_condition_number
                        | additional_data.compute_all_condition_numbers
                        | additional_data.compute_eigenvalues;
  double eigen_beta_alpha = 0;

                                   //! vectors used for eigenvalue
                                   //! computations
  std::vector<double> diagonal;
  std::vector<double> offdiagonal;

  try {
                                     //! define some aliases for simpler access
    VECTOR& g  = *Vr;
    VECTOR& h  = *Vp;
    VECTOR& d  = *Vz;
    VECTOR& Ad = *VAp;
                                     //! resize the vectors, but do not set
                                     //! the values since they'd be overwritten
                                     //! soon anyway.
    g.reinit(x/*, true*/);
    h.reinit(x/*, true*/);
    d.reinit(x/*, true*/);
    Ad.reinit(x/*, true*/);
                                     //! Implementation taken from the DEAL
                                     //! library
    int  it=0;
    double res,gh,alpha,beta;

                                     //! compute residual. if vector is
                                     //! zero, then short-circuit the
                                     //! full computation
    if (!x.all_zero())
      {
        A.vmult(g,x);
        g.sadd(-1.,1.,b);
      }
    else
      g = b;

    res = g.l2_norm();

    std::cout << "b is: "<< b.l2_norm() << std::endl;
    //g.print(std::cout);

    conv = this->control().check(0,res);
    if (conv)
      {
        cleanup();
        return;
      }

    g *= -1.;
    precondition.vmult(h,g);


    //d.equ(-1.,h);
    d = h;
    d *= -1;

    gh = g*h;

    while (conv == dealii::SolverControl::iterate)
      {
        it++;
        A.vmult(Ad,d);

        alpha = d*Ad;
        alpha = gh/alpha;

        g.add(alpha,Ad);
        x.add(alpha,d );
        res = g.l2_norm();

        print_vectors(it, x, g, d);

        conv = this->control().check(it,res);
        if (conv)
          break;

        precondition.vmult(h,g);

        beta = gh;
        gh   = g*h;
        beta = gh/beta;

        if (additional_data.log_coefficients)
          dealii::deallog << "alpha-beta:" << alpha << '\t' << beta << std::endl;
                                         //! set up the vectors
                                         //! containing the diagonal
                                         //! and the off diagonal of
                                         //! the projected matrix.
        if (do_eigenvalues)
          {
            diagonal.push_back(1./alpha + eigen_beta_alpha);
            eigen_beta_alpha = beta/alpha;
            offdiagonal.push_back(std::sqrt(beta)/alpha);
          }

        if (additional_data.compute_all_condition_numbers && (diagonal.size()>1))
          {
            dealii::TridiagonalMatrix<double> T(diagonal.size(), true);
            for (unsigned int i=0;i<diagonal.size();++i)
              {
                T(i,i) = diagonal[i];
                if (i< diagonal.size()-1)
                  T(i,i+1) = offdiagonal[i];
              }
            T.compute_eigenvalues();
            dealii::deallog << "Condition number estimate: " <<
              T.eigenvalue(T.n()-1)/T.eigenvalue(0) << std::endl;
          }

        d.sadd(beta,-1.,h);
      }
  }
  catch (...)
    {
      cleanup();
      throw;
    }

                                   //! Write eigenvalues or condition number
  if (do_eigenvalues)
    {
      dealii::TridiagonalMatrix<double> T(diagonal.size(), true);
      for (unsigned int i=0;i<diagonal.size();++i)
        {
          T(i,i) = diagonal[i];
          if (i< diagonal.size()-1)
            T(i,i+1) = offdiagonal[i];
        }
      T.compute_eigenvalues();
      if (additional_data.compute_condition_number
          && ! additional_data.compute_all_condition_numbers
          && (diagonal.size() > 1))
        dealii::deallog << "Condition number estimate: " <<
          T.eigenvalue(T.n()-1)/T.eigenvalue(0) << std::endl;
      if (additional_data.compute_eigenvalues)
        {
          for (unsigned int i=0;i<T.n();++i)
            dealii::deallog << ' ' << T.eigenvalue(i);
          dealii::deallog << std::endl;
        }
    }

                                   //! Deallocate Memory
  cleanup();
                                   //! in case of failure: throw
                                   //! exception
  if (this->control().last_check() != dealii::SolverControl::success)
    throw dealii::SolverControl::NoConvergence (this->control().last_step(),
                                        this->control().last_value());
                                   //! otherwise exit as normal
}


#endif //! SOLVER_CG_H
