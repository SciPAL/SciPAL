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


#ifndef SOLVER_QMRCGSTAB_H
#define SOLVER_QMRCGSTAB_H
#include<QDebug>
#include <float.h>

#include <lac/solver_cuda.h>
#include <deal.II/lac/tridiagonal_matrix.h>
// wieder einbauen spaeter
//#include <lac/PreconditionerAINV.h>


template <class VECTOR = dealii::Vector<double> >
                         class SolverQMRCGSTAB : public Solver<VECTOR>
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

//    struct PreconditionAdapters {

//    template<typename PC>
//        class Impl : public PC {



//        void vmult(dst, src)
//        {

//            // prec(A)*p = p_
//            if (add_data->right_prec)
//            {
//                this->vmult(p_, p);
//                this->A_ptr.vmult(v, p_);
//            }
//            else if (left)
//            {
//               this->A_ptr.vmult(v, p_);
//                this->vmult(p_, p);
//            }

//        }

//        };

//        typedefs Impl<dealii::PreconditionSSOR< > > SSOR;


//    };


    /**
                                      * Constructor.
                                      */
    SolverQMRCGSTAB (dealii::SolverControl        &cn,
                     dealii::VectorMemory<VECTOR> &mem,
                     const AdditionalData &data = AdditionalData());

    /**
                                      * Constructor. Use an object of
                                      * type GrowingVectorMemory as
                                      * a default to allocate memory.
                                      */
    SolverQMRCGSTAB (dealii::SolverControl        &cn,
                     const AdditionalData &data=AdditionalData());

    /**
                                      * Virtual destructor.
                                      */
    virtual ~SolverQMRCGSTAB ();

    /**
                                      * Solve the linear system $Ax=b$
                                      * for x.
                                      */
    template <class MATRIX, class PRECONDITIONER>
            void
            solve (const MATRIX         &A,
                   VECTOR               &x,
                   const VECTOR         &b,
                   const PRECONDITIONER &precondition_right
                   );

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
    VECTOR *Vr_0;
    VECTOR *Vp;
    VECTOR *Vp_;
    VECTOR *Vv;
    VECTOR *Vd;
    VECTOR *Vs;

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
    bool is_zero(double x){ return ( (x<10.0*DBL_MIN) && (x>-10.0*DBL_MIN) ? true : false) ;}

};

/*@}*/

/*------------------------- Implementation ----------------------------*/


template <class VECTOR>
        inline
        SolverQMRCGSTAB<VECTOR>::AdditionalData::
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
        SolverQMRCGSTAB<VECTOR>::SolverQMRCGSTAB (dealii::SolverControl        &cn,
                                                  dealii::VectorMemory<VECTOR> &mem,
                                                  const AdditionalData &data)
                                                      :
                                                      Solver<VECTOR>(cn,mem),
                                                      additional_data(data)
{}



template <class VECTOR>
        SolverQMRCGSTAB<VECTOR>::SolverQMRCGSTAB (dealii::SolverControl        &cn,
                                                  const AdditionalData &data)
                                                      :
                                                      Solver<VECTOR>(cn),
                                                      additional_data(data)
{}



template <class VECTOR>
        SolverQMRCGSTAB<VECTOR>::~SolverQMRCGSTAB ()
{}



template <class VECTOR>
        double
        SolverQMRCGSTAB<VECTOR>::criterion()
{
    return std::sqrt(res2);
}



template <class VECTOR>
        void
        SolverQMRCGSTAB<VECTOR>::cleanup()
{


    this->memory.free(Vr);
    this->memory.free(Vr_0);
    this->memory.free(Vp);
    this->memory.free(Vp_);
    this->memory.free(Vv);
    this->memory.free(Vd);
    this->memory.free(Vs);


    this->memory.free(VAp);
    dealii::deallog.pop();
}



template <class VECTOR>
        void
        SolverQMRCGSTAB<VECTOR>::print_vectors(const unsigned int,
                                               const VECTOR&,
                                               const VECTOR&,
                                               const VECTOR&) const
{}



template <class VECTOR>
        template <class MATRIX, class PRECONDITIONER>
        void
        SolverQMRCGSTAB<VECTOR>::solve (const MATRIX         &A,
                                        VECTOR               &x,
                                        const VECTOR         &b,
                                        const PRECONDITIONER &precondition)
{
    dealii::SolverControl::State conv=dealii::SolverControl::iterate;

    dealii::deallog.push("qmrcgstab");



    // Memory allocation
    Vr    = this->memory.alloc();
    Vr_0  = this->memory.alloc();
    Vp    = this->memory.alloc();
    Vp_   = this->memory.alloc();
    Vv    = this->memory.alloc();
    Vd    = this->memory.alloc();
    Vs    = this->memory.alloc();

    VAp   = this->memory.alloc();



    try {

        VECTOR& r = *Vr;
        VECTOR& r_0 = *Vr_0;
        VECTOR& p = *Vp;
        VECTOR& p_ = *Vp_;
        VECTOR& v = *Vv;
        VECTOR& d = *Vd;
        VECTOR& s = *Vs;
        VECTOR& Ap = *VAp;

        // resize the vectors, but do not set
        // the values since they'd be overwritten
        // soon anyway.


        r.reinit(x/*, false*/);
        r_0.reinit(x/*, false*/);
        p.reinit(x/*, false*/);
        p_.reinit(x/*, false*/);
        v.reinit(x/*, false*/);
        d.reinit(x/*, false*/);
        s.reinit(x/*, false*/);

        Ap.reinit(x/*, false*/);



        // Implementation taken from the DEAL
        // library
        int  it=0;
        double rho, alpha, eta, theta, c, gam, beta, tau, res;

        // compute residual. if vector is
        // zero, then short-circuit the
        // full computation
        if (!x.all_zero())
        {
            A.vmult(r,x);
            r.sadd(-1.,1.,b);
        }
        else
            r = b;

        r_0 = r;
        p = r;

        rho = tau = r*r;
        res = sqrt(tau);

        if( is_zero(tau) ) { printf("tau is zero: %f\n", tau); throw;}

        conv = this->control().check(0,res);
        if (conv)
        {
            cleanup();
            return;
        }

        while (conv == dealii::SolverControl::iterate)
        {


            // prec(A)*p = p_
            precondition.vmult(p_, p);



            A.vmult(v, p_);

//            A.vmult(v,p);

            alpha = v*r_0;
            if( is_zero(alpha) ) { printf("alpha is zero: %f\nin step: %d\n", alpha,this->control().last_step()); throw;}

            alpha = rho / alpha;
//            alpha = rho/(v*r_0);

            if( is_zero(alpha) ) { printf("alpha is zero: %f\nin step: %d\n", alpha,this->control().last_step()); throw;}


            // s = r - alpha*v
            //s.equ(1.0, r, -alpha, v);
            s=r;
            s.sadd(1.0, -alpha, v);

            // d = p_ + (theta*eta/alpha) * d
            //d.sadd((theta*eta/alpha), p_);
            d *= (theta*eta/alpha);

            d.add(p_);
//            d.add(p);
            //d.equ(1.0, p_, (theta*eta/alpha), d);

            theta = s*s/tau;

            c = 1.0/(1.0 + theta);
            eta = c * alpha;
            tau *= theta * c;

            // x = x + eta * d
            x.add(eta, d);


            // prec(A)*p_ = s
            precondition.vmult(p_, s);


            // r = A*p_
            A.vmult(r,p_);

//            A.vmult(r,s);

            // gam = <r,s>/<r,r>
            gam = (r*s)/(r*r);

            if( is_zero(gam) ) { printf("gam is zero: %f\n", gam); throw;}


            // r = s - gam*r
            //r.sadd(-gam, s);
            r *= -gam;
            r.add(s);

            // d = p_ + (theta*eta/gam) * d
            //d.sadd((theta*eta/gam), p_);
            d *= (theta*eta/gam);

            d.add(p_);
//            d.add(s);

            theta = (r*r)/tau;

            c = 1.0/(1.0 + theta);
            eta = c * gam;
            tau *= theta * c;

            // x = x + eta * d
            x.add(eta, d);

            it++;
            res = sqrt( (2.0 * it + 1.0 ) *tau);

            conv = this->control().check(it,res);

            if (conv) break;

            if( is_zero(gam*rho) ) { printf("gam*rho is zero: %f\n", gam*rho); throw;}


            beta = alpha/ (gam*rho);
            rho = r*r_0;
            beta *= rho;

            p.add(-gam, v);
            //p.equ( 1.0,r, beta, p );
            p *= beta;
            p.add(r);
        }

    }
    catch (...)
    {
        cleanup();
//        throw;
    }


    // Deallocate Memory
    cleanup();
    // in case of failure: throw
    // exception



//    if (this->control().last_check() != dealii::SolverControl::success)
//        throw dealii::SolverControl::NoConvergence (this->control().last_step(),
//                                                    this->control().last_value());
    // otherwise exit as normal
}


#endif // SOLVER_QMRCGSTAB_H
