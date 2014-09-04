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


#ifndef SOLVER_H
#define SOLVER_H


#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector_memory.h>


template <class VECTOR = dealii::Vector<double> >
class Solver : public dealii::Subscriptor
{
  public:
                                     /**
                                      * Constructor. Takes a control
                                      * object which evaluates the
                                      * conditions for convergence,
                                      * and an object to provide
                                      * memory.
                                      *
                                      * Of both objects, a reference is
                                      * stored, so it is the user's
                                      * responsibility to guarantee that the
                                      * lifetime of the two arguments is at
                                      * least as long as that of the solver
                                      * object.
                                      */
    Solver (dealii::SolverControl        &solver_control,
            dealii::VectorMemory<VECTOR> &vector_memory);

                                     /**
                                      * Constructor. Takes a control
                                      * object which evaluates the
                                      * conditions for convergence. In
                                      * contrast to the other
                                      * constructor, this constructor
                                      * denotes an internal object of
                                      * type GrowingVectorMemory to
                                      * allocate memory.
                                      *
                                      * A reference to the control
                                      * object is stored, so it is the
                                      * user's responsibility to
                                      * guarantee that the lifetime of
                                      * the two arguments is at least
                                      * as long as that of the solver
                                      * object.
                                      */
    Solver (dealii::SolverControl        &solver_control);

                                     /**
                                      * Access to object that controls
                                      * convergence.
                                      */
    dealii::SolverControl & control() const;

  protected:
                                     /**
                                      * A static vector memory object
                                      * to be used whenever no such
                                      * object has been given to the
                                      * constructor.
                                      */
    mutable dealii::PrimitiveVectorMemory<VECTOR> static_vector_memory;

                                     /**
                                      * Control structure.
                                      */
    dealii::SolverControl &cntrl;

                                     /**
                                      * Memory for auxilliary vectors.
                                      */
    dealii::VectorMemory<VECTOR> &memory;
};

/*-------------------------------- Inline functions ------------------------*/

template<class VECTOR>
inline
Solver<VECTOR>::Solver (dealii::SolverControl        &solver_control,
                        dealii::VectorMemory<VECTOR> &vector_memory)
                :
                cntrl(solver_control),
                memory(vector_memory)
{}



template<class VECTOR>
inline
Solver<VECTOR>::Solver (dealii::SolverControl        &solver_control)
                :
                cntrl(solver_control),
                memory(static_vector_memory)
{}



template <class VECTOR>
inline
dealii::SolverControl &
Solver<VECTOR>::control() const
{
  return cntrl;
}


#endif //! SOLVER_H
