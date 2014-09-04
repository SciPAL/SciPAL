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


//---------------------------------------------------------------------------
//    $Id: solver_selector.h 23876 2011-06-28 18:21:51Z kanschat $
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__solver_selector_h
#define __deal2__solver_selector_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
//#include <deal.II/lac/solver.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
//#include <deal.II/lac/solver_cg.h>
//#include <deal.II/lac/solver_bicgstab.h>
//#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector_memory.h>
//#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition.h>

#include <lac/solver_cuda.h>
#include <lac/solver_cg_cuda.h>
#include <lac/solver_bicgstab_cuda.h>
#include <lac/solver_gmres_cuda.h>
#include <lac/solver_qmrcgstab_cuda_mod.h>

#include <string>
#include <algorithm>

/*!@addtogroup Solvers */
/*@{*/

/**
 * Selects a solver by changing a parameter.
 * 
 * By calling the @p solve function of this @p SolverSelector, it selects
 * the @p solve function of that @p Solver that was specified in the constructor
 * of this class.
 *
 * <h3>Usage</h3>
 * The simplest use of this class is the following:
 * @verbatim
 *                                  // generate a @p SolverControl and
 *                                  // a @p VectorMemory
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // Line 3:
 *                                  //
 *                                  // generate a @p SolverSelector that
 *                                  // calls the @p SolverCG
 * SolverSelector<Vector<double> > 
 *   solver_selector("cg", control, memory);
 *                                  // generate e.g. a @p PreconditionRelaxation
 * PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
 *   preconditioning(A, &SparseMatrix<double>
 *                   ::template precondition_SSOR<double>,0.8);
 *                                  // call the @p solve function with this
 *                                  // preconditioning as last argument
 * solver_selector.solve(A,x,b,preconditioning);
 * @endverbatim
 * But the full usefulness of the @p SolverSelector class is not
 * clear until the presentation of the following example that assumes
 * the user using the @p ParameterHandler class and having declared a
 * "solver" entry, e.g. with
 * @verbatim
 * Parameter_Handler prm;
 * prm.declare_entry ("solver", "none",
 *                    Patterns::Selection(SolverSelector<>::get_solver_names()));
 * ...
 * @endverbatim
 * Assuming that in the users parameter file there exists the line
 * @verbatim
 * set solver = cg
 * @endverbatim
 * then `Line 3' of the above example reads
 * @verbatim
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector(prm.get("solver"), control, memory);
 * @endverbatim
 *
 *
 * If at some time there exists a new solver "xyz" then the user does not need
 * to change his program. Only in the implementation of the @p SolverSelector
 * the calling of this solver has to be added and each user with program lines
 * quoted above only needs to 'set solver = xyz' in his parameter file to get
 * access to that new solver.  :-)
 *
 * (By the way, thanks to Wolfgang for implementing the @p ParameterHandler.)
 * 
 * @author Ralf Hartmann, 1999
 */
template <class VECTOR = dealii::Vector<double> >
class SolverSelector : public dealii::Subscriptor
{
  public:

				     /**
				      * Constructor, filling in
				      * default values
				      */
    SolverSelector ();
  
				     /**
				      * @deprecated Use deafult
				      * constructor and select.
				      *
				      * Constructor. Use the arguments
				      * to initialize actual solver
				      * objects. The VectorMemory
				      * argument is ignored.
				      */
    SolverSelector (const std::string    &solvername,
            dealii::SolverControl        &control
            );
  
				     /**
				      * Destructor
				      */
    ~SolverSelector();

				     /**
				      * Solver procedure. Calls the @p solve
				      * function
				      * of the @p solver whose @p SolverName
				      * was specified in the constructor.
				      * 
				      */
    template<class Matrix, class Preconditioner>
    void solve(const Matrix &A,
	       VECTOR &x,
	       const VECTOR &b,
	       const Preconditioner &precond) const;

				     /**
				      * Select a new solver. Note that
				      * all solver names used in this
				      * class are all lower case.
				      */
    void select(const std::string& name);

				     /**
				      * Set a new SolverControl. This needs to
				      * be set before solving.
				      */
    void set_control(dealii::SolverControl & ctrl);
    
				     /**
				      * Set the additional data. For more
				      * info see the @p Solver class.
				      */
    void set_data(const typename SolverQMRCGSTAB<VECTOR>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the @p Solver class.
				      */
    void set_data(const typename SolverCG<VECTOR>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the @p Solver class.
				      */
    void set_data(const typename SolverBicgstab<VECTOR>
		  ::AdditionalData &data); 

				     /**
				      * Set the additional data. For more
				      * info see the @p Solver class.
				      */
    void set_data(const typename SolverGMRES<VECTOR>
		  ::AdditionalData &data);

				     /**
				      * Set the additional data. For more
				      * info see the @p Solver class.
				      */
    void set_data(const typename SolverFGMRES<VECTOR>
		  ::AdditionalData &data);

				     /**
				      * Get the names of all implemented
				      * solvers.
				      */
    static std::string get_solver_names ();

				     /**
				      * Exception.
				      */
    DeclException1 (ExcSolverDoesNotExist,
		    std::string, << "Solver " << arg1 << " does not exist. Use one of "
		    << std::endl << get_solver_names());

    

  protected:
				     /**
				      * Stores the @p SolverControl that is
				      * needed in the constructor of each @p
				      * Solver class. This can be changed with
				      * @p set_control().
				      */
    dealii::SmartPointer< dealii::SolverControl, SolverSelector< VECTOR > > 	control;

				     /**
				      * Stores the name of the solver.
				      */
    mutable std::string solver_name;
    
  private:


				     /**
				      * Stores the additional data.
				      */
    typename SolverCG<VECTOR>::AdditionalData cg_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverBicgstab<VECTOR>::AdditionalData bicgstab_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverGMRES<VECTOR>::AdditionalData gmres_data;

				     /**
				      * Stores the additional data.
				      */
    typename SolverFGMRES<VECTOR>::AdditionalData fgmres_data;

    /**
     * Stores the additional data.
     */
typename SolverFGMRES<VECTOR>::AdditionalData qmrcgstab_data;
};

/*@}*/
/* --------------------- Inline and template functions ------------------- */


template <class VECTOR>
SolverSelector<VECTOR>::SolverSelector()
{}


template <class VECTOR>
SolverSelector<VECTOR>::SolverSelector(const std::string    &solver_name,
                       dealii::SolverControl        &control) :
		control(&control),
		solver_name(solver_name)
{}


template <class VECTOR>
SolverSelector<VECTOR>::~SolverSelector()
{}


template <class VECTOR>
void
SolverSelector<VECTOR>::select(const std::string& name)
{
  solver_name = name;
}


template <class VECTOR>
template<class Matrix, class Preconditioner>
void
SolverSelector<VECTOR>::solve(const Matrix &A,
			      VECTOR &x,
			      const VECTOR &b,
			      const Preconditioner &precond) const
{

    std::transform(solver_name.begin(), solver_name.end(), solver_name.begin(), ::tolower);
  
  if (solver_name=="qmrcgstab")
    {
      ::SolverQMRCGSTAB<VECTOR> solver(*control);
      solver.solve(A,x,b,precond);
    }       
  else if (solver_name=="cg")
    {
      ::SolverCG<VECTOR> solver(*control, cg_data);
      solver.solve(A,x,b,precond);
    }
  else if (solver_name=="bicgstab")
    {
      ::SolverBicgstab<VECTOR> solver(*control, bicgstab_data);
      solver.solve(A,x,b,precond);
    }
  else if (solver_name=="gmres")
    {
      ::SolverGMRES<VECTOR> solver(*control, gmres_data);
      solver.solve(A,x,b,precond);
    }
  else if (solver_name=="fgmres")
    {
      ::SolverFGMRES<VECTOR> solver(*control,fgmres_data);
      solver.solve(A,x,b,precond);
    }
  else
    Assert(false,ExcSolverDoesNotExist(solver_name));
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_control(
  dealii::SolverControl & ctrl)
{
  control=&ctrl;
}


template <class VECTOR>
std::string SolverSelector<VECTOR>::get_solver_names()
{
  return "qmrcgstab|cg|bicgstab|gmres|fgmres";
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_data(
  const typename SolverGMRES<VECTOR>::AdditionalData &data)
{ 
  gmres_data=data; 
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_data(
  const typename SolverFGMRES<VECTOR>::AdditionalData &data)
{ 
  fgmres_data=data; 
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_data(
  const typename SolverQMRCGSTAB<VECTOR>::AdditionalData &data)
{ 
  qmrcgstab_data=data;
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_data(
  const typename SolverCG<VECTOR>::AdditionalData &data) 
{ 
  cg_data=data; 
}


template <class VECTOR>
void SolverSelector<VECTOR>::set_data(
  const typename SolverBicgstab<VECTOR>::AdditionalData &data) 
{ 
  bicgstab_data=data; 
}


#endif
