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


#ifndef SIM_PARAMETER_HH
#define SIM_PARAMETER_HH

#include <deal.II/base/parameter_handler.h>
#include <step-1/SimParams.h>


// @sect4{Function: declare}
//
void
step1::SimParams::declare(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    prm.declare_entry("Run directory", "./test_me",
                      dealii::Patterns::Anything(),
                      "Specify a directory where results of "
                      "the test are to be stored. This can be either "
                      "an absolute path or path relative to the directory "
                      "where the program has been started. The default is "
                      "subdir called test_me-<date> where <date> will be replaced "
                      "by the date at which the program has been started. "
                      "this simplifies keeping the projects directory clean "
                      "");

    prm.leave_subsection();


    prm.enter_subsection("CUDA parameters");


    prm.declare_entry("Device", "0",
                      dealii::Patterns::Integer(),
                      "which CUDA-enabled GPU should be used");

    prm.declare_entry("Shared Memory", "true",
                      dealii::Patterns::Bool(),
                      "Whether shared (true) or L1 (false) memory should be used.");


    prm.leave_subsection();

    prm.enter_subsection("Testcase parameters");

    prm.declare_entry("Double-Precision","true",
                      dealii::Patterns::Bool(),
                      "Decide between double (true) or float (false) precision");

    prm.declare_entry("Matrix size - lower limit", "256",
                      dealii::Patterns::Integer(),
                      "Start value of the range of matrix sizes tested by the simulation");


    prm.declare_entry("Matrix size - upper limit", "513",
                      dealii::Patterns::Integer(),
                      "End value of the range of matrix sizes tested by the simulation");

    prm.declare_entry("Matrix size - step size", "1024",
                      dealii::Patterns::Integer(),
                      "Increment for the size of the test matrices");

    prm.declare_entry("Average - runs", "10",
                      dealii::Patterns::Integer(),
                      "Number of runs being used for averaging");

    prm.leave_subsection();

}


// @sect4{Function: get}
//
void
step1::SimParams::get(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    run_dir.setPath(prm.get("Run directory").c_str());
    run_dir.makeAbsolute();

    if (!run_dir.exists())
        run_dir.mkpath(".");

    prm.leave_subsection();

    prm.enter_subsection("CUDA parameters");

    device        = prm.get_integer("Device");

    prm.leave_subsection();

    prm.enter_subsection("Testcase parameters");


    use_double   = prm.get_bool("Double-Precision");

    matrix_low   = prm.get_integer("Matrix size - lower limit");

    matrix_high  = prm.get_integer("Matrix size - upper limit");

    step_size    = prm.get_integer("Matrix size - step size");

    average_runs = prm.get_integer("Average - runs");


    prm.leave_subsection();


}


#endif
