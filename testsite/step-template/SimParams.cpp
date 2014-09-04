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


#include <SimParams.h>

// @sect4{Function: declare}
//
// Declare yet another set of parameters
void
steptemplate::SimParams::declare(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    // The most generic parameter for a simulation is a string which indicates
    // where the simulation should store its results. By default, it contains
    // the current date and uses the time as subdirectory.
    prm.declare_entry("Run directory",
                      QString(
                          QString("./")
                          +
                          QString("test_me-"
                                  +
                                  QDateTime::currentDateTime().toString("ddd-yyyy-MM-dd/hh_mm_ss")
                                  ).remove(".")).toStdString(),
                      dealii::Patterns::Anything(),
                      "Specify a directory where the results of "
                      "the test are to be stored. This can be either "
                      "an absolute path or path relative to the directory "
                      "where the program has been started. The default is "
                      "subdir called test_me-<date> where <date> will be replaced "
                      "by the date at which the program has been started. "
                      "this simplifies keeping the projects directory clean");

    prm.leave_subsection();

}


// @sect4{Function: get}
//
// The get function is just a mirror of the declare function.
void
steptemplate::SimParams::get(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    this->run_dir.setPath(prm.get("Run directory").c_str() );
    this->run_dir.makeAbsolute();

    if(!this->run_dir.exists())
        this->run_dir.mkpath(".");

    prm.leave_subsection();

}
