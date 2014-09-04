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


#ifndef SIM_PARAMETER_H
#define SIM_PARAMETER_H

#include <deal.II/base/parameter_handler.h>

#include <QDir>

namespace step1 {

// @sect3{Class: SimParams}
//
// This structure contains all parameters necessary for controling
// the global test properties, i.e. precision and what BLAS to use.
// The documentation strings given in the declare function provide
// a detailed documentation of the individual attributes of this class
// and are available at run-time as well.
struct SimParams {

    SimParams() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    int
    device,
    matrix_low,
    matrix_high,
    step_size,
    average_runs;

    bool use_double;

    QDir run_dir;

private:
    // As usual, inhibit automatic generation of copy ctor and assignment operator.
    SimParams(const SimParams& /*other*/) {}

    SimParams& operator= (const SimParams& /*other*/)
    {
        return *this;
    }

};
}


#endif
