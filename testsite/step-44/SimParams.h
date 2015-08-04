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


#ifndef SIMPARAMS_H
#define SIMPARAMS_H

#include <QDir>
#include <QTime>


// deal.II
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>

namespace step44 {


    // @sect3{Class: SimParams}
    //
    // This class contains parameters to steer the simulation.
    struct SimParams {

        SimParams() {}


        static void declare(dealii::ParameterHandler & prm);

        void get(dealii::ParameterHandler & prm);

        // One useful thing a default implementation can provide
        // is the location where potential results should be stored.
        // This becomes important when the parametric dependence of a
        // physical problem gets investigated.
        QDir run_dir;

        // Secondly, we need a directory where the parameter settings of a run should be logged.
        // This will always be @p run_dir + "/log".
        QDir prm_log_dir;

    private:
        SimParams (const SimParams & /*other*/) {}

        SimParams & operator = (const SimParams & /*other*/) { return *this; }
    };

}


#endif // SIMPARAMS_H
