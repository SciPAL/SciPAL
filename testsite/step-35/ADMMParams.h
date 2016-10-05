//@sect3{File: ADMMParams.h}
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

Copyright  Lutz Künneke, Jan Lebert 2014
*/

#ifndef ADMMPARAMS_H
#define ADMMPARAMS_H

//std
#include <iostream>

//Qt
#include <QDir>
#include <QTime>

// deal.II
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/solver_control.h>

namespace step35 {
//@sect4{enum: regularisationType}
//@param haar regularization by haar wavelet sparsity
//@param sparse regularization by direct space sparsity
//@param quadratic Minimize 2-Norm
// enum regularisationType {haar, sparse, quadratic, TV};


// @sect4{Class: ADMMParams}
//

// @sect4{Class: ADMMParams}
//
// This class contains parameters read from control file, default: step-35.prm
class ADMMParams : public dealii::Subscriptor {

public:

    // (Initial) timestep for Heun
    double pseudo_timestep;

    //Maximum number of Heun steps, number of iteration steps to compute the reconstructed image for given noise.
    int n_max_heun_steps;

    // In each time step of the heun algorithm we adapt the time step this many times.
    int n_max_heun_dt_adaption_steps;

    // Tolerance of the norm of the difference of two successive Heun steps.
    double heun_Tol;

    // Tolerance of the norm of the difference of two successive time step attempts within a Heun step.
    double dt_adaption_Tol;




    //Report interval of ADMM
    int report_interval;



    //Wether or not to use approximative projections
    // bool do_approx; // FIXME: give list of different Dykstra strategies




    static void declare(dealii::ParameterHandler &prm);

    void get(dealii::ParameterHandler &prm);

    dealii::SolverControl solver_control;

private:
   // ADMMParams (const ADMMParams & /*other*/) {}

    ADMMParams & operator = (const ADMMParams & /*other*/) {
        return *this;
    }
};

//@sect5{Function: declare}
//@brief Declaration our parameters
void ADMMParams::declare(dealii::ParameterHandler &prm) {


    prm.enter_subsection ("Solver control");

    prm.declare_entry ("Report interval",
                       "1",
                       dealii::Patterns::Integer(),
                       "Reporting progress in intervals of ... Iterations");



    dealii::SolverControl::declare_parameters(prm);

    prm.enter_subsection("Heun");
    {
        prm.declare_entry ("N max time steps",
                           "1000",
                           dealii::Patterns::Integer(),"Maximum number of Heun steps, number of iteration steps to compute the reconstructed image for given noise.");

        prm.declare_entry ("N max time step adaptions",
                           "10",
                           dealii::Patterns::Integer(),"");

        prm.declare_entry ("Initial pseudo-timestep",
                           "1e-4",
                           dealii::Patterns::Double(),
                           "Initial value for the timestep in the very first Heun step. "
                           "Afterwards the Heun solver uses the timestep from the previous iteration.");


        prm.declare_entry ("Tolerance",
                           "1e-4",
                           dealii::Patterns::Double(),
                           "Tolerance of the norm of the difference of two successive Heun steps.");

        prm.declare_entry ("Time step adaption tolerance",
                           "1e-4",
                           dealii::Patterns::Double(),
                           "Tolerance of the norm of the difference of two successive time step attempts within a Heun step.");

    }
    prm.leave_subsection();


    prm.leave_subsection ();
}

//@sect5{Function: get}
//@brief Read the parameters from file or get default value
void ADMMParams::get(dealii::ParameterHandler &prm) {


    prm.enter_subsection ("Solver control");

    solver_control.parse_parameters(prm);

    this->report_interval = prm.get_integer("Report interval");


    prm.enter_subsection("Heun");
    {
        this->n_max_heun_steps = prm.get_integer("N max time steps");

        this->n_max_heun_dt_adaption_steps = prm.get_integer ("N max time step adaptions");

        this->pseudo_timestep = prm.get_double("Initial pseudo-timestep");

        this->heun_Tol = prm.get_double("Tolerance");

        this->dt_adaption_Tol = prm.get_double("Time step adaption tolerance");
    }
    prm.leave_subsection();


    prm.leave_subsection ();
}


} // namespace step35 END

#endif // ADMMPARAMS_H
