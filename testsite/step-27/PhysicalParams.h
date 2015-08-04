#ifndef PHYSICALPARAMS_H
#define PHYSICALPARAMS_H

#include <iostream>
#include <string>
#include <deal.II/base/parameter_handler.h>

namespace step27
{

// @sect4{struct PhysicalParams}
//
// This structure collects all the parameters which tune the physical behavior of the system.
// In the FEM-BEM Poisson problem most of them are not needed. However, to simplify the extension of the model
// we collect them all here.
// These are the knobs for your Phys Rev Lett.
struct PhysicalParams {

    double potential_difference;

    // These two numbers are the reason for the complicated nature of the problem.
    // We have two different dielectric constants in the cavity and in the aqueous exterior
    // which gives rise to a dielectric interface, the protein's surface.
    double EPSILON_S;
    double EPSILON_P;

    // The properties of the electrodes enter via the redox rates.
    double k_O;
    double k_R;

    // Because of the conservative nature of the DRS problem we have to prescribe the average densities of the ion species.
    double cations_av_density;
    double anions_av_density;
    double neutrals_av_density;

    PhysicalParams()
        :
          potential_difference(+1.),
          EPSILON_S (80.),
          EPSILON_P (4.),
          k_O (1.0e+0),
          k_R (1.0e+0),
          cations_av_density(0.1),
          anions_av_density(0.1),
          neutrals_av_density(0.1)
    {}

    static std::string subsection_name() { return "Dielectric relaxation spectroscopy setup"; }

    static void declare( dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler &prm);

    PhysicalParams( const PhysicalParams &/*other*/) {}

    PhysicalParams & operator=( PhysicalParams &/*other*/) { return *this; }
};

}

#endif // PHYSICALPARAMS_H
