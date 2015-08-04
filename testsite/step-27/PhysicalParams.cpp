#include <step-27/PhysicalParams.h>


void step27::PhysicalParams::declare( dealii::ParameterHandler & prm)
{
    prm.enter_subsection(subsection_name());

    prm.declare_entry("Potential difference", "1.",
                      dealii::Patterns::Double(),
                      "The potential difference applied to the electrodes.");

    prm.declare_entry("Dielectric constant of the solvent", "80.",
                      dealii::Patterns::Double(),
                      "For water use 80.");

    prm.declare_entry("Dielectric constant of the protein", "4.",
                      dealii::Patterns::Double(),
                      "Reasonable values are 2 - 4.");


    prm.declare_entry("Reduction rate at the cathode", ".1",
                      dealii::Patterns::Double(),
                      "");

    prm.declare_entry("Oxidation rate at the anode", ".1",
                      dealii::Patterns::Double(),
                      "");

    prm.declare_entry("Average cation density", ".1",
                      dealii::Patterns::Double(),
                      "");

    prm.declare_entry("Average density of neutral particles", ".1",
                      dealii::Patterns::Double(),
                      "");

    prm.declare_entry("Average anion density", ".1",
                      dealii::Patterns::Double(),
                      "");

    prm.leave_subsection();
}



void step27::PhysicalParams::get(dealii::ParameterHandler &prm)
{
    prm.enter_subsection(subsection_name());

    potential_difference = prm.get_double("Potential difference");


    EPSILON_S = prm.get_double("Dielectric constant of the solvent");

    EPSILON_P = prm.get_double("Dielectric constant of the protein");

    std::cout << "EPS_H2O : " << EPSILON_S << ", EPSILON_P : " << EPSILON_P << std::endl;

    k_R = prm.get_double("Reduction rate at the cathode");

    k_O = prm.get_double("Oxidation rate at the anode");


    cations_av_density = prm.get_double("Average cation density");

    neutrals_av_density = prm.get_double("Average density of neutral particles");

    anions_av_density = prm.get_double("Average anion density");

    prm.leave_subsection();
}


