// Dummy pragma to trick doxygen
#pragma
// <h2>step-27.cpp</h2>
// @sect3{File: step-27.cpp}

// TODO: secondary intro


// This is standard C++.
#include <iostream>
#include <vector>


// Include the actual dielectric relaxation spectroscopy simulation,
// courtesy of Stephan Kramer
#include <drs_simulation.hh>

// For the timing
#include <timing.h>

// Include the SimulationHandler class
#include <simulation_handler.h>

//deal.II includes
#include <deal.II/lac/vector.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>


//for Qt
#include <QtCore>
#include <QDebug>
#include <QVector>
#include <QTextStream>


#include <step-27/assemble_system.hh>


// @sect4{Function: main}
//
// The main function is pretty boring.
int main(int argc, char *argv[])
{
    using namespace step27;

    // In contrast to the other examples, this program is entirely CPU-based.
    // Hence, we do not have to care about which GPUs are available as in the other programs at this place.

    // Instead, we set up the DRS simulation right away ...
    SimulationManager< step27::FEMBEMPoissonProblem<3>, step27::SimParams > simulation(argc, argv);

    // and run it
    return simulation.run();
}
