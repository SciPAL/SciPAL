#ifndef TIMING_H
#define TIMING_H

#include <QDebug>
#include <QFile>

// To ease the reuse of things implemented in this project,
// each step gets its own namespace:
namespace step27 {


// @sect4{Class: Timing}
//
// A timer class, which is used to measure the times spent assembling the
// matrices on the GPU using CUDA and on the CPU. The real timing is done by
// the surrounding code, only the measured times are then assigned to members
// of this class. The template parameter is necessary because of NumericsParams<dim>.
template<int dim>
class Timing {

public:

    Timing( const step27::NumericsParams<dim> &numerics_params,
            const uint n_dofs_per_be,
            const uint n_cells,
            const uint n_bem_points,
            const uint cycle);


    // to get the average iteration time, the times for each iteration is
    // summed and eventuelly diveded by the value of this counter to get the
    // average
    uint n_iterations;

    double total_time;
    double average_time_per_iteration;

    // A function to add the current numerics and timing values to a
    // @p dealii::ConvergenceTable. This table is used to store the timings and
    // parameters over all iterations and is, in this case, owned by the
    // @p Step16::FEMBEMPoissonProblem in drs_simulation.h
    void add_values_to_table(dealii::ConvergenceTable &table);

private:
    QFile timing_file;
    const step27::NumericsParams<dim> & numerics_params;

    const uint n_dofs_per_be;
    const uint n_cells;
    const uint n_bem_points;
    const uint cycle;
};

} // namespace step27

// @sect4{Class: Timing}
// @sect5{Constructor: Timing}
//
// Sets up the Timing object, used to keep track of the GPU or CPU execution
// times
//
// @param table : a deal.II Convergence table. This is used to hold the parameters as well as
// the execution times and easily write them to a file
// @param numerics_params : the used numerics parameters for the current run
// @param dofs_per_be, n_cells, n_bem_quad_points, n_bem_points:
// more parameters describing the current problem size
//
// @param cycle: the current refinement step
template<int dim>
step27::Timing<dim>::Timing(const step27::NumericsParams<dim> &numerics_params,
                       uint dofs_per_be, uint n_cells,
                       uint n_bem_points, uint cycle)
    : n_iterations(0),
      total_time(0),
      average_time_per_iteration(0),
      timing_file(numerics_params.timing_file.c_str()),
      numerics_params(numerics_params),
      n_dofs_per_be(dofs_per_be),
      n_cells(n_cells),
      n_bem_points(n_bem_points),
      cycle(cycle)
{}


// @sect5{Function: Timing::writeTimes()}
//
// Member function of the Timing class, used to write the measured times to a file.
// @param table: a @p ConvergenceTable into which the parameters and the measured times
// are written.
template<int dim>
void step27::Timing<dim>::add_values_to_table(dealii::ConvergenceTable &table)
{
    this->average_time_per_iteration /= this->n_iterations;
    qDebug() << "Average: " << this->average_time_per_iteration << "s"
             << "\tTotal: " << this->total_time << "s";


    table.add_value("#cycle", cycle);

    table.add_value("n_dofs_per_be", n_dofs_per_be);
    table.add_value("n_cells", n_cells);
    table.add_value("n_bem_inner_quad_points", numerics_params.bem_quad_order_1);
    table.add_value("n_bem_outer_quad_points", numerics_params.bem_quad_order_2);
    table.add_value("n_bem_points", n_bem_points);
    table.add_value("FE_Degree", numerics_params.fe_degree);
    table.add_value("FE_Mapping_Degree", numerics_params.fe_mapping_degree);

    table.add_value("avg_iter_time", average_time_per_iteration);
    table.set_precision("avg_iter_time", 8);
    table.add_value("total_Time", total_time);
    table.set_precision("total_Time", 8);
}


#endif // TIMING_H
