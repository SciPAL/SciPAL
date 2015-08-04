#ifndef LAPLACEFORM_HH
#define LAPLACEFORM_HH

#include <step-27/bemforms.h>

namespace step27 {



// Evaluation of the particular expression
// is done in a separate class inheriting the FormData.
template<int dim>
struct LaplaceForm : public ::step27::FormData<dim, ::step27::CellIntTraits<dim> > {

    typedef step27::FormData<dim, ::step27::CellIntTraits<dim> > Base;

    dealii::FullMatrix<double>   cell_matrix;
    dealii::FullMatrix<double>   cell_mass_matrix;
    dealii::Vector<double>       cell_rhs;

    static dealii::UpdateFlags update_flags()
    {
        using namespace dealii;
        return update_values | update_support_points| update_gradients |
                update_quadrature_points | update_JxW_values;
    }

    LaplaceForm(const dealii::Mapping<dim> & mapping,
                const dealii::FiniteElement<dim> & fe,
                const dealii::Quadrature<dim> & q_rule)
        :
          Base(mapping, fe, q_rule, update_flags() ),
          cell_matrix (this->dofs_per_cell, this->dofs_per_cell),
          cell_mass_matrix (this->dofs_per_cell, this->dofs_per_cell),
          cell_rhs (this->dofs_per_cell)
    {}

    // The assemble() function assembles the contribution to the
    // left- and right-hand side on a cell.
    template<typename cellIterator>
    void assemble(const cellIterator & cell);

    void assemble_tdep_rhs_contrib(const double timestep,
                                   const std::vector<double> & sol_values,
                                   dealii::Vector<double> & tdep_cell_rhs);

    // Passing the constraints together with the global
    // matrix and right-hand side allows to distribute the local contributions
    // into different global matrices which differ in their constraints.
    // For instance in the location of Dirichlet boundary values.
    void distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                    dealii::SparseMatrix<double> & matrix,
                                    dealii::Vector<double> & rhs);

    void distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                    dealii::SparseMatrix<double> & matrix);
};

} // namespace step27 END





// TO DO: each FORM gets its own header
// The assemble() function assembles the contribution to the
// left- and right-hand side on a cell.
template <int dim>
template<typename cellIterator>
void step27::LaplaceForm<dim>::assemble(const cellIterator & cell) {

    this->Base::reinit(cell);

    cell_mass_matrix = 0.;
    cell_matrix = 0;
    cell_rhs = 0;

    for (uint q_pt=0; q_pt<this->n_q_points; ++q_pt){
        for (uint i=0; i<this->dofs_per_cell; ++i)
        {
            // assemble Laplacian
            for (uint j=0; j<this->dofs_per_cell; ++j)
            {
                cell_matrix(i,j) += // 1.0 * //this->phys_data.DK_H2O *
                        this->shape_grad(i,q_pt) *
                        this->shape_grad(j,q_pt) *
                        this->JxW[q_pt];

                cell_mass_matrix(i,j) +=
                        this->shape_value(i, q_pt)
                        *
                        this->shape_value(j, q_pt)
                        *
                        this->JxW[q_pt];
            }

            // assemble nonlinearity
        }
    }
}



template <int dim>
void step27::LaplaceForm<dim>::assemble_tdep_rhs_contrib(const double timestep,
                                                                       const std::vector<double> & sol_values,
                                                                       dealii::Vector<double> & tdep_cell_rhs)
{
    tdep_cell_rhs = 0.;
    for (uint q_pt=0; q_pt<this->n_q_points; ++q_pt)
    {
        double u_dot_JxW_q = sol_values[q_pt] * this->JxW[q_pt] / timestep;

        for (uint j=0; j<this->dofs_per_cell; ++j)
            tdep_cell_rhs(j) +=
                    this->shape_value(j, q_pt)
                    *
                    u_dot_JxW_q;
    }
}


// Passing the constraints together with the global
// matrix and right-hand side allows to distribute the local contributions
// into different global matrices which differ in their constraints.
// For instance in the location of Dirichlet boundary values.
template <int dim>
void step27::LaplaceForm<dim>::distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                                                        dealii::SparseMatrix<double> & matrix,
                                                                        dealii::Vector<double> & rhs)
{
    constraints.distribute_local_to_global (cell_matrix,
                                            cell_rhs,
                                            this->local_dof_indices,
                                            matrix,
                                            rhs);
}


template <int dim>
void step27::LaplaceForm<dim>::distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                                                        dealii::SparseMatrix<double> & matrix)
{
    constraints.distribute_local_to_global (cell_matrix,
                                            this->local_dof_indices,
                                            matrix);
}





#endif // LAPLACEFORM_HH
