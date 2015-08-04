#include <step-27/NumericsParams.h>

#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>

template<int dim>
step27::NumericsParams<dim>::NumericsParams()
    :
fem_q_rule(0),
  bem_q_rule_1(0),
  bem_q_rule_2(0),
  arch(step27::cpu)
{}

template<int dim>
step27::NumericsParams<dim>::NumericsParams(const step27::NumericsParams<dim>& /*other*/) {}

template<int dim>
step27::NumericsParams<dim>::~NumericsParams()
{
    std::cout << __FUNCTION__ << std::endl;
    if (fem_q_rule)
        delete fem_q_rule;

    if (bem_q_rule_1)
        delete bem_q_rule_1;

    if (bem_q_rule_2)
        delete bem_q_rule_2;
}


template<int dim>
void step27::NumericsParams<dim>::declare(dealii::ParameterHandler &prm)
{

    prm.enter_subsection("Numerics");


    prm.declare_entry("FE mapping degree", "2",
                      dealii::Patterns::Integer(),
                      "Mapping degree for the finite element method. For FE Degree >= 2, mapping degree and fe degree should be equal."
                      );

    prm.declare_entry("FE degree", "1",
                      dealii::Patterns::Integer(),
                      "Degree of the finite elements.");


    std::string q_rule_doc =
            std::string("Select the type of quadrature rule to use for the finite element part."
            "Possible values are:\n") +
            dealii::QuadratureSelector<dim>::get_quadrature_names();

    prm.declare_entry("Quadrature rule for FE part", "gauss",
                      dealii::Patterns::Selection( dealii::QuadratureSelector<dim>::get_quadrature_names() ),
                      q_rule_doc);

    prm.declare_entry("FE quad order", "-1",
                      dealii::Patterns::Integer(),
                      "Quadrature order for the finite element method. If set to -1, the order is set automatically "
                      "to 2*fe_degree+1 for rules, where the order can be adjusted.");


    prm.declare_entry("BEM quad order inner integral", "4",
                      dealii::Patterns::Integer(),
                      "Quadrature order for the inner integral of"
                      " the boundary element method.");


    prm.declare_entry("Quadrature rule for outer integral of BE part", "gauss",
                      dealii::Patterns::Selection( dealii::QuadratureSelector<dim>::get_quadrature_names() ),
                      "For collocation methods this parameter is ignored, because the outer integral does not exist. "
                      "As q rules you can choose from the same set as for the FE part.");

    prm.declare_entry("BEM quad order outer integral", "-1",
                      dealii::Patterns::Integer(),
                      "Quadrature order for the outer integral of"
                      " the boundary element method (integral over test function in weak formulation!). ");


    prm.declare_entry("N init refinements", "0",
                      dealii::Patterns::Integer(),
                      "After parsing, the sceleton mesh is refined this many times for generating the coarse mesh "
                      "which is the starting point for the adaptove solutino process.");

    prm.declare_entry("N mesh refinements", "3",
                      dealii::Patterns::Integer(),
                      "During mesh refinement the coarse mesh may be refined "
                      "this many times in order to compute an improved solution.");


    prm.declare_entry("N Schur Iterations", "15",
                      dealii::Patterns::Integer(),
                      "Adds to global settings");

    prm.declare_entry("N threads", "4",
                      dealii::Patterns::Integer(),
                      "Number of threads for the cell loop in assemble_system(). Each thread starts its own CUDA kernels for the BEM part.");



    prm.declare_entry("N pseudo time steps", "5",
                      dealii::Patterns::Integer(),
                      "Number steps to iterate the nonlinearity before moving to the next mesh.");


    prm.declare_entry("Pseudo timestep", "1e+1",
                      dealii::Patterns::Double(),
                      "Stabilization parameter for the pure von Neumann subproblems.");


    prm.enter_subsection("Multigrid parameters");

    prm.declare_entry("N smoothing steps", "2",
                      dealii::Patterns::Integer(),
                      "");

    prm.declare_entry("Smoothing method",
                      "relaxation",
                      dealii::Patterns::Selection("relaxation|chebyshev"),
                      "The classical choice of smoother is a relaxation method like Jacobi or Gauss-Seidel. "
                      "However, Krylov-like methods like the iterative inverse based on Chebyshev polynomials "
                      "is an option as well. The advantage of the Chebyshev smoother is "
                      "that it only requires matrix-vector products and thus can be parallelized easily.");

    prm.leave_subsection();


    prm.leave_subsection();

}

template<int dim>
void step27::NumericsParams<dim>::get(dealii::ParameterHandler &prm)
{

    prm.enter_subsection("Numerics");


    fe_mapping_degree = prm.get_integer("FE mapping degree");

    fe_degree = prm.get_integer("FE degree");

    std::string fe_q_rule_name = prm.get("Quadrature rule for FE part");
    int fe_q_order = prm.get_integer("FE quad order");



    bem_quad_order_1 = prm.get_integer("BEM quad order inner integral");

    std::string bem_q_rule_name_2 = prm.get("Quadrature rule for outer integral of BE part");
    bem_quad_order_2 = prm.get_integer("BEM quad order outer integral");


    n_init_refinements = prm.get_integer("N init refinements");

    n_mesh_refinements = prm.get_integer("N mesh refinements");


    n_schur_iter = prm.get_integer("N Schur Iterations");

    n_threads= prm.get_integer("N threads");

    n_time_steps = prm.get_integer("N pseudo time steps");

    pseudo_time_step = prm.get_double ("Pseudo timestep");



    prm.enter_subsection("Multigrid parameters");

    n_smoothing_steps = prm.get_integer("N smoothing steps");

    prm.leave_subsection();

    prm.leave_subsection();

    if (fe_q_rule_name == "gauss")
        this->fem_q_rule = new dealii::QuadratureSelector<dim>(fe_q_rule_name, (fe_q_order == -1 ? fe_degree+1 : fe_q_order) );
    else
        this->fem_q_rule = new dealii::QuadratureSelector<dim>(fe_q_rule_name);

    // The inner integrals of the BEM part are always evaluated by Gauss rules.

     this->bem_q_rule_1 = new dealii::QGauss<dim-1>( bem_quad_order_1 );

    // TO DO: if bem_quad_order_2 = 1: use collocation
    if (bem_q_rule_name_2 == "gauss")
        this->bem_q_rule_2 = new dealii::QuadratureSelector<dim-1>(bem_q_rule_name_2, (bem_quad_order_2 == -1 ? fe_degree+1 : bem_quad_order_2) );
    else
        this->bem_q_rule_2 = new dealii::QuadratureSelector<dim-1>(bem_q_rule_name_2);

}



template class step27::NumericsParams<2>;

template class step27::NumericsParams<3>;

