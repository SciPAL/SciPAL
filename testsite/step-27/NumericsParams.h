#ifndef NUMERICSPARAMS_H
#define NUMERICSPARAMS_H

#include <deal.II/base/parameter_handler.h>

#include <deal.II/base/quadrature.h>

#include <step-27/Architecture.h>

// @sect3{struct NumericsParams}
//
namespace step27 {

template<int dim>
struct NumericsParams {

    NumericsParams();


     ~NumericsParams();

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    // Our numerical parameters can be categorized as follows.
    // First, the parameters for the geometry of the computational domain.
    // Since the shape of the domain is fixed to a cylinder with a void,
    // these are only the linear dimensions.
    double cavity_radius, cyl_radius, cyl_half_height;

    // Next, the parameters for the mesh refinement and the finite and boundary element methods.
    uint n_init_refinements, n_mesh_refinements;

    uint fe_mapping_degree, fe_degree;

    int bem_quad_order_1;
    int bem_quad_order_2;

    // Last but not least, the parameters for the multigrid method.
    int n_smoothing_steps;

    // The last attribute only minimizes the number of arguments of the constructor of the FEMBEMPoissonProblem class.
    std::string timing_file;

    // For later steps extending this one to the full DRS problem we still need a few spare parameters.
    uint n_schur_iter, n_threads, n_time_steps;

    double pseudo_time_step;

    dealii::Quadrature<dim> * fem_q_rule;

    // In contrast to finite element methods boundary element methods may result in double integrals
    // when based on the weak formulation. To remove the necessity to later add the quadrature rule
    // for the outer integral by subclassing this parameter structure we declare it already here.
    dealii::Quadrature<dim-1> * bem_q_rule_1;

    dealii::Quadrature<dim-1> * bem_q_rule_2;

    // Parameters set by the application and not by the parameter handler.
    std::string prm_drs_cell_filepath;

    step27::Architecture arch;

private:
    // Because of the pointers to the quadrature rules we define
    // the copy constructor as private, so that no-one can accidentally use it.
    NumericsParams(const NumericsParams& /*other*/);
};

// @sect3{struct: LaplaceKernel}
//
//
template <int dim>
struct LaplaceKernel
{

    static double single_layer(const dealii::Point<dim> &R)
    {
        switch(dim)
        {
        case 2:
            return (-std::log(R.norm()) / (2*dealii::numbers::PI) );

        case 3:
//            return (1./( R.norm()*4*numbers::PI ) );
            return (1./( R.norm() ) );

        default:
            Assert(false, dealii::ExcInternalError());
            return 0.;
        }
    }



    static const dealii::Point<dim> double_layer(const dealii::Point<dim> &R)
    {
        switch(dim)
        {
        case 2:
            return R / ( +2*dealii::numbers::PI * R.square());
        case 3:
//            return R / ( +4*numbers::PI * R.square() * R.norm() );
            return R / ( R.square() * R.norm() );

        default:
            Assert(false, dealii::ExcInternalError());
            return dealii::Point<dim>();
        }
    }

    static double hyper(const dealii::Point<dim> &R, const dealii::Point<dim> &na, const dealii::Point<dim> &ni)
    {
        switch(dim)
        {
        //TODO: 2d
        case 3:
            return -3*(na * R)*(ni * R) / (R.square() * R.square() * R.norm()) + (na*ni) / (R.square() * R.norm());

        default:
            Assert(false, dealii::ExcInternalError());
            return 0;
        }
    }
};



} // namespace step27 END

#endif // NUMERICSPARAMS_H
