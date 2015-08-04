#ifndef BEM_FORMS_S27_H
#define BEM_FORMS_S27_H

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe.h>

#include <lac/FullMatrixAccessor.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_values.h>

#include <iostream>

#include <step-27/NumericsParams.h>


namespace step27
{
// To reduce the numbers of headers to include we use some forward declarations.
template <int dim> struct LaplaceKernel;
template<int dim> struct MGData;
template<int dim> struct GeoMGPC;
template <int dim> class FEMBEMPoissonProblem;

// @sect{FormTraits}
//
// With the proper abstraction the assembly of FEM and BEM matrices in FEM-BEM methods
// is virtually the same.
// This is reflected by the fact that with proper templatization both assembly processes can share
// large amounts of code.
// The basis for this code unification are the following two classes which
// provide unfied naems for the spatial dimensiob of the domain of integration
// and the type of evaluator for test and trial functions. For FE methods we need
// deal.II's FeValues class and for BE part we need FeFaceValues.

// Traits for integration over a cell.
template<int dim>
struct CellIntTraits {

    // We provide a flag for the dimension of the domain of integration.
    // For a surface integral this is different from the dimension of the cell.
    static const int dim_int_domain = dim;

    // Last but not least the type implementing the evaluation of the finite element.
    typedef dealii::FEValues<dim> FeValues;
};

// Traits for performing a surface integral.
// The dim template parameter has to be the dimension of the embedding dimension.
// The dimension of the domain of integration for the durface integral then is dim-1.
template<int dim>
struct FaceIntTraits {

    // For a surface integral this is different from the dimension of the cell.
    static const int dim_int_domain = dim-1;

    // Last but not least the type implementing the evaluation of the finite element.
    typedef dealii::FEFaceValues<dim> FeValues;
};


// @sect3{struct FormData}
//
// To shorten the assembly routines we collect the ubiquitous data
// like the FEValues object, the short cut references to e.g. the JxW values, etc ...
// in a class which serves as base for the per-cell evaluation of the weak forms.
// The important thing is, that it subclasses
// the FEeValues class so that we do not have to write e.g.
// <code>fe_values.shape_value(i,q)</code> in order to get
// the value of the ith shape function at quadrature point q, but rather just
// <code>shape_value(i,q)</code>. Since we need similar things
// for the surface integrals we pass a template parameter @p IntegrationTraits
// in order to avoid code duplication.
template<int dim, typename IntegrationTraits>
struct FormData
        :
        public IntegrationTraits::FeValues {

    typedef typename IntegrationTraits::FeValues FeV;

    const unsigned int n_q_points;

    const unsigned int dofs_per_cell;


    const std::vector<dealii::Point<dim> > & q;

    const std::vector<double> & JxW;


    std::vector<uint> local_dof_indices;

    // Since the FEValues object requires an initialized quadrature rule
    // we have to instantiate the quadrature externally. Note that FEValues actually
    // has an attribute of type dealii::Quadrature and not only a pointer.
    // This allows for a very easy combination of weak forms and quadrature rules
    // in case the approximation properties are to be tested as well.
    FormData(const dealii::Mapping<dim> & mapping,
             const dealii::FiniteElement<dim> & fe,
             const dealii::Quadrature<IntegrationTraits::dim_int_domain> & q_rule,
             const dealii::UpdateFlags update_flags)
        :
          FeV (mapping, fe, q_rule, update_flags),
          n_q_points (q_rule.size()),
          dofs_per_cell (fe.dofs_per_cell),
          q (this->FeV::get_quadrature_points() ),
          JxW (this->FeV::get_JxW_values() ),
          local_dof_indices(dofs_per_cell)
    {}

    // Besides collecting all the convenience variables we still provide a function
    // doing all the reinitializations needed on a cell.
    // We have to distinguish active from general, possibly in-active cells.
    template<typename cellIterator>
    inline void reinit (const cellIterator & cell)
    {
        this->FeV::reinit(cell);

        if (cell->active())
            cell->get_dof_indices (local_dof_indices);
        else
            cell->get_mg_dof_indices (local_dof_indices);
    }
};

// @sect4{struct BemFormBase}
//
// Base class for BemCollocationForm and BemWeakForm
// Provides quadrature points, normals, and global to
// local boundary point mapping.
template<int dim, typename IntegrationtTraits>
struct BemFormBase : public FormData<dim, IntegrationtTraits> {

    typedef FormData<dim, IntegrationtTraits> Base;

    static dealii::UpdateFlags update_flags()
    {
        using namespace dealii;
        return update_values|update_normal_vectors|update_quadrature_points|update_JxW_values;
    }

    // Return the normal vectors at the quadrature points.
    // For a face, these are the outward normal vectors to the cell.
    // For a cell of codimension one, the orientation is
    // given by the numbering of vertices.
    const std::vector<dealii::Point<dim> > & normals;

    std::vector<uint> global_bc_dofs;

    static const bool is_col_major = true;

    // const
    std::vector<unsigned int>	dof_to_boundary_mapping;

    BemFormBase(const dealii::Mapping<dim> & mapping,
                const dealii::FiniteElement<dim> & fe,
                const dealii::Quadrature<dim-1> & q_rule,
                const dealii::UpdateFlags update_flags)
        :
          Base(mapping, fe, q_rule, update_flags | this->update_flags() ),
          normals (this->get_normal_vectors () ),
          global_bc_dofs(this->dofs_per_cell)
    {}


    // Besides collecting all the convenience variables we still provide a function
    // doing all the reinitializations needed on a cell.
    template<typename cellIterator>
    void reinit (const cellIterator & cell,
                 const uint face_id,
                 const std::vector<uint>& fem_local_dof_indices)
    {
        this->Base::FeV::reinit(cell, face_id);

        cell->get_dof_indices (this->local_dof_indices);

        for (uint j = 0; j <this->dofs_per_cell; ++j)
            global_bc_dofs[j] = dof_to_boundary_mapping[fem_local_dof_indices[j]];
    }

    template<typename cellIterator>
    inline void reinit (const cellIterator & cell,
                        const uint face_id)
    {
        this->Base::FeV::reinit(cell, face_id);
    }

};



// @sect4{struct DummyCUDADriver}
//
// Before we can define the BemCollocationForm we need a dummy class which provides
// the interface of the CUDADriver of the next step such that we can use it as template argument
// of the BemCollocationForm class (see below). The bodies of the functions
// of this dummy class are just empty. If you want to work with CUDA you have to reimplement this
// class.
// We let the functions throw an exception in case of an accidental use of this implementation.
template <int dim>
class DummyCUDADriver {

public:
    void  reinit(const uint DL_SL_matrix_n_rows,
                 const uint DL_SL_matrix_n_cols) { dealii::ExcNotImplemented(); }

    ~DummyCUDADriver () {}

    template<typename FullMatrixAccessor>
    double assemble_bem_matrices(FullMatrixAccessor &DL_matrix_k,
                                 FullMatrixAccessor &SL_matrix_k,
                                 const std::vector<double> &x_k,
                                 const std::vector<double> &q_k,
                                 const std::vector<double> &n_k,
                                 const FullMatrixAccessor &W_k)
    {
        dealii::ExcNotImplemented();
        return 0;
    }

    void point_to_double_vector(std::vector<double> &out,
                                const std::vector<dealii::Point<dim> > &in)
    {
        dealii::ExcNotImplemented();
    }

};


// @sect4{struct BemCollocationForm}

// This structure adds a set of collocation points to BemFormBase.
// It also provides the per-cell, collocation-based assembly for SL and DL matrices.
// To switch between a cpu- or gpu-based assembly we introduce a template parameter
// @p arch which allows to eliminate the unneeded code sections at compile-time, such
// that we do not have to execute the same switch statement over and over again when looping
// over the faces of the boundary cells.
// The assembly on a face is driven by the assemble_sl_dl member function which delegates
// the actual work to two functions, each for a specific hardware (i.e. cpu or gpu).
// The function arguments are chosen so that they are capable of collocation and weak forms
// of the boundary integrals.
// In the constructor of the FemBem problem we then choose which implementation to take,
// but at run-time and only once.
template<int dim, typename IntegrationTraits, typename CudaDriver=step27::DummyCUDADriver<dim> >
struct BemCollocationForm : public BemFormBase<dim, IntegrationTraits>, protected CudaDriver {

    typedef BemFormBase<dim, IntegrationTraits> Base;

    static dealii::UpdateFlags update_flags()
    {
        using namespace dealii;
        return update_values|update_normal_vectors|update_quadrature_points|update_JxW_values;
    }

    // Collocation points
    std::map<unsigned char, std::vector<dealii::Point<dim> > > support_points;

    /*const*/ uint n_bem_points;

    // W stores the values (shape_value[j,a] * JxW[a]) as a matrix.
    FullMatrixAccessor<double> W;

    BemCollocationForm(const unsigned char bc_id,
                       const dealii::Mapping<dim> & mapping,
                       const dealii::FiniteElement<dim> & fe,
                       const NumericsParams<dim> & numerics_prm)
        :
          Base(mapping, fe, *numerics_prm.bem_q_rule_1, update_flags() ),
          bc_id(bc_id),
          // We have only one component
          component_select(1, true),
          W(numerics_prm.bem_q_rule_1->size(), this->dofs_per_cell, this->is_col_major),
          n_bem_points (0) // x.size() )
    {
        this->interface[bc_id] = 0;

        this->boundary_indicators.insert(bc_id); // DRS::BCids::cavity());

        // can only be done after distribute_dofs:   this->reinit_collocation(fem_dof_handler, mapping);
    }

    // Calls Base::reinit and additionally calculates the W values
    template<typename cellIterator>
    inline void reinit (const cellIterator & cell,
                        const uint face_id,
                        const std::vector<uint>& fem_local_dof_indices)
    {
        this->Base::reinit(cell, face_id, fem_local_dof_indices);

        for (uint j = 0; j <this->dofs_per_cell; ++j)
            for (uint a=0; a < this->n_q_points; ++a)
                W(a,j) = this->shape_value(j, a) * this->JxW[a];

    }

    template<typename cellIterator>
    void assemble_sl_dl(const cellIterator & cell,
                        const uint face_id,
                        dealii::FullMatrix<double> & DL_matrix,
                        dealii::FullMatrix<double> & SL_matrix,
                        const Architecture arch);

    template<typename cellIterator>
    void assemble_sl_dl_cpu(const cellIterator & cell,
                            const uint face_id,
                            dealii::FullMatrix<double> & DL_matrix,
                            dealii::FullMatrix<double> & SL_matrix);

    template<typename cellIterator>
    void assemble_sl_dl_gpu(const cellIterator & cell,
                            const uint face_id,
                            dealii::FullMatrix<double> & DL_matrix,
                            dealii::FullMatrix<double> & SL_matrix);


    const unsigned char bc_id;

    // The BEM- and collocation-specific data
    std::vector< bool > component_select;

    typename dealii::FunctionMap<dim>::type interface;

    unsigned int n_fem_dofs;
    unsigned int n_interior_dofs;
    unsigned int n_bc_dofs;

    std::set<unsigned char>      boundary_indicators;

    std::vector<bool>            dof_flags;

    // This setup is to some extent specific for collocation methods.
    // Since we change the numbering of dofs, we have to pass the dof handler as non-const.
    void reinit_bem_data(/*const*/ dealii::DoFHandler<dim>& fem_dof_handler,
                         const dealii::Mapping<dim>& mapping)
    {

        this->n_fem_dofs  = fem_dof_handler.n_dofs();
        this->n_bc_dofs   = fem_dof_handler.n_boundary_dofs(interface);
        this->n_interior_dofs = this->n_fem_dofs - this->n_bc_dofs;

        n_bem_points = n_bc_dofs;

        dof_flags.resize( n_fem_dofs );

        // TO DO: move following 2 fct calls into bem_form.reinit_collocation
        dealii::DoFTools::extract_boundary_dofs	(fem_dof_handler,
                                                 component_select,
                                                 dof_flags,
                                                 boundary_indicators);

        dealii::DoFRenumbering::sort_selected_dofs_back (fem_dof_handler,
                                                         dof_flags);


        // cavity interface dofs have been sorted to the back of the dof index array
        this->dof_to_boundary_mapping.resize(this->n_fem_dofs,
                                             dealii::DoFHandler<dim>::invalid_dof_index);
        // std::cout << " dof to stbd mapping : " << std::endl;

        for (unsigned int i = 0; i < n_bc_dofs; i++)
            this->dof_to_boundary_mapping[ n_interior_dofs + i] = i;

        // Get support points on FEM-BEM boundary for later evaluation in the integral operators
        BoundaryDoFTools::map_boundary_dofs_to_support_points(mapping,
                                                              fem_dof_handler,
                                                              boundary_indicators,
                                                              this->dof_to_boundary_mapping,
                                                              support_points);

        const unsigned char bc_id = *boundary_indicators.begin();

        typename std::vector<dealii::Point<dim> >::const_iterator
                p = this->support_points[bc_id].begin(),
                endp = this->support_points[bc_id].end();

#ifndef QT_NO_DEBUG
        std::ofstream cavity_support_points_out("cavity_support_points.out");

        for (; p != endp; ++p)
            cavity_support_points_out << *p << std::endl;
        cavity_support_points_out << std::endl << std::endl;
#endif
    }
};
}

#endif // DRS_SIMULATION_H
