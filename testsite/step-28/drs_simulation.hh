/* $Id: step-16.cc 24291 2011-09-09 03:17:50Z bangerth $ */
/* Author: Guido Kanschat, University of Heidelberg, 2003  */
/*         Baerbel Janssen, University of Heidelberg, 2010 */
/*         Wolfgang Bangerth, Texas A&M University, 2010   */

/*    $Id: step-16.cc 24291 2011-09-09 03:17:50Z bangerth $       */
/*                                                                */
/*    Copyright (C) 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

// As discussed in the introduction, most of
// this program is copied almost verbatim
// from step-6, which itself is only a slight
// modification of step-5. Consequently, a
// significant part of this program is not
// new if you've read all the material up to
// step-6, and we won't comment on that part
// of the functionality that is
// unchanged. Rather, we will focus on those
// aspects of the program that have to do
// with the multigrid functionality which
// forms the new aspect of this tutorial
// program.

// @sect3{Include files}

// Again, the first few include files
// are already known, so we won't
// comment on them.

// Qt
#include <QDebug>

// This is C++:
#include <fstream>
#include <sstream>
#include <iomanip>

// For creating directories
#include <sys/stat.h>
#include <sys/types.h>

// The physical parameters have already been declared in step-27.
#include <step-27/PhysicalParams.h>
#include <step-27/NumericsParams.h>

// For choosing the runtime architecture
#include <step-27/Architecture.h>
#include <step-28/coupling_type.h>

#include <step-27/drs_cell.hh>

#include <step-27/BoundaryDoFTools.hh>

#include <step-27/UnitTestProblemData.hh>

// Taylor-made linear algebra classes
// #include <step-27/drs_lac.hh>

#include <step-27/geomgpc.h>

#include <step-27/bemforms.h>
// GPU_PART needed in step-27/bemforms.hh
#define GPU_PART
#include <step-27/bemforms.hh>
#include <step-28/cuda_driver_step-28.h>
#include <step-28/cuda_driver_step-28.hh>
#undef GPU_PART

#include <step-27/drs_simulation.hh>

#define USE_ADAP
#undef USE_ADAP

#define USE_FEM
// #undef USE_FEM


#define USE_FEM_BEM
//#undef USE_FEM_BEM

// #define USE_OLD_STUFF

// #define FEM_BEM_UNIT_TEST



// The last step is as in all
// previous programs:
namespace step28
{

// @sect3{The <code>DRSProblem</code> class template}

// This main class is basically the same
// class as in step-6. As far as member
// functions is concerned, the only addition
// is the <code>assemble_multigrid</code>
// function that assembles the matrices that
// correspond to the discrete operators on
// intermediate levels:
template <int dim_>
class DRSProblem
{

    typedef step27::FemBemData<dim_> FemBemData;

public:
    DRSProblem (step27::NumericsParams<dim_> &numerics_prms,
                const std::string & prm_path,
                const QDir &log_dir,
                step27::TestExpr testcase,
                step27::SubModel sm);


    ~DRSProblem()
    {
        if (bem_form)
            delete bem_form;
        bem_form = 0;
    }

    void run ();

    void block_solve_fix_point(const uint cycle);

#ifdef USE_NEWTON
    void block_solve_Newton();
#endif

    static void point_to_double_vector(std::vector<double> &out,
                                       const std::vector<dealii::Point<dim_> > &in);

    static const int dim = dim_;

protected:
    void setup_discretization ();
    void initialize_solution();





    // Evaluation of the particular expression
    // is done in a separate class inheriting the FormData.
public:
    // Before defining the evaluators for the cell contributions of the nonlinear terms we
    // define some abbreviations for the basic type names.
    typedef step27::FormData<dim, step27::CellIntTraits<dim> > CellFormData;

    typedef step27::FormData<dim, step27::FaceIntTraits<dim> > FaceFormData;


    struct DriftForm : public CellFormData {

        typedef CellFormData Base;

        dealii::FullMatrix<double>   cell_matrix_cations;
        dealii::FullMatrix<double>   cell_matrix_anions;
        dealii::FullMatrix<double>   cell_matrix_drift;


        dealii::Vector<double>       cell_rhs;

        std::vector<double> cation_density_local_values;
        std::vector<double> anion_density_local_values;
        std::vector<dealii::Tensor<1,dim> > Phi_local_values;


        const dealii::Vector<double> & cation_density;
        const dealii::Vector<double> & anion_density;
        const dealii::Vector<double> & potential;

        static dealii::UpdateFlags update_flags()
        {
            using namespace dealii;
            return update_values | update_support_points| update_gradients |
                    update_quadrature_points | update_JxW_values;
        }

        DriftForm(const dealii::Mapping<dim> & mapping,
                  const dealii::FiniteElement<dim> & fe,
                  const dealii::Quadrature<dim> & q_rule,
                  const dealii::Vector<double> & cations_fe_field,
                  const dealii::Vector<double> & anions_fe_field,
                  const dealii::Vector<double> & pot_fe_field)
            :
              Base(mapping, fe, q_rule, update_flags() ),
              cell_matrix_cations (this->dofs_per_cell, this->dofs_per_cell),
              cell_matrix_anions (this->dofs_per_cell, this->dofs_per_cell),
              cell_matrix_drift (this->dofs_per_cell, this->dofs_per_cell),
              cell_rhs (this->dofs_per_cell),
              cation_density_local_values  (this->n_q_points),
              anion_density_local_values  (this->n_q_points),
              Phi_local_values(this->n_q_points),
              cation_density(cations_fe_field),
              anion_density(anions_fe_field),
              potential(pot_fe_field)
        {}


        // The assemble() function assembles the contribution to the
        // left- and right-hand side on a cell.
        template<typename cellIterator>
        void assemble(const cellIterator & cell);


        void distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                        dealii::SparseMatrix<double> & cat_drift_matrix,
                                        dealii::SparseMatrix<double> & an_drift_matrix,
                                        dealii::SparseMatrix<double> & pot_drift_matrix);
    };




    struct FluxForm : public FaceFormData {

        typedef FaceFormData Base;

        static dealii::UpdateFlags update_flags()
        {
            using namespace dealii;
            return update_values|update_quadrature_points|update_JxW_values;
        }

        // The values of the integrand are the values of one solution component.
        std::vector<double> itg_values;

        FluxForm(const dealii::Mapping<dim> & mapping,
                 const dealii::FiniteElement<dim> & fe,
                 const dealii::Quadrature<dim-1> & q_rule)
            :
              Base(mapping, fe, q_rule, update_flags() ),
              itg_values(this->n_q_points)
        {}


        // Besides collecting all the convenience variables we still provide a function
        // doing all the reinitializations needed on a cell.
        template<typename cellIterator>
        inline double operator() (const cellIterator & cell,
                                const uint face_id,
                                const dealii::Vector<double>& integrand)
        {
            this->Base::FeV::reinit(cell, face_id);

            this->Base::FeV::get_function_values(integrand, itg_values);

            double cell_integral = 0.;

            for (uint a=0; a < this->n_q_points; ++a)
                cell_integral += itg_values[a] * this->JxW[a];

            return cell_integral;
        }
    };




    struct PureVNBoundaryForm : public step27::FormData<dim, step27::FaceIntTraits<dim> > {

        typedef step27::FormData<dim, step27::FaceIntTraits<dim> > Base;

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

        // The von Neumann values are the gradients.
        std::vector<dealii::Tensor<1,dim> >vN_values;

        dealii::Vector<double>       cell_rhs;

        const dealii::SmartPointer<const dealii::Function<dim> > bc;

        PureVNBoundaryForm(const dealii::Mapping<dim> & mapping,
                           const dealii::FiniteElement<dim> & fe,
                           const dealii::Quadrature<dim-1> & q_rule,
                           const dealii::Function<dim> & bc_eval)
            :
              Base(mapping, fe, q_rule, update_flags() ),
              normals (this->get_normal_vectors () ),
              vN_values(this->n_q_points),
              cell_rhs(this->dofs_per_cell),
              bc(&bc_eval)
        {}


        // Besides collecting all the convenience variables we still provide a function
        // doing all the reinitializations needed on a cell.
        template<typename cellIterator>
        void assemble (const cellIterator & cell,
                       const uint face_id);
    };

    // surface curl of trial function:
    // this->normals[a] "X" this->shape_grad(i, a) * this->normals[b] "X" this->shape_grad(j, b)

protected:

    //void assemble_system (uint cycle);

    void assemble_system_generic (uint cycle, step27::Architecture);

    void assemble_boundary_mass_matrices();

    void reassemble_nonlinearity_cpu ();



    void assemble_system_multithreaded (uint cycle);
    void assemble_multigrid ();

    void assemble_bem_unit_test();

    void setup_preconditioner();


    void compute_bem_assembly_errors();


    void refine_grid_global ();
    void refine_grid_adaptive ();

    // For debugging purposes we add the option to give an additional counter
    // for the individual iteration step.
    void output_results (const unsigned int cycle, const int iteration=-1) const;



    //Function to generate matrix values, to be run in a thread
    void generate_bem_matrix_values(const unsigned int n_bem_points, const unsigned int n_bem_quad_points,
                                    const unsigned int dofs_per_be,
                                    const std::vector<dealii::Point<dim> > x, const std::vector<dealii::Point<dim> > q,
                                    const dealii::FullMatrix<double> W,
                                    const std::vector<unsigned int> global_bc_dofs,
                                    const std::vector<double> normals,
                                    const std::vector<dealii::Tensor<1, dim> > u_gradients, const std::vector<double> JxW);
    // Copy std::vector<Point> to std::vector<double>


    step27::NumericsParams<dim> & numerics_prms;

    dealii::ConvergenceTable     timing_table;

    dealii::Triangulation<dim>   triangulation;

    dealii::MappingQ<dim>        mapping;

    dealii::FE_Q<dim>            fe; // for DG this is the space Q

    // FE_DGQ<dim> fe_dgq;

    dealii::MGDoFHandler<dim>    mg_dof_handler;


    // The boundary elements for the interface
    // need a different quadrature rule which only integrates over surfaces and not over volumes.
    // The matrix entries are assembled by the BEMForm. In this step it is based on collocation.
  typedef step27::BemCollocationForm<dim, step27::FaceIntTraits<dim>, step28::CUDADriver<dim> > BemCollForm;

   BemCollForm * bem_form;


    // Similar to the @p SolutionSet we pool all sparsity patterns
    // of the FEM part of the DRS problem in a structure in order
    // to centralize and simplify the memory management.

    // For simplified access we define a list of Ids for the different structures which we will encounter.



    typedef step27::PatternInfo PatternInfo;

    // TO DO: subclass step27::SparsityPatternSet and add boundary mass matrices in the two flavors that occur:
    // - as sparse-row-sparse-column matrix (n x n)
    // - as sparse boundary mass matrix (m x m)
    // (somehow done)


    // Unlike in the simple FEM-BEM problem we need additional mass matrices
    // for the flux boundary conditions for the ions.
    // Therefore we subclass the SparsityPatternSet and overwrite the reinit_FEM function of the base class.
    struct SparsityPatternSet28 : public step27::SparsityPatternSet<PatternInfo>
    {
        typedef step27::SparsityPatternSet<PatternInfo> Base;

        SparsityPatternSet28() : Base() {}

        void reinit_FEM(const dealii::DoFHandler<dim> & dof_handler)
        {
            this->Base::reinit_FEM(dof_handler);

            // kill any old content. This is more bullet-proof than calling clear().
            {
                std::map<dealii::types::boundary_id, std::map<uint, uint> >  tmp;
                std::swap(tmp, electrode_bc_local_2_global_map);
            }

            // We need the matrices once as sparse-row-sparse-column matrix (n x n)
            dealii::SparsityPattern & cathode_mass = (*this)[PatternInfo::FEMCathodeBoundaryMass];
            dealii::SparsityPattern & anode_mass = (*this)[PatternInfo::FEMAnodeBoundaryMass];

            // and once as sparse boundary mass matrix (m x m)
            dealii::SparsityPattern & cathode_mass_bc_dofs = (*this)[PatternInfo::BEMCathodeBoundaryMass];
            dealii::SparsityPattern & anode_mass_bc_dofs = (*this)[PatternInfo::BEMAnodeBoundaryMass];

            build_boundary_pattern(cathode_mass, cathode_mass_bc_dofs, dof_handler, DRS::BCids::cathode() );

            build_boundary_pattern(anode_mass, anode_mass_bc_dofs, dof_handler, DRS::BCids::anode() );
        }
    };

    typedef SparsityPatternSet28 PatternSet;
    // typedef step27::SparsityPatternSet<PatternInfo> PatternSet;

    PatternSet sparsity_patterns;

    // @sect3{Class: MtxInfo}
    // We need additional matrices for the flux terms on the sub-boundaries representing the electrodes.
    // Therefore, we subclass step27::MtxInfo and add the necessary matrices.
    struct MtxInfo : public step27::MtxInfo {

        typedef step27::MtxInfo Base;

        typedef step27::MtxInfo::Ids Ids;

        typedef Base::SPInfo SPInfo;

        MtxInfo()
            :
              Base()
        {
            // Extend the pattern map of the base class
            // At the electrodes we need FE mass matrices for the boundary conditions.
            // These mass matrices are formed by the surface DoFs. And therefore are rather BEM matrices (TODO: rename enums?)
            pattern_types[CathodeBoundaryMass] = PatternInfo::FEMCathodeBoundaryMass;
            pattern_types[AnodeBoundaryMass] = PatternInfo::FEMAnodeBoundaryMass;
            // The drift terms in the Poisson-Nernst-Planck equations are part of the FEM problem
            // and do not have any Dirichlet boundary conditions.
            pattern_types[DriftCations] = PatternInfo::FEMwvNBC;
            pattern_types[DriftAnions] = PatternInfo::FEMwvNBC;
            pattern_types[DriftPot] = PatternInfo::FEMwvNBC;



            qDebug() << "n sparsity patterns : " << pattern_types.size();
        }
    };


    //typedef step27::MtxInfo MtxInfo27;

   // typedef MtxInfo27::Ids MatrixIds28;

    // SparseMatrixSet28 matrices;
    // For the set of sparse matrices of the blocks of the PDE system we reuse the one from step-27.
    typedef step27::SparseMatrixSet<MtxInfo> SparseMatrixSet28;
    SparseMatrixSet28 matrices;




    step27::BlockPatternInfo drs_block_pattern;  

    step27::BlockPatternInfo drs_pc_block_pattern;

    // We need an additional object for the
    // hanging nodes constraints. They are
    // handed to the transfer object in the
    // multigrid. Since we call a compress
    // inside the multigrid these constraints
    // are not allowed to be inhomogeneous so
    // we store them in different ConstraintMatrix
    // objects.
    dealii::ConstraintMatrix     hanging_node_constraints;
    dealii::ConstraintMatrix     constraints;

    dealii::ConstraintMatrix     constraints_fem_bem_unit_test;


    typedef step27::SolComps SolComps;
   // enum SolComps {Cations=1, Anions=2, Neutral=3, Phi=4, Phi_derivative=5 };

    // static const SolComps drs_sol_comps[];
    static const int n_drs_components = 5;
    const std::vector<SolComps> drs_sol_comps;

    const std::vector<SolComps> drs_sol_comps_fem_part;

    const std::vector<SolComps> drs_sol_comps_fem_part_all_vN;

    const std::vector<SolComps> drs_sol_comps_bem_part;



    // We need quite few solution vectors. To simplify administrative things
    // like resizing we create a little structure which does all these things
    // and which holds a few references so that we can address the individual
    // solution components by name. We implement it based on a container
    // from the standard library.
    typedef step27::SolutionSet SolutionSet;


    SolutionSet old_solution;
    SolutionSet new_solution;

    SolutionSet vN_unit_test_new_solution;
    SolutionSet vN_unit_test_old_solution;

    dealii::Vector<double>       system_rhs;

    dealii::Vector<double> system_rhs_nl_contrib_cats;
    dealii::Vector<double> system_rhs_nl_contrib_ans;

    dealii::Vector<double>       system_rhs_pure_vN_bc;




    struct LinearSystem {

        // Before we set up the blocked system of algebraic equations we define
        // some named constants for the block positions of the individual components.
        // In case of a negative index the corresponding component is omitted.
        // This is primarily for testing and development.

        // mapping of physical quantities to rows of blocks


//        static const int fem_pot_idx = 0;
//        static const int bem_pot_idx = 1;

//        static const int cats_idx = 2; // 2;

//        static const int neut_idx = 3;

//        static const int ans_idx = 4;

// TODO: make non-public
         std::map<SolComps, uint> component_names2index;

        const unsigned int n_blocks() const { return component_names2index.size(); }

        const uint block_id(const SolComps c) const
        {
            std::map<SolComps, uint>::const_iterator b = component_names2index.find(c);

            if (b != component_names2index.end())
                return b->second;
            else
                Assert(false, dealii::ExcMessage("unknown solution component. Did you already pass them to the linear system?"));
        }


        dealii::BlockMatrixArray<double> A;



        dealii::BlockMatrixArray<double> A_nl;

        dealii::BlockMatrixArray<double> A_2dt;

        dealii::PreconditionIdentity bem_preconditioner;


        LinearSystem(/*TODO: pass block_pattern and put setup_block_pattern into CTOR*/) : block_sizes_fem_bem_unit_test(2) {}

        void clear() {
            std::cout << __FUNCTION__ << std::endl;
            A.clear();
            A_nl.clear();
            A_2dt.clear();
            block_prec.clear();
        }

        ~LinearSystem() {

            this->clear();

            std::cout << __FUNCTION__ << " DONE" << std::endl;
        }

        std::vector<unsigned int> block_sizes, block_sizes_fem_bem_unit_test;


        dealii::BlockVector<double> rhs;


        dealii::BlockVector<double> sol;
        dealii::BlockVector<double> old_sol;


        dealii::BlockTrianglePrecondition<double> block_prec;


        dealii::BlockVector<double> rhs_fem_bem_unit_test;
        dealii::BlockVector<double> sol_fem_bem_unit_test;


        void setup_block_system(const SparseMatrixSet28 & matrices,
                                const step27::BlockPatternInfo& block_pattern,
                                const SolutionSet & sol_set
                                );

        void setup_pc (const step27::BlockPatternInfo & pc_block_pattern,
                       const step27::GeoMGPC<dim> & mg_pc_Phi,
                       const step27::GeoMGPC<dim> & mg_pc_Ions,
                       const SparseMatrixSet28 & matrices);

        void reinit (const unsigned int n_fem_dofs,
                     const unsigned int n_bem_dofs)
        {
            __n_fem_dofs = n_fem_dofs;
            __n_bem_dofs = n_bem_dofs;

            {
                std::vector<unsigned int> tmp;
                std::swap(tmp, block_sizes);
            }


            block_sizes.resize(n_blocks(), n_fem_dofs);

            block_sizes_fem_bem_unit_test[0] = n_fem_dofs;
            block_sizes_fem_bem_unit_test[1] = n_bem_dofs;


            if (this->component_names2index.find(step27::Phi_derivative) != this->component_names2index.end() )
                block_sizes[  this->component_names2index[step27::Phi_derivative]  ] = n_bem_dofs;


            rhs.reinit(block_sizes);
            sol.reinit(block_sizes);
            old_sol.reinit(block_sizes);

            rhs_fem_bem_unit_test.reinit(block_sizes_fem_bem_unit_test);
            sol_fem_bem_unit_test.reinit(block_sizes_fem_bem_unit_test);
        }

        void set_solution(const step27::SolutionSet& s_set)
        {
             std::map<SolComps, uint>::const_iterator
                     s = component_names2index.begin(),
                     end_s = component_names2index.end();

             for (; s != end_s; ++s)
             {
                 const dealii::Vector<double> & src = s_set(s->first); // access by name
                 dealii::Vector<double> & dst  = sol.block(s->second);  // access by index

                 if (dst.size() == src.size())
                     dst = src;
                 else {
                    // qDebug() << "size mismatch in comp" << s->second;
                   const dealii::VectorSlice<const dealii::Vector<double> > slice(src, (__n_fem_dofs - __n_bem_dofs), __n_bem_dofs );
                     std::copy(slice.begin(), slice.end(), dst.begin());
                 }
             }

        }

        void compute_incr_and_push_sol_to( std::map<SolComps, double> & l2_incr, step27::SolutionSet& s_set)
        {
            std::map<SolComps, uint>::const_iterator
                    s = component_names2index.begin(),
                    end_s = component_names2index.end();

            for (; s != end_s; ++s)
            {
                dealii::Vector<double> & dst = s_set(s->first); // access by name
                const dealii::Vector<double> & src  = sol.block(s->second);  // access by index

                if (dst.size() == src.size())
                {
                    dst -= src;
                    l2_incr[s->first] = dst.l2_norm();
                    dst = src;
                }
                else {
                    // qDebug() << "size mismatch in comp" << s->second;
                    dealii::VectorSlice<dealii::Vector<double> > slice(dst, (__n_fem_dofs - __n_bem_dofs), __n_bem_dofs );

                    double l2_incr_sqr = 0.;
                    for (uint i = 0; i < slice.size(); i++)
                    {
                        slice[i] -= src(i);
                    l2_incr_sqr += slice[i]*slice[i];
                    }

                    l2_incr[s->first] = std::sqrt(l2_incr_sqr);
                    std::copy(src.begin(), src.end(), slice.begin());
                }
            }

        }


        void push_sol_to(step27::SolutionSet& s_set)
        {
            std::map<SolComps, uint>::const_iterator
                    s = component_names2index.begin(),
                    end_s = component_names2index.end();

            for (; s != end_s; ++s)
            {
                dealii::Vector<double> & dst = s_set(s->first); // access by name
                const dealii::Vector<double> & src  = sol.block(s->second);  // access by index

                if (dst.size() == src.size())
                {
                    dst = src;
                }
                else {
                    // qDebug() << "size mismatch in comp" << s->second;
                    dealii::VectorSlice<dealii::Vector<double> > slice(dst, (__n_fem_dofs - __n_bem_dofs), __n_bem_dofs );
                    std::copy(src.begin(), src.end(), slice.begin());
                }
            }

        }



        void solve();

        // Temporary functions until BlockVectors can be fed into SolutionTransfer
        void copy_sol_set_into_block_vec()
        {


        }
private:
        uint __n_fem_dofs;
    uint __n_bem_dofs;
    };
    // Put fem_bem data before lsys, because lsys connects its internal matrices to the fem_bem SL matrix,
    // so fem_bem must be destroyed after lsys! Otherwise, dealii catches fire.

    FemBemData fem_bem;

    LinearSystem lsys;


    const unsigned int degree;

    // Although ther different problems can share their sparsity patterns
    // we need separate matrices because the potential has Dirichlet boundary
    // conditions while the ions have not.
    // Therefore, we need two independent data sets and eventually preconditioners.
    step27::MGData<dim> mg_Phi;
    step27::MGData<dim> mg_Ions;


    step27::PhysicalParams phys_data;

    step27::Architecture arch;

    step27::SubModel sub_model;

    mutable dealii::ConvergenceTable unit_test_table;

    mutable dealii::ConvergenceTable measured_currents_table;


    void setup_drs_block_pattern();

    void setup_fem_bem_unit_test_block_pattern();
protected:
    void compute_electrode_fluxes() const;
};


// template <int dim>
// const step27::SolComps DRSProblem<dim>::drs_sol_comps[] = { step27::Phi, step27::Phi_derivative, step27::Cations, step27::Anions, step27::Neutral };




// The assemble() function assembles the contribution to the
// left- and right-hand side on a cell.
template <int dim>
template<typename cellIterator>
void DRSProblem<dim>::DriftForm::assemble(const cellIterator & cell) {

    this->Base::reinit(cell);

    this->Base::FeV::get_function_values(cation_density, cation_density_local_values);
    this->Base::FeV::get_function_values(anion_density, anion_density_local_values);

    this->Base::FeV::get_function_gradients(potential, Phi_local_values);

    cell_matrix_cations = 0.;
    cell_matrix_anions = 0.;
    cell_matrix_drift = 0.;
    cell_rhs = 0;

    for (uint q_pt=0; q_pt<this->n_q_points; ++q_pt)
    {
        // assemble Laplacian
        for (uint i=0; i<this->dofs_per_cell; ++i)
        {
            double shape_grad_dot_Phi_grad_JxW_q
                    =
                    this->shape_grad(i,q_pt) * this->Phi_local_values[q_pt]
                    *
                    this->JxW[q_pt];

            for (uint j=0; j<this->dofs_per_cell; ++j)
            {
                double grad_dot_grad_JxW_q
                        =
                        (this->shape_grad(i,q_pt) * this->shape_grad(j,q_pt) )
                        *
                        this->JxW[q_pt];

                cell_matrix_cations(i,j)
                        +=
                        cation_density_local_values[q_pt] * grad_dot_grad_JxW_q;


                cell_matrix_anions(i,j)
                        // anions move in the opposite direction, therefore: MINUS!
                        -=
                        anion_density_local_values[q_pt] * grad_dot_grad_JxW_q;

                cell_matrix_drift(i,j) +=
                        shape_grad_dot_Phi_grad_JxW_q
                        +
                        this->shape_value(j, q_pt);

            }

            // assemble nonlinearity
        }
    }
}


template <int dim>
void DRSProblem<dim>::DriftForm::distribute_local_to_global(const dealii::ConstraintMatrix & constraints,
                                                            dealii::SparseMatrix<double> & cat_drift_matrix,
                                                            dealii::SparseMatrix<double> & an_drift_matrix,
                                                            dealii::SparseMatrix<double> & pot_drift_matrix)
{

    constraints.distribute_local_to_global (cell_matrix_cations,
                                            this->local_dof_indices,
                                            cat_drift_matrix);

    constraints.distribute_local_to_global (cell_matrix_anions,
                                            this->local_dof_indices,
                                            an_drift_matrix);

    constraints.distribute_local_to_global (cell_matrix_drift,
                                            this->local_dof_indices,
                                            pot_drift_matrix);
}





// Besides collecting all the convenience variables we still provide a function
// doing all the reinitializations needed on a cell.
template <int dim>
template<typename cellIterator>
void DRSProblem<dim>::PureVNBoundaryForm::assemble (const cellIterator & cell,
                                                    const uint face_id)
{
    this->Base::FeV::reinit(cell, face_id);

    cell->get_dof_indices (this->local_dof_indices);

    this->bc->gradient_list(this->q, vN_values);

    for (uint a=0; a < this->n_q_points; ++a)
        for (uint j = 0; j <this->dofs_per_cell; ++j)
            cell_rhs(j) += this->shape_value(j, a)
                    *
                    (this->normals[a]
                     *
                     this->vN_values[a])
                    *
                    this->JxW[a];
}



// @sect3{The <code>DRSProblem</code> class implementation}

// @sect4{DRSProblem::DRSProblem}

// The constructor is left mostly
// unchanged. We take the polynomial degree
// of the finite elements to be used as a
// constructor argument and store it in a
// member variable.
//
// By convention, all adaptively refined
// triangulations in deal.II never change by
// more than one level across a face between
// cells. For our multigrid algorithms,
// however, we need a slightly stricter
// guarantee, namely that the mesh also does
// not change by more than refinement level
// across vertices that might connect two
// cells. In other words, we must prevent the
// following situation:
//
// @image html limit_level_difference_at_vertices.png ""
//
// This is achieved by passing the
// Triangulation::limit_level_difference_at_vertices
// flag to the constructor of the
// triangulation class.
template <int dim_>
DRSProblem<dim_>::DRSProblem (step27::NumericsParams<dim_> &numerics_prms,
                              const std::string &prm_path,
                              const QDir &log_dir,
                              step27::TestExpr testcase,
                              step27::SubModel sm)
    :
      numerics_prms(numerics_prms),
      triangulation (dealii::Triangulation<dim_>::
                     limit_level_difference_at_vertices),
      mapping (numerics_prms.fe_mapping_degree),
      fe (numerics_prms.fe_degree),
      mg_dof_handler (triangulation),
      old_solution(),
      new_solution(),
      degree (numerics_prms.fe_degree),
      arch(numerics_prms.arch),
      sub_model(sm),
      fem_bem(testcase, prm_path, log_dir),
      matrices(sparsity_patterns)
{
    dealii::ParameterHandler prm_handler;
    step27::PhysicalParams::declare(prm_handler);

    std::string phys_prm_path(prm_path + "/physical_properties.prm");
    prm_handler.read_input(phys_prm_path);

    phys_data.get(prm_handler);

    std::ofstream log_phys_prop_prm( (log_dir.absolutePath() + QDir::separator() + "physical_properties.prm.log").toStdString().c_str() );
    prm_handler.print_parameters (log_phys_prop_prm,
                                  dealii::ParameterHandler::Text);

    // setup boundary element and form evaluation

    bem_form = new BemCollForm  (DRS::BCids::cavity(), // TO DO: pass id of a dielectric interface, not of a particular one
                                 mapping,
                                 fe,
                                 numerics_prms);

    // Before we can set up any vectors or function spaces we have to instantiate
    // the skeleton of our PDE system.
    // The block structure of the PDE system does not change over the course of the simulation.
    // Hence we instantiate it once and for all in the constructor.
    switch (this->sub_model) {
    case step27::fem_bem_unit_test:
        setup_fem_bem_unit_test_block_pattern();
        break;
    case step27::full_DRS:
        setup_drs_block_pattern();
        break;
    case step27::pure_von_Neumann:
        AssertThrow(false, dealii::ExcNotImplemented() );
    default:
        AssertThrow(false, dealii::ExcMessage("Unknown submodel specified. "
                                              "Please select either fem_bem_unit_test or full_DRS. "
                                              "Check your parameter file."));
        break;
    }

    std::vector<uint> block_sizes(n_drs_components, 1); // Dummy sizes.
    old_solution.reinit(drs_sol_comps, block_sizes);

    new_solution.reinit(drs_sol_comps, block_sizes);

   // vN_unit_test_new_solution.reinit(step27::vN_sol_comps, std::vector<uint>(1,1));
   // vN_unit_test_old_solution.reinit(step27::vN_sol_comps, std::vector<uint>(1,1));



    // Initialize the simulation
    {

        DRS::MeshGenerator::DRSCellParams drs_cell_params;

        // Before we set up the triangulation we have to read the parameters
        // which determine its size.
        dealii::ParameterHandler prm_handler;

        DRS::MeshGenerator::DRSCellParams::declare(prm_handler);
        prm_handler.read_input(numerics_prms.prm_drs_cell_filepath);

        drs_cell_params.get(prm_handler);


        std::ofstream log_drs_cell_prm( (log_dir.absolutePath() + QDir::separator() + QFileInfo(numerics_prms.prm_drs_cell_filepath.c_str()).fileName() + ".log").toStdString().c_str() );
        prm_handler.print_parameters (log_drs_cell_prm,
                                      dealii::ParameterHandler::Text);

        DRS::MeshGenerator::drs_cell(triangulation, drs_cell_params, //numerics_prms.prm_drs_cell_filepath,
                                     numerics_prms.n_init_refinements);


        // Prepare finite element space on the initial coarse grid.
        this->mg_dof_handler.distribute_dofs (fe);

        setup_discretization ();

        lsys.setup_block_system(this->matrices, drs_block_pattern, new_solution);


        lsys.reinit(fem_bem.n_fem_dofs, bem_form->n_bc_dofs);
    }

}


template<int dim>
void DRSProblem<dim>::setup_drs_block_pattern() {

    const double EPS_S = phys_data.EPSILON_S;
    const double EPS_P = phys_data.EPSILON_P;

    const double K_O = phys_data.k_O;

    const double K_R = phys_data.k_R;

   // const double inv_2dt = 2. / dt;

    // First of all, setup the vector with the names of the components
   const SolComps sol_comps [] = { step27::Phi, step27::Phi_derivative, step27::Cations, step27::Anions, step27::Neutral };


   const_cast<std::vector<SolComps>& >(drs_sol_comps).assign(sol_comps, sol_comps + 5);

   const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part_all_vN).assign(sol_comps+2, sol_comps + 5);

   const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part).resize(4);

   const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part)[0] = sol_comps [0];

   std::copy(sol_comps+2, sol_comps+5,  const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part).begin()+1);

   const_cast<std::vector<SolComps>& >(drs_sol_comps_bem_part).assign(sol_comps+1, sol_comps + 2);

    // Set up of the block patterns of the linear systems

    // diagonal, bulk
    SolComps ions [] = { step27::Cations, step27::Neutral, step27::Anions };

    for (int c = 0; c < 3; c++)
        drs_block_pattern[step27::IJ(ions[c], ions[c])].push_back({ MtxInfo::LaplaceIons, 1., MtxInfo::LaplaceAllvNPC, 1. });

    drs_block_pattern[step27::IJ(step27::Phi, step27::Phi)].push_back({ MtxInfo::LaplacePhi, 1., MtxInfo::LaplaceDRSPhiPC, 1. });
    drs_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi_derivative)].push_back({ MtxInfo::SLFullMatrix, +1.*(EPS_S/EPS_P), MtxInfo::Identity, 1. });

    // off-diagonals
    drs_block_pattern[step27::IJ(step27::Cations, step27::Phi)].push_back({ MtxInfo::DriftCations, 1. });
    drs_block_pattern[step27::IJ(step27::Phi, step27::Cations)].push_back({ MtxInfo::BulkMass, -1./EPS_S });

    drs_block_pattern[step27::IJ(step27::Anions, step27::Phi)].push_back({ MtxInfo::DriftAnions, 1. });
    drs_block_pattern[step27::IJ(step27::Phi, step27::Anions)].push_back({ MtxInfo::BulkMass, +1./EPS_S });

    // boundary terms
    drs_block_pattern[step27::IJ(step27::Cations, step27::Cations)].push_back({ MtxInfo::CathodeBoundaryMass, +K_R });
    drs_block_pattern[step27::IJ(step27::Neutral, step27::Cations)].push_back({ MtxInfo::CathodeBoundaryMass, -K_R });

    drs_block_pattern[step27::IJ(step27::Cations, step27::Neutral)].push_back({ MtxInfo::AnodeBoundaryMass, -K_O });
    drs_block_pattern[step27::IJ(step27::Neutral, step27::Neutral)].push_back({ MtxInfo::AnodeBoundaryMass, +K_O });

    drs_block_pattern[step27::IJ(step27::Phi, step27::Phi_derivative)].push_back({MtxInfo::FEMBEMUpperOffBlock, -1.*(EPS_P/EPS_S)});
    drs_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi)].push_back({MtxInfo::BEMFEMLowerOffBlock, +1., MtxInfo::BEMFEMLowerOffBlock, 1. });

#ifndef nUSE_PC_PATTERN
    // For the ions we can use the generic MG preconditioner which uses von Neumann boundary conditions on all boundaries.
    for (int c = 0; c < 3; c++)
        drs_pc_block_pattern[step27::IJ(ions[c], ions[c])].push_back({MtxInfo::LaplaceAllvNPC, 1. });

    // For the potential we need a MG preconditioner with problem-specific boundary conditions.
    drs_pc_block_pattern[step27::IJ(step27::Phi, step27::Phi)].push_back({ MtxInfo::LaplaceDRSPhiPC, 1. });

    // for the time being the BEM part is preconditioned with the identity matrix.
    drs_pc_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi_derivative)].push_back({ MtxInfo::Identity, 1. });

    // Now, do the off-diagonals
    drs_pc_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi)].push_back({ MtxInfo::BEMFEMLowerOffBlock, 1. });

#endif
}



template<int dim>
void DRSProblem<dim>::setup_fem_bem_unit_test_block_pattern() {

    const double EPS_S = 1; // phys_data.EPSILON_S;
    const double EPS_P = 1; // phys_data.EPSILON_P;

    const double K_O = 0; // phys_data.k_O;

    const double K_R = 0; // phys_data.k_R;

   // const double inv_2dt = 2. / dt;

    // First of all, setup the vector with the names of the components
   const SolComps sol_comps [] = { step27::Phi, step27::Phi_derivative };


   const_cast<std::vector<SolComps>& >(drs_sol_comps).assign(sol_comps, sol_comps + 2);



   const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part).resize(1);

   const_cast<std::vector<SolComps>& >(drs_sol_comps_fem_part)[0] = sol_comps [0];

   const_cast<std::vector<SolComps>& >(drs_sol_comps_bem_part).assign(sol_comps+1, sol_comps + 2);

    // Set up of the block patterns of the linear systems
    drs_block_pattern[step27::IJ(step27::Phi, step27::Phi)].push_back({ MtxInfo::LaplacePhi, 1., MtxInfo::LaplaceDRSPhiPC, 1. });
    drs_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi_derivative)].push_back({ MtxInfo::SLFullMatrix, +1.*(EPS_S/EPS_P), MtxInfo::Identity, 1. });

    // boundary terms
    drs_block_pattern[step27::IJ(step27::Phi, step27::Phi_derivative)].push_back({MtxInfo::FEMBEMUpperOffBlock, -1.*(EPS_P/EPS_S)});
    drs_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi)].push_back({MtxInfo::BEMFEMLowerOffBlock, +1., MtxInfo::BEMFEMLowerOffBlock, 1. });

#ifndef nUSE_PC_PATTERN
    // For the potential we need a MG preconditioner with problem-specific boundary conditions.
    drs_pc_block_pattern[step27::IJ(step27::Phi, step27::Phi)].push_back({ MtxInfo::LaplaceDRSPhiPC, 1. });

    // for the time being the BEM part is preconditioned with the identity matrix.
    drs_pc_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi_derivative)].push_back({ MtxInfo::Identity, 1. });

    // Now, do the off-diagonals
    drs_pc_block_pattern[step27::IJ(step27::Phi_derivative, step27::Phi)].push_back({ MtxInfo::BEMFEMLowerOffBlock, 1. });

#endif
}





// @sect4{DRSProblem::setup_discretization}

// The following function extends what
// setup_system one in deal.II's step-16 did.
// We setup the finite element space, extrat the dofs for the boundary elements,
// create the sparsity patterns, the MG preconditioners and some initial values.
// This function has to be called before the linera system is set up, otherwise
// the matrix set is empty and the linear system cannot be built.
template <int dim>
void DRSProblem<dim>::setup_discretization ()
{
    // Here we output not only the
    // degrees of freedom on the finest
    // level, but also in the
    // multilevel structure
    dealii::deallog << "Number of degrees of freedom: "
                    << mg_dof_handler.n_dofs();

    dealii::DoFHandler<dim> & dh = static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler);

    // Before we deal with the dof enumeration of FEM-BEM coupling
    // we cope with the dofs on the electrodes.

    std::set< dealii::types::boundary_id > cathode_boundary_indicator;
    std::vector<bool> cathode_dof_flags;

    std::set< dealii::types::boundary_id > anode_boundary_indicator;
    std::vector<bool> anode_dof_flags;

    cathode_boundary_indicator.insert(DRS::BCids::cathode());
    anode_boundary_indicator.insert(DRS::BCids::anode());
    /* this is probably useful in order to get from the small boundary mass matrix back to the large global matrix
 std::vector< unsigned int >	cathode_mapping;
    std::vector< unsigned int >	anode_mapping;
 DoFTools::map_dof_to_boundary_indices(mg_dof_handler,
                                          cathode_boundary_indicator,
                                          cathode_mapping
                                          );*/

    // We have only one component. So the following is trivial but necessary.
    std::vector< bool > component_select(1 /*n_components*/);
    component_select[0] = true;

    cathode_dof_flags.resize(mg_dof_handler.n_dofs());
    dealii::DoFTools::extract_boundary_dofs	(dh,
                                             component_select,
                                             cathode_dof_flags,
                                             cathode_boundary_indicator);
    // Now sort back. At the end ALL boundary dofs will be selected back but their relative order
    // is given by the order in which we tackle the different subboundaries.
    dealii::DoFRenumbering::sort_selected_dofs_back (dh,
                                                     cathode_dof_flags);

    anode_dof_flags.resize(mg_dof_handler.n_dofs());
    dealii::DoFTools::extract_boundary_dofs	(dh,
                                             component_select,
                                             anode_dof_flags,
                                             anode_boundary_indicator);

    dealii::DoFRenumbering::sort_selected_dofs_back (static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                     anode_dof_flags);


    //! For the FEM-BEM coupling we sort the dofs of the potential
    //! on the surface of the cavity to the back.
    //! This leads to a dense block in the lower right of the stiffness matrix.

    // Only after the final dof numbering scheme is established
    // we can setup the BEM problem
    bem_form->reinit_bem_data(mg_dof_handler, this->mapping);

    fem_bem.reinit(DRS::BCids::cavity(), mg_dof_handler);

    // After the FEM do the sparsity patterns
    // for the off-diagonal coupling matrices
    {
        // Using operator[] automatically adds an entry to the underlying std::map of the collection
        // of sparsity patterns, if it does not exist, yet.
        this->sparsity_patterns[PatternInfo::FEMBEM].reinit(bem_form->n_bc_dofs,
                                                            bem_form->n_bc_dofs,
                                                            mg_dof_handler.max_couplings_between_dofs());


        dealii::DoFTools::make_boundary_sparsity_pattern(
                    static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                    bem_form->interface,
                    bem_form->dof_to_boundary_mapping,
                    this->sparsity_patterns(PatternInfo::FEMBEM)
                    );

        this->sparsity_patterns(PatternInfo::FEMBEM).compress();

        sparsity_patterns.reinit_FEM_BEM(mg_dof_handler, fem_bem);



    }

    // Only after all the dof sorting is done we can setup the sparsity patterns.
    // In addition to the standard setup of the FEM part we need boundary mass matrices
    // for the Robin-like boundary conditions for the ion densities.
    sparsity_patterns.reinit_FEM (mg_dof_handler);



    for (unsigned int l=0;l<triangulation.n_levels();++l)
        dealii::deallog << "   " << 'L' << l << ": "
                        << mg_dof_handler.n_dofs(l);
    dealii::deallog  << std::endl;



    system_rhs.reinit (mg_dof_handler.n_dofs());
    system_rhs_pure_vN_bc.reinit (mg_dof_handler.n_dofs());

    system_rhs_nl_contrib_cats.reinit (mg_dof_handler.n_dofs());
    system_rhs_nl_contrib_ans.reinit (mg_dof_handler.n_dofs());

    //system_rhs_pure_bem_unit_test.reinit(mg_dof_handler.n_boundary_dofs());


    // But it starts to be a wee bit different
    // here, although this still doesn't have
    // anything to do with multigrid
    // methods. step-6 took care of boundary
    // values and hanging nodes in a separate
    // step after assembling the global matrix
    // from local contributions. This works,
    // but the same can be done in a slightly
    // simpler way if we already take care of
    // these constraints at the time of copying
    // local contributions into the global
    // matrix. To this end, we here do not just
    // compute the constraints due to hanging
    // nodes, but also due to zero boundary
    // conditions. We will
    // use this set of constraints later on to
    // help us copy local contributions
    // correctly into the global linear system
    // right away, without the need for a later
    // clean-up stage:
    constraints.clear ();
    hanging_node_constraints.clear ();
    dealii::DoFTools::make_hanging_node_constraints (mg_dof_handler, hanging_node_constraints);
    dealii::DoFTools::make_hanging_node_constraints (mg_dof_handler, constraints);


    dealii::VectorTools::interpolate(mapping, static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                     *fem_bem.u_ref, fem_bem.solution_ref);

    typename dealii::FunctionMap<dim>::type       Phi_dirichlet_boundary;



    dealii::VectorTools::interpolate(mapping, static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                     *fem_bem.newton_potential,
                                     fem_bem.newton_potential_rhs);

    fem_bem.newton_potential_rhs *= 4 * M_PI;

// #ifndef FEM_BEM_UNIT_TEST
    if (this->sub_model == step27::fem_bem_unit_test)
    {
        typename dealii::FunctionMap<dim>::type & Phi_dirichlet_boundary_unit_test
                = Phi_dirichlet_boundary;

        Phi_dirichlet_boundary_unit_test[DRS::BCids::anode()]   = fem_bem.u_ref; //&harmonic_test_bc; // &homogeneous_dirichlet_bc;
        Phi_dirichlet_boundary_unit_test[DRS::BCids::cathode()] = fem_bem.u_ref; //&harmonic_test_bc; // &cathode_dirichlet_bc;

        // For debugging purposes it is helpful to have the exact solution also on the cavity surface as BC.
        // dirichlet_boundary[DRS::BCids::cavity()] = fem_bem.u_ref; // &harmonic_test_bc; // &protein_surface_bc;

        Phi_dirichlet_boundary_unit_test[DRS::BCids::cylinder_surface()] = fem_bem.u_ref; // &harmonic_test_bc;

        dealii:: VectorTools::interpolate_boundary_values (mapping,
                                                           static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                           Phi_dirichlet_boundary_unit_test,
                                                           constraints/*_fem_bem_unit_test*/);
    }
    else
    {
        // The potential has non-zero values at the electrodes
        // but otherwise it is zero.

        dealii::ConstantFunction<dim>           cathode_dirichlet_bc(
                    + phys_data.potential_difference
                    , 1);

        dealii::ConstantFunction<dim>           anode_dirichlet_bc(
                    - phys_data.potential_difference
                    , 1);

        Phi_dirichlet_boundary[DRS::BCids::cathode()] = &cathode_dirichlet_bc;
        Phi_dirichlet_boundary[DRS::BCids::anode()] = &anode_dirichlet_bc;

        dealii::VectorTools::interpolate_boundary_values (mapping,
                                                          static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                          Phi_dirichlet_boundary,
                                                          constraints);
    }
    constraints.close ();

    hanging_node_constraints.close ();

    hanging_node_constraints.distribute (fem_bem.solution_ref);

    // The total set of constraints, i.e. those due to hanging nodes and due to the Dirichlet boundary conditions
    // have to be distributed to the sparsity of the potential.
    constraints.condense (sparsity_patterns(PatternInfo::FEMwDirBC));

    // For the ions we only need the constraints which take care of the hanging nodes.
    hanging_node_constraints.condense (sparsity_patterns(PatternInfo::FEMwvNBC) );

    // After all constraints have been applied we can compress the sparsity patterns ...
    sparsity_patterns.compress();

    // ... and reinitialize the matrices.
    matrices.reinit();


    mg_Phi.reinit_Dir_bc(mg_dof_handler, Phi_dirichlet_boundary);

    // For the ion densities there no Dirichlet boundary conditions.
    // At least not those of the potential.
    mg_Ions.reinit_vN_bc(mg_dof_handler);


    // At the end, print some info about the dof distribution.
    std::cout << "   Number of degrees of freedom: "
              << mg_dof_handler.n_dofs()
              << " (by level: ";
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
        std::cout << mg_dof_handler.n_dofs(level)
                  << (level == triangulation.n_levels()-1 ? ")" : ", ");
    std::cout << std::endl;
}



// @sect4{DRSProblem::initialize_solution}
// Resizes the solution vector to the right size. Overwrites any previous content.
template <int dim>
void DRSProblem<dim>::initialize_solution()
{
    this->new_solution.reinit(mg_dof_handler.n_dofs());
    this->vN_unit_test_new_solution.reinit(mg_dof_handler.n_dofs());
}

// TODO: inherit the following from step27::FEMBEMPoissonProblem
template <int dim>
void DRSProblem<dim>::LinearSystem::setup_block_system(const SparseMatrixSet28 &matrices,
                                                       const step27::BlockPatternInfo& block_pattern,
                                                       const SolutionSet & sol_set
                                                       )
{
    // Before setting up the PDE system we copy the mapping of the names of the solution components to their indices ...
    this->component_names2index = sol_set.component_names2index;

    A.reinit(n_blocks(), n_blocks());
    A_nl.reinit(n_blocks(), n_blocks());
    A_2dt.reinit(n_blocks(), n_blocks());
    block_prec.reinit(n_blocks());

    typename step27::BlockPatternInfo::const_iterator ij = block_pattern.begin(),
            end_ij = block_pattern.end();

    for (; ij != end_ij; ++ij) {
        step27::IJ key = ij->first;

        uint row = sol_set.block_id(key.first);
        uint col = sol_set.block_id(key.second);

       typename step27::BlockPatternInfo::mapped_type::const_iterator
        m = ij->second.begin(),
                end_m = ij->second.end();
        for (; m != end_m; ++m)
            if (m->m_tag != MtxInfo::SLFullMatrix)// would be more appropriate: MtxInfo::SPInfo::Dense)
            {
                A.enter(matrices(m->m_tag), row, col, m->mtx_coeff);
                if (m->pc_tag != MtxInfo::None)
                   qDebug() << "pc name : " <<  m->pc_tag;
                   // block_prec.enter(precs(m->pc_tag), row, col, m->pc_coeff);
            }
            else {
                A.enter(matrices.fm(m->m_tag), row, col, m->mtx_coeff);
                if (m->pc_tag != MtxInfo::None)
                   qDebug() << "pc name : " <<  m->pc_tag;
                // for a future version: block_prec.enter(precs(m->pc_tag), row, col, m->pc_coeff);
            }
    }
}



// For the FEM-BEM unit test we need a 2x2 block matrix of the form
//
// \f{equation}
//
//   {\bf A} = \left(
//    \begin{array}{cc}
//         ( \psi \nabla, \nabla \phi )_{\Omega} &   -( \psi, t_P )_{\Gamma \subset \partial \Omega}  \\
//         2\pi I - K  & V
//     \end{array}
// \f}
// At the beginning we multiplied the BEM part of the system by $\pi$ in order to save some divisions and to make the
// shorter. Therefore, at this point all matrices are taken as they are, except for the upper off-diagonal block
// for which we have to reverse the sign.
// To stress which block keeps its sign and which one not we explicitly give all coefficients including their signs.





template <int dim>
void DRSProblem<dim>::LinearSystem::setup_pc (const step27::BlockPatternInfo & pc_block_pattern,
                                              const step27::GeoMGPC<dim> & mg_pc_Phi,
                                              const step27::GeoMGPC<dim> & mg_pc_Ions,
                                              const /*step27::*/SparseMatrixSet28 & matrices)
{

    this->block_prec.clear();

    typename step27::BlockPatternInfo::const_iterator ij = pc_block_pattern.begin(),
            end_ij = pc_block_pattern.end();

    for (; ij != end_ij; ++ij) {
        step27::IJ key = ij->first;

        uint row = this->block_id(key.first);
        uint col = this->block_id(key.second);

        typename step27::BlockPatternInfo::mapped_type::const_iterator
                m = ij->second.begin(),
                end_m = ij->second.end();
        for (; m != end_m; ++m)
            switch(m->m_tag) {
            case MtxInfo::LaplaceAllvNPC:
                block_prec.enter(mg_pc_Ions.preconditioner, row, col, m->mtx_coeff);
                break;
            case MtxInfo::LaplaceDRSPhiPC:
                block_prec.enter(mg_pc_Phi.preconditioner, row, col, m->mtx_coeff);
                break;
            case MtxInfo::Identity:
                block_prec.enter(bem_preconditioner, row, col, m->mtx_coeff);
                break;
            case MtxInfo::BEMFEMLowerOffBlock:
              block_prec.enter(matrices(m->m_tag), row, col, m->mtx_coeff);
            break;
            default:
                break;
            }
    }
}


template <int dim>
void DRSProblem<dim>::LinearSystem::solve()
{
    dealii::SolverControl solver_control (100000, 1e-6);

    // The problem is not symmetric anymore. Thus we use GMRES
    typedef dealii::SolverGMRES<dealii::BlockVector<double> > SolverType;
    SolverType::AdditionalData a_data(50 /*size of ortho basis*/, false /*right pc?*/ );
    SolverType  solver (solver_control); //, a_data);

    //===================
    //
    //       SOLVE
    //
    //===================
    // std::cout << __FUNCTION__ << " start solving" << std::endl;

    solver.solve (A, sol, rhs, block_prec);

    std::cout << "   " << solver_control.last_step()
              << " GMRES iterations needed to obtain convergence. "
              << std::endl;
}




#ifdef USE_NEWTON
template <int dim>
void DRSProblem<dim>::block_solve_Newton()
{
    std::cout << __FUNCTION__ << " begin" << std::endl;

    // The last step is the application of the constraints of the FEM problem
    std::vector<bool> all_constraints(SolutionSet::n_components, false);
    std::vector<bool> hnc_constraints(SolutionSet::n_components, true);


    // The potential needs Dirichlet data and not only hnc.
    all_constraints[SolutionSet::Phi_id] = true;
    hnc_constraints[SolutionSet::Phi_id] = false;


    if (LinearSystem::cats_idx >= 0)
        lsys.old_sol.block(LinearSystem::cats_idx)   = this->old_solution(step27::Cations);

    if (LinearSystem::neut_idx >= 0)
        lsys.old_sol.block(LinearSystem::neut_idx)   = this->old_solution(step27::Neutral);

    if (LinearSystem::ans_idx >= 0)
        lsys.old_sol.block(LinearSystem::ans_idx)    = this->old_solution(step27::Anions);

    if (LinearSystem::fem_pot_idx >= 0) {
        lsys.old_sol.block(LinearSystem::fem_pot_idx) = this->old_solution(step27::Phi);
    }


    // well ...
    if (LinearSystem::bem_pot_idx >= 0)
        for (unsigned int K = 0; K < bem_form->n_bc_dofs; ++K)
        {
            unsigned int k_fem_index = fem_bem.dof_to_boundary_mapping.size() -  bem_form->n_bc_dofs + K;

            lsys.old_sol.block(LinearSystem::bem_pot_idx)(K) = this->new_solution(step27::Phi_derivative)(k_fem_index);
        }


    // Newton rhs
    {
        lsys.rhs = 0.;
        lsys.A_nl.vmult(lsys.rhs, lsys.old_sol);
        // Newton-specific
        lsys.rhs *= -1.;

        lsys.A_2dt.vmult_add(lsys.rhs, lsys.sol);


        if (LinearSystem::fem_pot_idx >= 0)
            lsys.rhs.block(LinearSystem::fem_pot_idx) += this->system_rhs;


        if (LinearSystem::bem_pot_idx >= 0)
            for (unsigned int K = 0; K < bem_form->n_bc_dofs; ++K)
            {
                unsigned int k_fem_index = fem_bem.dof_to_boundary_mapping.size() -  bem_form->n_bc_dofs + K;
                lsys.rhs.block(LinearSystem::bem_pot_idx)(K) += fem_bem.newton_potential_rhs(k_fem_index);
            }
    }

#ifdef USE_OLD_MGPC
    // The MG preconditioners need a fully initialized dof handler.
    GeoMGPC<dim>  mg_pc_Phi (mg_Phi, mg_dof_handler,
                             hanging_node_constraints,
                             this->mg_params.n_smoothing_steps);

    GeoMGPC<dim> mg_pc_Ions (mg_Ions, mg_dof_handler,
                             hanging_node_constraints,
                             this->mg_params.n_smoothing_steps);

    // Unlike the matrices the preconditioners have to be set over and over again.
    lsys.setup_pc(mg_pc_Phi, mg_pc_Ions,
                  this->matrices);
#endif
    // A few steps before we go to the next mesh.
    for (uint nl_iter = 0; nl_iter < 2; nl_iter++)
    {
        lsys.solve();
        double l2_cats_incr = 0;
        double l2_neut_incr = 0;
        double l2_ans_incr = 0;
        double l2_pot_incr = 0;

        // After solving copy the solution components back
        if (LinearSystem::cats_idx >= 0)
        {
            l2_cats_incr =  lsys.sol.block(LinearSystem::cats_idx).l2_norm();
        }

        if (LinearSystem::neut_idx >= 0)
        {
            l2_neut_incr = lsys.sol.block(LinearSystem::neut_idx).l2_norm();
        }

        if (LinearSystem::ans_idx >= 0)
        {
            l2_ans_incr = lsys.sol.block(LinearSystem::ans_idx).l2_norm();
        }

        if (LinearSystem::fem_pot_idx >= 0)
        {
            l2_pot_incr =  lsys.sol.block(LinearSystem::fem_pot_idx).l2_norm();
        }

        lsys.old_sol += lsys.sol;



        std::cout << " l2 incr cast : " << l2_cats_incr << ", neut : " << l2_neut_incr
                  << ", ans : " << l2_ans_incr << ", pot : " << l2_pot_incr << std::endl;

        reassemble_nonlinearity_cpu();

        // Newton rhs
        {
            lsys.rhs = 0.;
            lsys.A_nl.vmult(lsys.rhs, lsys.old_sol);
            // Newton-specific
            lsys.rhs *= -1.;

            lsys.A_2dt.vmult_add(lsys.rhs, lsys.sol);


            if (LinearSystem::fem_pot_idx >= 0)
                lsys.rhs.block(LinearSystem::fem_pot_idx) += this->system_rhs;


            if (LinearSystem::bem_pot_idx >= 0)
                for (unsigned int K = 0; K < bem_form->n_bc_dofs; ++K)
                {
                    unsigned int k_fem_index = fem_bem.dof_to_boundary_mapping.size() -  bem_form->n_bc_dofs + K;
                    lsys.rhs.block(LinearSystem::bem_pot_idx)(K) += fem_bem.newton_potential_rhs(k_fem_index);
                }
        }


    }
    lsys.solve();

    // final update
    lsys.sol += lsys.old_sol;

    lsys.block_prec.clear();


    // After solving copy the solution components back
    if (LinearSystem::cats_idx >= 0)
        this->new_solution(step27::Cations) = lsys.sol.block(LinearSystem::cats_idx);

    if (LinearSystem::neut_idx >= 0)
        this->new_solution(step27::Neutral) = lsys.sol.block(LinearSystem::neut_idx);

    if (LinearSystem::ans_idx >= 0)
        this->new_solution(step27::Anions) = lsys.sol.block(LinearSystem::ans_idx);

    if (LinearSystem::fem_pot_idx >= 0)
        this->new_solution(step27::Phi) = lsys.sol.block(LinearSystem::fem_pot_idx);

    if (LinearSystem::bem_pot_idx >= 0) {
        fem_bem.solution_num = lsys.sol.block(LinearSystem::bem_pot_idx);

        for (unsigned int K = 0; K < bem_form->n_bc_dofs; ++K)
        {
            new_solution(step27::Phi_derivative)(fem_bem.dof_to_boundary_mapping.size()
                                         -
                                         bem_form->n_bc_dofs + K) = lsys.sol.block(LinearSystem::bem_pot_idx)(K);
        }
    }

    new_solution.distribute_constraints(constraints, all_constraints);
    new_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);


    std::cout << "l2 norm cats : " << this->new_solution(step27::Cations).l2_norm()
              << ", neuts : " << this->new_solution(step27::Neutral).l2_norm()
              << ", ans : " << this->new_solution(step27::Anions).l2_norm()
              << std::endl;

    std::cout << __FUNCTION__ << " end" << std::endl;
}
#endif



template <int dim>
void DRSProblem<dim>::block_solve_fix_point(const uint cycle)
{
    std::cout << __FUNCTION__ << " begin" << std::endl;

    // TODO: hnc constraints flags are part of the PDE setup -> setup_block_system
    // The last step is the application of the constraints of the FEM problem
    std::vector<bool> all_constraints(new_solution.n_blocks(), false);
    std::vector<bool> hnc_constraints(new_solution.n_blocks(), true);

    // The potential needs Dirichlet data and not only hnc.
    all_constraints[new_solution.block_id(step27::Phi)] = true;
    hnc_constraints[new_solution.block_id(step27::Phi)] = false;

    double inv_pseudo_dt = 1./this->numerics_prms.pseudo_time_step;


   /*new_*/old_solution.distribute_constraints(constraints, all_constraints);
    /*new_*/old_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);

    // The new solution contains the most recent solution after its interpolation ono the current mesh.
    // The old_solutionwill be used for measuring the increment within an iteration.
    // Therefore, we have to copy it.
    //old_solution = new_solution;

    // We need a few subsets of our solution components for the different parts of the final assembly.
    const std::vector<SolComps> & ions = drs_sol_comps_fem_part_all_vN;
    const uint n_ions = ions.size();

    // Use the last known solution as initial value for the Krylov solver ...
    lsys.set_solution(this->/*new_*/old_solution);

    // ... and set the right hand side for the potential.
    // The other right-hand side in the FEM part are zero due to the DRS model.
    lsys.rhs.block(this->new_solution.block_id(step27::Phi)) = this->system_rhs;


    // Add the contribution due to the pseudo-time stepping.
    for (uint c = 0; c < n_ions; c++)
    {
        uint comp = this->old_solution.block_id(ions[c]);
        this->matrices(MtxInfo::BulkMass).vmult(  lsys.rhs.block(comp),
                                                  this->old_solution(ions[c]) );
        lsys.rhs.block(comp) *= inv_pseudo_dt;
    }

    // well ...
    if (lsys.component_names2index.find(step27::Phi_derivative) != lsys.component_names2index.end() )
        // fem_bem.newton_potential_rhs.extract_subvector_to(index_range.begin(), index_range.end(), lsys.rhs.block(bem_pot_idx).begin());
        for (unsigned int K = 0; K < bem_form->n_bc_dofs; ++K)
        {
            unsigned int k_fem_index = bem_form->dof_to_boundary_mapping.size() -  bem_form->n_bc_dofs + K;

            lsys.rhs.block(lsys.block_id(step27::Phi_derivative)/*LinearSystem::bem_pot_idx*/)(K) = fem_bem.newton_potential_rhs(k_fem_index);
            lsys.sol.block(lsys.block_id(step27::Phi_derivative)/*LinearSystem::bem_pot_idx*/)(K) = this->/*new_*/old_solution(step27::Phi_derivative)(k_fem_index);
        }


    // The MG preconditioners need a fully initialized dof handler and a thousand other things.
    // Because of their complex design it is best to re-instantiate a pristine object whenever it is needed.
    step27::GeoMGPC<dim>  mg_pc_Phi (mg_Phi, mg_dof_handler,
                                     hanging_node_constraints,
                                     this->numerics_prms.n_smoothing_steps);

    step27::GeoMGPC<dim> mg_pc_Ions (mg_Ions, mg_dof_handler,
                                     hanging_node_constraints,
                                     this->numerics_prms.n_smoothing_steps);

    // Unlike the matrices the preconditioners have to be set over and over again.
    lsys.setup_pc(drs_pc_block_pattern, mg_pc_Phi, mg_pc_Ions,
                  this->matrices);

    lsys.sol.collect_sizes(); lsys.rhs.collect_sizes();
    // A few steps before we go to the next mesh.
    for (int nl_iter = 0; nl_iter < this->numerics_prms.n_time_steps; nl_iter++)
    {
        lsys.solve();

        std::map<SolComps, double> l2_incr;

        // After solving copy the solution components back and compute the l2 norm of the increments.
        lsys.compute_incr_and_push_sol_to(l2_incr, old_solution);

        // Make new solution conform to the constraints due to hanging nodes and boundary conditions.
        old_solution.distribute_constraints(constraints, all_constraints);
        old_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);

        // Set new linearization point for the linear system. Basically, this the current contetn of lsys.sol but
        // modified such that it fulfills all constrainst due to BCs or HNCs.
        lsys.set_solution(old_solution);

        new_solution = old_solution;
        if (this->sub_model == step27::fem_bem_unit_test)
            output_results(cycle, nl_iter);

        compute_electrode_fluxes();

        qDebug() << "time step : " << this->numerics_prms.pseudo_time_step << ", l2 incr cats : " << l2_incr[step27::Cations] << ", neut : " << l2_incr[step27::Neutral]
                 << ", ans : " << l2_incr[step27::Anions] << ", pot : " << l2_incr[step27::Phi] << ", Dphi : " << l2_incr[step27::Phi_derivative];

        if (!this->drs_sol_comps_fem_part_all_vN.empty())
            reassemble_nonlinearity_cpu();

        // Add the contribution due to the pseudo-time stepping.
        for (uint c = 0; c < n_ions; c++)
        {
            uint comp = this->old_solution.block_id(ions[c]);

            this->matrices(MtxInfo::BulkMass).vmult(lsys.rhs.block(comp), this->old_solution(ions[c]) );
            lsys.rhs.block(comp) *= inv_pseudo_dt;
        }
    }

    if (this->numerics_prms.n_time_steps > 0)
        lsys.solve();


    lsys.block_prec.clear();


    // After solving copy the solution components back
    lsys.push_sol_to(new_solution);

    // Make new solution conform to the constraints due to hanging nodes and boundary conditions.
    new_solution.distribute_constraints(constraints, all_constraints);
    new_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);
  if (this->sub_model == step27::fem_bem_unit_test)
    output_results(cycle, this->numerics_prms.n_time_steps);

    if (lsys.component_names2index.find(step27::Phi_derivative) != lsys.component_names2index.end() )
        fem_bem.solution_num = lsys.sol.block(lsys.component_names2index[step27::Phi_derivative]);

    if (!drs_sol_comps_fem_part_all_vN.empty())
        std::cout << "l2 norm cats : " << this->new_solution(step27::Cations).l2_norm()
                  << ", neuts : " << this->new_solution(step27::Neutral).l2_norm()
                  << ", ans : " << this->new_solution(step27::Anions).l2_norm() << std::endl;

    std::cout << __FUNCTION__ << " end" << std::endl;
}



// @sect4{Postprocessing}

// The following two functions postprocess a
// solution once it is computed. In
// particular, the first one refines the mesh
// at the beginning of each cycle while the
// second one outputs results at the end of
// each such cycle. The functions are almost
// unchanged from those in step-6, with the
// exception of two minor differences: The
// KellyErrorEstimator::estimate function
// wants an argument of type DoFHandler, not
// MGDoFHandler, and so we have to cast from
// derived to base class; and we generate
// output in VTK format, to use the more
// modern visualization programs available
// today compared to those that were
// available when step-6 was written.
template <int dim_>
void DRSProblem<dim_>::refine_grid_global ()
{
    // flag some cells for refinement, e.g.
    // dealii::GridRefinement::refine_and_coarsen_fixed_fraction(*tria, error_indicators, 1., 0);
    // prepare the triangulation
    // for refinement,
    this->triangulation.set_all_refine_flags(); //
    this->triangulation.prepare_coarsening_and_refinement();


    // transfer solutions to new grid
    dealii::SolutionTransfer<dim, dealii::Vector<double> > soltrans(this->mg_dof_handler);

    SolutionSet & sol_on_old_grid = this->new_solution;



    // W.r.t. time the solution after interpolation is the "old" solution
    // which is needed in the tight-hand side of the time-steppign scheme.
    SolutionSet & sol_on_new_grid = this->old_solution;

    // tell the SolutionTransfer object
    // that we intend to do coarsening and refinement,
    soltrans.prepare_for_coarsening_and_refinement(sol_on_old_grid.data());

    // soltrans.prepare_for_coarsening_and_refinement(lsys.sol);

    // actually execute the refinement,
    this->triangulation.execute_coarsening_and_refinement();
    // and redistribute dofs.
    mg_dof_handler.distribute_dofs (fe);


    // resize solution vector to the correct
    // size, as the @p refine_interpolate
    // function requires the vectors to be
    // of right sizes
   // std::vector<typename SolutionSet::size_type> drs_block_sizes(n_drs_components, this->mg_dof_handler.n_dofs());
    sol_on_new_grid.reinit(mg_dof_handler.n_dofs());

    // lsys.old_sol.reinit();

    // and finally interpolate
    soltrans.interpolate(sol_on_old_grid.data(), sol_on_new_grid.data());

    // The last thing we do is to reinitialize the vector for the new solution we are going to compute.
    this->new_solution.reinit(mg_dof_handler.n_dofs()); // drs_block_sizes);
    vN_unit_test_new_solution.reinit(mg_dof_handler.n_dofs());
}

template <int dim_>
void DRSProblem<dim_>::refine_grid_adaptive ()
{
    dealii::Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    dealii::KellyErrorEstimator<dim>::estimate (static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                dealii::QGauss<dim-1>(3),
                                                typename dealii::FunctionMap<dim>::type(),
                                                new_solution(step27::Phi),
                                                estimated_error_per_cell);
    dealii::GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                             estimated_error_per_cell,
                                                             0.3, 0.03);
    // triangulation.execute_coarsening_and_refinement ();

    // Only refine cells at cavity surface
    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();


    // ------------------------- CELL LOOP BEGIN -----------------------
    //    unsigned int n_cells_to_refine = 0;
    //    for (; cell != endc; ++cell)
    //    {
    //        if(cell->at_boundary() )
    //            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
    //                if(cell->face(f)->boundary_indicator() == DRS::BCids::cavity() )
    //                {
    //                    cell->set_refine_flag();
    //                    n_cells_to_refine++;
    //                }
    //    }

    //std::cout << "n cells to refine : " << n_cells_to_refine << std::endl;

    this->triangulation.prepare_coarsening_and_refinement();


    // transfer solutions to new grid
    dealii::SolutionTransfer<dim, dealii::Vector<double> > soltrans(this->mg_dof_handler);

    SolutionSet & sol_on_old_grid = this->new_solution;

    // W.r.t. time the solution after interpolation is the "old" solution
    // which is needed in the tight-hand side of the time-steppign scheme.
    SolutionSet & sol_on_new_grid = this->old_solution;

    // tell the SolutionTransfer object
    // that we intend to do coarsening and refinement,
    soltrans.prepare_for_coarsening_and_refinement(sol_on_old_grid);

    // actually execute the refinement,
    this->triangulation.execute_coarsening_and_refinement();
    // and redistribute dofs.
    mg_dof_handler.distribute_dofs (fe);


    // resize solution vector to the correct
    // size, as the @p refine_interpolate
    // function requires the vectors to be
    // of right sizes
    sol_on_new_grid.reinit(mg_dof_handler.n_dofs());

    // and finally interpolate
    soltrans.interpolate(sol_on_old_grid, sol_on_new_grid);

    // after interpolation we have to resize the future new solution as well.
    this->new_solution.reinit(mg_dof_handler.n_dofs());
    this->vN_unit_test_new_solution.reinit(mg_dof_handler.n_dofs());
}



template <int dim_>
void DRSProblem<dim_>::output_results (const unsigned int cycle, const int iteration) const
{

    dealii::DataOut<dim> data_out;

    const dealii::DoFHandler<dim> & dh = mg_dof_handler;

    data_out.attach_dof_handler (dh/*mg_dof_handler*/);

   // for (uint c = 0; c < drs_sol_comps_fem_part.size(); c++)


    data_out.add_data_vector (this->new_solution(step27::Phi), "potential");
    //
    data_out.add_data_vector (fem_bem.solution_ref, "pot_ref_sol");

    if (drs_sol_comps_fem_part_all_vN.size()>0)
    {
        data_out.add_data_vector (this->new_solution(step27::Neutral), "neutrals");
        data_out.add_data_vector (this->new_solution(step27::Anions), "anions");
        data_out.add_data_vector (this->new_solution(step27::Cations), "cations");
    }

    data_out.add_data_vector (this->new_solution(step27::Phi_derivative), "n_x_grad_Phi_BEM_part");

    data_out.build_patches(mapping, this->degree);

    std::ostringstream filename;
    filename << "solution-"
             << fem_bem.testcase_name.c_str()
             << "-q" << numerics_prms.fe_degree
             << "-m" << numerics_prms.fe_mapping_degree
             << "-p_bem" << numerics_prms.bem_quad_order_1
             << "-" << cycle;
    if (iteration > -1)
        filename << "-" << iteration;

    filename << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);


    // error of FEM solution
    if (true)
    {
        dealii::Vector<float> difference_per_cell (this->triangulation.n_active_cells());
        double L2_error;
        double H1_error;

        if (this->sub_model == step27::fem_bem_unit_test)
        {
        // potential
        dealii::VectorTools::integrate_difference (mapping,
                                                   mg_dof_handler,
                                                   this->new_solution(step27::Phi),
                                                   *fem_bem.u_ref,
                                                   difference_per_cell,
                                                   dealii::QGauss<dim>(2*fe.degree+1),
                                                   dealii::VectorTools::L2_norm);
        double L2_error = difference_per_cell.l2_norm();

        dealii::VectorTools::integrate_difference (mapping,
                                                   mg_dof_handler,
                                                   this->new_solution(step27::Phi),
                                                   *fem_bem.u_ref,
                                                   difference_per_cell,
                                                   dealii::QGauss<dim>(2*fe.degree+1),
                                                   dealii::VectorTools::H1_norm);
        double H1_error = difference_per_cell.l2_norm();

        std::cout << "FEM, dealii::VectorTools::integrate_difference : final Pot ||u - u_h ||_L2 : " << L2_error
                  << " ||u - u_h ||_H1 : " << H1_error << std::endl;
}

        // von Neumann problem
        difference_per_cell = 0.;

    //    std::cout << "FEM, dealii::VectorTools::integrate_difference : final von Neumann ||u - u_h ||_L2 : " << L2_error_vN
      //            << "||u - u_h ||_H1 : " << H1_error_vN << std::endl;

        this->unit_test_table.add_value("L2 error pot FEM", L2_error );
        this->unit_test_table.set_scientific("L2 error pot FEM", true);
        this->unit_test_table.set_precision("L2 error pot FEM", 4);

        this->unit_test_table.add_value("H1 error pot FEM", H1_error );
        this->unit_test_table.set_scientific("H1 error pot FEM", true);
        this->unit_test_table.set_precision("H1 error pot FEM", 4);


        this->unit_test_table.add_value("L2 error pot BEM", fem_bem.L2_diff_u_h );
        this->unit_test_table.set_scientific("L2 error pot BEM", true);
        this->unit_test_table.set_precision("L2 error pot BEM", 4);

        this->unit_test_table.add_value("L2 error dnu BEM", fem_bem.L2_diff_dnu_h );
        this->unit_test_table.set_scientific("L2 error dnu BEM", true);
        this->unit_test_table.set_precision("L2 error dnu BEM", 4);

    //    std::cout << "BEM: ||u - u_h ||_L2_surf : " << fem_bem.L2_diff_u_h
      //            << ", ||t - t_h ||_L2_surf : " << fem_bem.L2_diff_dnu_h << std::endl;
    }


    // Compute the fluxes through electrodes.
    this->compute_electrode_fluxes();
}


// @sect4{DRSProblem::run}

// Like several of the functions above, this
// is almost exactly a copy of of the
// corresponding function in step-6. The only
// difference is the call to
// <code>assemble_multigrid</code> that takes
// care of forming the matrices on every
// level that we need in the multigrid
// method.
template <int dim_>
void DRSProblem<dim_>::run ()
{

    // initialize
    {
        std::cout << "Cycle " << 0 << ':' << std::endl;

        // The setup of the solution vectors can only be done after the setup of the system as we need the fem bem data.
        std::vector<typename SolutionSet::size_type> drs_block_sizes(n_drs_components, this->mg_dof_handler.n_dofs());

        drs_block_sizes[  this->old_solution.block_id(step27::Phi_derivative) ] = bem_form->n_bc_dofs;

        this->old_solution.reinit(this->mg_dof_handler.n_dofs()); // drs_block_sizes);
        this->new_solution.reinit(this->mg_dof_handler.n_dofs()); // drs_block_sizes);

       //  std::vector<typename SolutionSet::size_type> vN_block_sizes(1, this->mg_dof_handler.n_dofs() );
        this->vN_unit_test_new_solution.reinit(this->mg_dof_handler.n_dofs()); // vN_block_sizes);
        this->vN_unit_test_old_solution.reinit(this->mg_dof_handler.n_dofs()); // vN_block_sizes);


        // template <int dim>
        // void DRSProblem<dim>::set_initial_conditions()
        if (this->sub_model == step27::full_DRS)
        {
            // At the beginning ions are homogeneously distributed.

            dealii::ConstantFunction<dim> cations_ic(this->phys_data.cations_av_density, 1);
            dealii::ConstantFunction<dim> anions_ic(this->phys_data.anions_av_density, 1);
            dealii::ConstantFunction<dim> neutrals_ic(this->phys_data.neutrals_av_density, 1);

            const dealii::DoFHandler<dim> & dh = this->mg_dof_handler;

            // Initial conditions for the ionic species.
            {
                dealii::VectorTools::interpolate(this->mapping, dh,
                                                 cations_ic,
                                                 this->old_solution(step27::Cations));

                dealii::VectorTools::interpolate(this->mapping, dh,
                                                 neutrals_ic,
                                                 this->old_solution(step27::Neutral));

                dealii::VectorTools::interpolate(this->mapping, dh,
                                                 anions_ic,
                                                 this->old_solution(step27::Anions));
            }

            // Apply inhomogeneous Dirichlet conditions to initial condition of potential.
            // Later we will only add increments which fulfill homogeneous Dirichlet conditions.
            {
                // The potential has non-zero values at the electrodes
                // but otherwise it is zero.
                dealii::ConstantFunction<dim>           cathode_dirichlet_bc(
                            + phys_data.potential_difference, 1);

                dealii::ConstantFunction<dim>           anode_dirichlet_bc(
                            - phys_data.potential_difference, 1);

                dealii::ConstraintMatrix dirichlet_bc;


                typename dealii::FunctionMap<dim>::type       Phi_dirichlet_boundary;


                Phi_dirichlet_boundary[DRS::BCids::cathode()] = &cathode_dirichlet_bc;
                Phi_dirichlet_boundary[DRS::BCids::anode()] = &anode_dirichlet_bc;

                dealii::VectorTools::interpolate_boundary_values (mapping,
                                                                  static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                                  Phi_dirichlet_boundary,
                                                                  dirichlet_bc);
                dirichlet_bc.close();
                dirichlet_bc.distribute(this->old_solution(step27::Phi) );
            }

        }


        //check for matrix identity:
        assemble_system_generic (0/*cycle*/, this->arch);

        assemble_multigrid ();

        block_solve_fix_point (0);

        output_results (0/*cycle*/);
    }

#ifndef sdkgfalsdfa
    for (unsigned int cycle=1; cycle<this->numerics_prms.n_mesh_refinements; ++cycle)
    {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        refine_grid_global ();

        setup_discretization ();



        lsys.reinit(fem_bem.n_fem_dofs, bem_form->n_bc_dofs);



        dealii::Timer timer;  timer.restart();

        //  if(numerics_prms.n_threads > 1) {
        //      assemble_system_multithreaded(cycle);
        //    } else
        {
            assemble_system_generic (cycle, this->arch);
        }
        std::cerr << "assemble_system took " << timer.wall_time() << " seconds."
                  << std::endl;

        //compute_bem_assembly_errors();

        assemble_multigrid ();

        block_solve_fix_point (cycle);

        output_results (cycle);


    }
#endif

    // Write convergence table
#ifndef nUSE_FEM

    std::ostringstream filename_base;
    std::string folder = "conv_results";
    mkdir(folder.c_str(), 0755);

    filename_base << fem_bem.testcase_name.c_str()
             << "fem_"
             <<   "bem_conv_" // "harmonic"
               << "_q" << numerics_prms.fe_degree
               << "_m" << numerics_prms.fe_mapping_degree
               << "_p_bem" << numerics_prms.bem_quad_order_1
               << ".txt" << std::ends;

    std::ofstream ut_out (( folder + "/unit_test_" + filename_base.str()).c_str());

    ut_out << "# " <<           ", FE degree : " << numerics_prms.fe_degree
        << ", MAPPING order : " << numerics_prms.fe_mapping_degree
        << ", inner BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_1
        << ", outer BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_2
        << std::endl;

    this->unit_test_table.write_text(ut_out);

    ut_out  << std::endl;
    ut_out  << std::endl;

     std::ofstream currents_out (( folder + "/currents_" + filename_base.str()).c_str());

     this->measured_currents_table.write_text(currents_out);



    static const bool eval_red_rate = false;

    if (eval_red_rate) {
        fem_bem.convergence_table
                .evaluate_convergence_rates("L2(step27::Phi)", dealii::ConvergenceTable::reduction_rate_log2);
        //        fem_bem.convergence_table
        //                .evaluate_convergence_rates("L_infty(step27::Phi)", ConvergenceTable::reduction_rate_log2);

        fem_bem.convergence_table
                .evaluate_convergence_rates("L2(DnPhi)", dealii::ConvergenceTable::reduction_rate_log2);
        //        fem_bem.convergence_table
        //                .evaluate_convergence_rates("L_infty(DnPhi)", ConvergenceTable::reduction_rate_log2);

        fem_bem.convergence_table
                .evaluate_convergence_rates("L_infty(angles)", dealii::ConvergenceTable::reduction_rate_log2);

        fem_bem.convergence_table
                .evaluate_convergence_rates("L2(Phi-FEM)", dealii::ConvergenceTable::reduction_rate_log2);
    }


    dealii::deallog << std::endl;

    std::cout << "init ref : " << numerics_prms.n_init_refinements
              << ", FE degree : " << numerics_prms.fe_degree
              << ", MAPPING order : " << numerics_prms.fe_mapping_degree
              << ", inner BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_1
              << ", outer BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_2
                 // <<  ", SG_QUAD_ORDER : " << SQ_ORDER
              << std::endl;

    fem_bem.convergence_table.write_text(std::cout);

     std::ofstream conv_out ((folder + "/fem_bem_conv_" + filename_base.str()).c_str());

   // fem_bem.convergence_table.write_text(conv_out);
#else
    std::cout << "init ref : " << numerics_prms.n_init_refinements
              << ", FE degree : " << numerics_prms.fe_degree
              << ", MAPPING order : " << numerics_prms.fe_mapping_degree
              << ", BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order
                 // <<  ", SG_QUAD_ORDER : " << SQ_ORDER
              << std::endl;
#endif // USE_FEM


    std::ofstream timing_out( numerics_prms.timing_file
                              , std::ios_base::app | std::ios_base::out);

    timing_table.write_text(timing_out);

    timing_out.flush();
    timing_out.close();

}



// From a biophysical point-of-view the followign function is the most interesting part
// of this program
template <int dim_>
void DRSProblem<dim_>::compute_electrode_fluxes () const
{
    // For the integration over the surface of the electrodes
    // the quadrature rule has to be accurate enough to integrate the
    // polynomials used for the finite elements. We do not have to integrate
    // over products of polynomials as in the assembla of the bilinear forms.
    // Therefore, using a Gauss rule of order FE degree plus one should suffice.
    const dealii::QGauss<dim -1> q_rule(this->fe.degree+2);

    // We need only one flux form in order to evaluate
    // $K_O\int_{\Gamma_A}c_0d\Gamma_A$ and $K_R\int_{\Gamma_C}c_+d\Gamma_C$.
    FluxForm flux_at_electrode(this->mapping, this->fe, q_rule);

    double flux_A = 0.,flux_C = 0.;


    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();


     const unsigned char anode_id = DRS::BCids::anode();

     const unsigned char cathode_id = DRS::BCids::cathode();

      SolComps anode_comp = step27::Neutral;
      SolComps cathode_comp = step27::Cations;

      if (this->sub_model == step27::fem_bem_unit_test)
          anode_comp = cathode_comp = step27::Phi;

     // Loop over the cells and compute the surface integral
     // if a cell is at the right boundary.
    for (; cell!=endc; ++cell)
    {
        if (cell->at_boundary() )
            for (unsigned int f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
            {
                if(cell->face(f)->boundary_indicator() == anode_id )
                    // At the anode we have to integrate over the densities of the neutral particles.
                    flux_A += flux_at_electrode (cell, f, new_solution(anode_comp) );

                if(cell->face(f)->boundary_indicator() == cathode_id )
                    // At the cathode we have to integrate over the densities of the cations.
                    flux_C += flux_at_electrode (cell, f, new_solution(cathode_comp) );
            }
    } // cell loop END

    double rel_error = 2*(flux_A - flux_C)/(flux_A + flux_C);

    std::cout << "FLUXES : I_C : " << flux_C << ", I_A : " << flux_A << ", relative error : " << rel_error << std::endl;

    // Once the cell loop is done, we write everything into a convergence table.
    // At the current stage of modeling the redox processes
    // only enter via a global scaling of the computed fluxes.
    // This can be done outside of the program.
    this->measured_currents_table.add_value("I_C", flux_C);
    this->measured_currents_table.set_scientific("I_C", false);
    this->measured_currents_table.set_precision("I_C", 12);

    this->measured_currents_table.add_value("k_R", this->phys_data.k_R);
    this->measured_currents_table.set_scientific("k_R", true);
    this->measured_currents_table.set_precision("k_R", 6);


    this->measured_currents_table.add_value("I_A", flux_A);
    this->measured_currents_table.set_scientific("I_A", false);
    this->measured_currents_table.set_precision("I_A", 12);

    this->measured_currents_table.add_value("k_O", this->phys_data.k_O);
    this->measured_currents_table.set_scientific("k_O", true);
    this->measured_currents_table.set_precision("k_O", 6);
}



// @sect4{Class: Step16::DRSProblem}




} // END namespace Step16


