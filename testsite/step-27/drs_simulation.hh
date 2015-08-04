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
// comment on them:
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

#include<deal.II/base/parameter_handler.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>


#include <deal.II/grid/grid_out.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/block_vector.h>


// This is C++:
#include <fstream>
#include <sstream>
#include <iomanip>

// A bit of Qt for creating directories
#include <QDir>
#include <QDebug>

// Because of the size of the project, the source code is
// distributed over several files. They roughly represent
// the different problems one has to solve during development.

// Physical aspects comprise the computational domain
#include <step-27/drs_cell.hh>
// the type of FEM-BEM coupling and the definition of
// the corresponding weak forms for the boundary integrals
// For choosing the runtime architecture
#include <step-27/Architecture.h>
#include <step-27/coupling_type.h>


#include <step-27/BoundaryDoFTools.hh>

#include <step-27/UnitTestProblemData.hh>

#include <step-27/bemforms.h>
#include <step-27/bemforms.hh>

#include <step-27/LaplaceForm.hh>

// The preconditioner from deal.II's step-16 is modified and put into a separate file.
#include <step-27/geomgpc.h>

// Taylor-made linear algebra classes
#include <step-27/drs_lac.hh>

// The run-time parameters are declared in separate files.
// From a programming point of view they are not particularly interesting
// and therefore do not appear in the online documentation.
// Therefore, they are declared in separate files.
#include <step-27/SimParams.h>
#include <step-27/NumericsParams.h>
#include <step-27/PhysicalParams.h>

#define USE_ADAP
#undef USE_ADAP

#define USE_FEM
// #undef USE_FEM


#define USE_FEM_BEM
//#undef USE_FEM_BEM

// #define USE_OLD_STUFF

// The last step is as in all
// previous programs:
namespace step27
{

// Forward declarations.
template <int dim> struct LaplaceKernel;
template <int dim> struct MGData;
template <int dim> struct GeoMGPC;
template <int dim> class FEMBEMPoissonProblem;


// @sect3{struct FemBemData}
// TODO: Separate into a base class and derived classes for collocation and weak form.
// Additional base class for the profiling stuff (convergence, runtimes).
// TO DO: spin off error estimation stuff.
template <int dim>
struct FemBemData {

    typename dealii::FunctionMap<dim>::type interface;


    const dealii::Function<dim> * u_ref;

    const dealii::Function<dim> * newton_potential;


    // Vector containing the values of the FEM-BEM reference solution
    // on the whole domain, not only on the boundary.
    dealii::Vector<double>       solution_ref;
    dealii::Vector<double>       solution_num;

    dealii::Vector<double>       solid_angles;

    dealii::FullMatrix<double> DL_matrix;
    dealii::FullMatrix<double> alpha_plus_DL_matrix;

    // needed in future extensions:
#ifdef USE_WEAK_FORM
     dealii::SparsityPattern hypersingular_to_laplace_sparsity;
    dealii::FullMatrix<double> hypersingular_matrix;
    dealii::SparseMatrix<double> hypersingular_to_laplace;

    dealii::FullMatrix<double> SL_matrix_weak;
    dealii::FullMatrix<double> DL_matrix_weak;
    dealii::FullMatrix<double> DL_plus_mass_matrix_weak;


    this->hypersingular_matrix.reinit(n_bc_dofs, n_bc_dofs); // XXX
    //this->hypersingular_to_laplace.reinit(n_fem_dofs, n_fem_dofs); // XXX
    this->SL_matrix_weak.reinit(n_bc_dofs, n_bc_dofs); // XXX
    this->DL_matrix_weak.reinit(n_bc_dofs, n_bc_dofs); // XXX
    this->DL_plus_mass_matrix_weak.reinit(n_bc_dofs, n_bc_dofs); // XXX

    this->hypersingular_matrix = 0.;
    //this->hypersingular_to_laplace = 0.;
    this->SL_matrix_weak  = 0.;
    this->DL_matrix_weak  = 0.;
    this->DL_plus_mass_matrix_weak  = 0.;

#endif
#ifdef USE_MG_PC_4_BEM
    // XXX END
    dealii::MGLevelObject<FullMatrixAccessor<double> > SL_mg_matrices; // TO DO: make part of preconditioner
#endif


    dealii::Vector<double>       newton_potential_rhs;

    // Only useful for unit testing
    dealii::Vector<double>       rhs_for_dirichlet_test;

    dealii::ConvergenceTable	convergence_table;

    dealii::Vector<double>  u_h, dnu_h;


    double L2_diff_u_h;
    double L2_diff_dnu_h;

    double L8_diff_u_h;
    double L8_diff_dnu_h;

    std::string testcase_name;


    unsigned int n_fem_dofs;
    unsigned int n_bulk_dofs;
    unsigned int n_bc_dofs;

    FemBemData( const FemBemData &other){
        //            this = other;
    }

    FemBemData & operator=( FemBemData &other){
        //            *this = other;
        return *this;
    }

    FemBemData (const step27::TestExpr test_epxr_flag, const std::string & prm_path, const QDir &log_dir)
    {

        std::cout << __FUNCTION__ << std::endl;

        switch (test_epxr_flag) {
        case step27::sum:
            u_ref =  new const step27::UnitTests::HarmonicTest<dim, step27::UnitTests::Sum<dim> >(prm_path, log_dir);
            newton_potential = new const typename step27::UnitTests::Sum<dim>::NewtonPotential();
            testcase_name = "sum";
            break;
        case step27::prod:
            u_ref =  new const step27::UnitTests::HarmonicTest<dim, step27::UnitTests::Prod<dim> >(prm_path, log_dir);
            newton_potential = new const typename step27::UnitTests::Prod<dim>::NewtonPotential();
            testcase_name = "prod";
            break;
        case step27::dipole:
            u_ref =  new const step27::UnitTests::HarmonicTest<dim, step27::UnitTests::Dipole<dim> >(prm_path, log_dir);
            newton_potential = new const typename step27::UnitTests::Dipole<dim>::NewtonPotential(prm_path, log_dir);
            testcase_name = "dipole";
            break;
        case step27::sum_prod:
            u_ref =  new const step27::UnitTests::HarmonicTest<dim, step27::UnitTests::SumProd<dim> >(prm_path, log_dir);
            newton_potential = new const typename step27::UnitTests::SumProd<dim>::NewtonPotential();
            testcase_name = "sum_prod";
            break;
        case step27::sum_prod_dipole:
            u_ref =  new const step27::UnitTests::HarmonicTest<dim, step27::UnitTests::SumProdDipole<dim> >(prm_path, log_dir);
            newton_potential = new const typename step27::UnitTests::SumProdDipole<dim>::NewtonPotential(prm_path, log_dir);
            std::cout << "sum_prod_dipole";
            break;
        default:
            break;
        }

        std::cout << testcase_name.c_str() << " testcase chosen" << std::endl;
    }


    void reinit(const unsigned char bc_id,
                const dealii::DoFHandler<dim>& fem_dof_handler)
    {

        this->interface[bc_id] = 0;

        this->solution_ref.reinit(fem_dof_handler.n_dofs() );
        this->u_h.reinit(fem_dof_handler.n_dofs() );
        this->dnu_h.reinit(fem_dof_handler.n_dofs() );
        this->newton_potential_rhs.reinit(fem_dof_handler.n_dofs());

        this->n_fem_dofs  = fem_dof_handler.n_dofs();
        this->n_bc_dofs   = fem_dof_handler.n_boundary_dofs(interface);
        this->n_bulk_dofs = this->n_fem_dofs - this->n_bc_dofs;

        this->solid_angles.reinit(n_bc_dofs);

        this->DL_matrix.reinit(n_bc_dofs, n_bc_dofs);
        this->alpha_plus_DL_matrix.reinit(n_bc_dofs, n_bc_dofs);


        this->solution_num.reinit(n_bc_dofs);


        this->rhs_for_dirichlet_test.reinit(n_bc_dofs);

        this->solid_angles = 0.;

        this->DL_matrix = 0.;
        this->alpha_plus_DL_matrix = 0.;


        this->newton_potential_rhs = 0.;
    }

    void copy_fem2bem(dealii::Vector<double> & dst, const dealii::Vector<double> & src) const
    {
      //  src.extract_subvector_to(index_range.begin(), index_range.end(), dst.begin());
    }

    void copy_bem2fem(dealii::Vector<double> & dst, const dealii::Vector<double> & src) const
    {
        for (unsigned int K = 0; K < n_bc_dofs; ++K)
        {
        //    dst(dof_to_boundary_mapping.size()
          //                                -
            //                              n_bc_dofs + K) = src(K);
        }

    }
};

// This structure collects identifiers for the different types of sparsity patterns we are going to need.
// Basically, we have to distinguish three cases:
//
// - a pattern belongs to a matrix from the FE part of the PDE system. Then its name only contains "FEM".
// - a pattern represents a boundary element matrix. Then it only has n_bem_dofs rows and columns.
// - the third case are the patterns for the off-diagonal blocks in the PDE system coupling the FEM DoFs to the BEM DoFs.
// In contrast to the others these matrices are not square.
//
// Within the FEM matrices we have to distinguish whether a sparsity pattern has to take care of Dirichlet boundary values or not.
// Most of the patterns are for the different mass matrices on the different subboundaries.
struct PatternInfo {

    enum Ids { FEMwDirBC=111, FEMwvNBC=112, FEMBulkMass=113,
               FEMCathodeBoundaryMass=114, FEMAnodeBoundaryMass=115, // tag for the result of scattering make_boundary_sparsity_pattern to the full FEM system
               BEMCathodeBoundaryMass=116, BEMAnodeBoundaryMass=117, // tag for the result of make_boundary_sparsity_pattern
               FEMBEMUpperOffBlock=118, BEMFEMLowerOffBlock=119 /* patterns for coupling BEM to FEM*/,
               FEMBEM=120, // id for the sparsity pattern of the BEM part of the FEM-BEM problem.
               Dense=130 // id for dense matrices to signal that no sparsity pattern should be stored.
             };

    typedef Ids PatternIds;
};

struct MtxInfo {
    // For future extensions we add also the tags for the drift terms.
    enum Ids { None = -1,
               Identity=1,
               LaplacePhi=10,
               LaplaceIons=20,
               BulkMass=30,
               CathodeBoundaryMass=40,
               AnodeBoundaryMass=50,
               DriftCations=60,
               DriftAnions=70,
               DriftPot=80,
               FEMBEMUpperOffBlock=90,
               BEMFEMLowerOffBlock=100,
               FEMBEMMass=110,
               SLFullMatrix=120,
               LaplaceAllvNPC=130,
               LaplaceDRSPhiPC=140,

             };

    MtxInfo();

    // Placeholder for making PatternInfo a nested class.
    // At the end, it is the matrices which decide whoch pattern they need.
    typedef PatternInfo SPInfo;

    typedef std::map<Ids, PatternInfo::Ids> Id2Pattern;

    // Sinc we know our equation structure in advance, we can set up a simple array which maps the matrix id to its type of sparsity pattern (finite element with Dirichlet boundary conditions, only with von Neumann conditions, mass matrix on a boundary, i.e. a surface mass matrix)
    Id2Pattern pattern_types; // = { PatternInfo::FEMwDirBC, PatternInfo::FEMwvNBC, PatternInfo::FEMwvNBC, PatternInfo::CathodeBoundaryMass, PatternInfo::AnodeBoundaryMass, PatternInfo::FEMwvNBC, PatternInfo::FEMwvNBC, PatternInfo::FEMwvNBC };
};

MtxInfo::MtxInfo() {

    // We need two Laplace matrices because the potential
    // has Dirichlet BC values,
    // whereas the different ion species satisfy some
    // sort of Neumann BC on all boundaries.
    // The only difference in the patterns are the boundary conditions.
    // Thus, we can use a pattern several times, e.g. for the Laplacian for the ion species and the mass matrix.
    pattern_types[LaplacePhi] = PatternInfo::FEMwDirBC;
    pattern_types[LaplaceIons] = PatternInfo::FEMwvNBC;
    // For pseudo-time stepping and the ionic contributions
    // to the Poisson equation for the potential we need a mass matrix.
    pattern_types[BulkMass] = PatternInfo::FEMwvNBC;
    // TO DO: correct pattern assignment
    //  pattern_types[CathodeBoundaryMass] = PatternInfo::FEMCathodeBoundaryMass; // PatternInfo::CathodeBoundaryMass;
    //  pattern_types[AnodeBoundaryMass] = PatternInfo::FEMCathodeBoundaryMass; // PatternInfo::AnodeBoundaryMass;
    //  pattern_types[DriftCations] = PatternInfo::FEMwvNBC;
    //  pattern_types[DriftAnions] = PatternInfo::FEMwvNBC;
    //  pattern_types[DriftPot] = PatternInfo::FEMwvNBC;

    pattern_types[FEMBEMUpperOffBlock] = PatternInfo::FEMBEMUpperOffBlock;
    pattern_types[BEMFEMLowerOffBlock] = PatternInfo::BEMFEMLowerOffBlock;
    pattern_types[FEMBEMMass] = PatternInfo::FEMBEM;

    pattern_types[SLFullMatrix] = PatternInfo::Dense;
}




// We need quite few solution vectors. To simplify administrative things
// like resizing we create a little structure which does all these things
// and which holds a few references so that we can address the individual
// solution components by name. We implement it based on a container
// from the standard library.

// TO DO: redesign!!!
// TO DO: check whether BlockVEctor can be put into SolutionTransfer, if yes, derive from BlockVector, Avoids copying in block_solve()

// We need quite few solution vectors. To simplify administrative things
// like resizing we create a little structure which does all these things
// and which holds a few references so that we can address the individual
// solution components by name. We implement it based on a container
// from the standard library.
struct SolutionSet : public dealii::BlockVector< double >
{
    typedef dealii::BlockVector<double> Base;

    typedef Base::size_type size_type;

    // Surprisingly, dealii::SolutionTransfer cannot work on dealii::BlockVectors.
    // To fix this, we have to provide access to the protected data of dealii::BlockVectorBase.
    std::vector< dealii::Vector<double> > & data() { return this->components; }

    const std::vector< dealii::Vector<double> > & data() const { return this->components; }

    std::map<SolComps, uint> component_names2index;


    const unsigned int block_id(SolComps s) const {
        typename std::map<SolComps, uint>::const_iterator m = component_names2index.find(s);

        if (!(m != component_names2index.end()))
        {
            qDebug() << "n comps : " << component_names2index.size();

             std::map<SolComps, uint>::const_iterator i = component_names2index.begin(),
                     endc = component_names2index.end();
             for (; i != endc; ++i)
             {
                 std::cerr << "name : " << int(i->first) << ", idx : " << i->second << std::endl;
             }

        }

        Assert(m != component_names2index.end(), dealii::ExcMessage("Illegal component id") );

        unsigned int c = m->second;
        return c;
    }

    // Lookup of individual components by name.
    const dealii::Vector<double> & operator () (const SolComps s) const { return this->block(block_id(s)); }


    dealii::Vector<double> & operator () (const SolComps s) { return this->block(block_id(s)); }


    SolutionSet(const uint n_dofs_per_comp=1)
        :
          Base(1, n_dofs_per_comp)
    {}


    // For some unknown reason functions from the base class do not get resolved on OSX 10.9.
    void reinit(const std::vector<size_type>& block_sizes)
    {
           this->Base::reinit(block_sizes);
    }


    void reinit(const uint new_block_size) { this->Base::reinit(this->n_blocks(), new_block_size); }

    // It is assumed that the order of @p block_sizes is the same as in @p sc and that both have the same size.

    void reinit(const std::vector<SolComps>& sc, const std::vector<size_type>& block_sizes)
    {
      const uint n_components = /*block_sizes*/sc.size();

        this->Base::reinit(block_sizes);

        component_names2index.clear();
        for (int k = 0; k < n_components; k++)
        {
            component_names2index[sc[k]] = k;
           qDebug() << "size of comp " << k << ": " << this->block(k).size();
        }

        this->collect_sizes();
    }


// Distribute the given constraints on the selected components, i.e. blocks, of the solution block vector.
    void distribute_constraints(const dealii::ConstraintMatrix& constraints,
                                std::vector<bool> selected_comps=std::vector<bool>())
    {
         const uint n_components = this->n_blocks();

        if (selected_comps.empty())
        {
            for (int c = 0; c < n_components; c++)
                constraints.distribute( this->block(c) );
        }
        else {

            for (int c = 0; c < n_components; c++)
                if (selected_comps[c])
                    constraints.distribute( this->block(c) );

        }
    }
};



// @sect3{The <code>FEMBEMPoissonProblem</code> class template}

// This main class is basically the same
// class as in step-6. As far as member
// functions is concerned, the only addition
// is the <code>assemble_multigrid</code>
// function that assembles the matrices that
// correspond to the discrete operators on
// intermediate levels:
template <int dim>
class FEMBEMPoissonProblem
{
public:
    FEMBEMPoissonProblem (const step27::NumericsParams<dim> &numerics_prms,
                          const std::string & prm_path,
                          const QDir &log_dir,
                          step27::TestExpr testcase,
                          step27::SubModel /*sm*/);

    ~FEMBEMPoissonProblem()
    {
        if (bem_form)
            delete bem_form;
        bem_form = 0;
    }

    void run ();

    void block_solve_fix_point();


protected:
    void setup_system ();
    void initialize_solution();

protected:
    //void assemble_system (uint cycle);

    void assemble_system_generic (uint cycle, step27::Architecture);

    void assemble_boundary_mass_matrices();


    void assemble_multigrid ();

    void assemble_bem_unit_test();

    void setup_preconditioner();


    void compute_bem_assembly_errors();

    void solve_bem_unit_test();


    void compute_L2_bem_error();

    void refine_grid_global ();
    void refine_grid_adaptive ();

    void output_results (const unsigned int cycle) const;


    //Function to generate matrix values, to be run in a thread
    void generate_bem_matrix_values(const unsigned int n_bem_points, const unsigned int n_bem_quad_points,
                                    const unsigned int dofs_per_be,
                                    const std::vector<dealii::Point<dim> > x, const std::vector<dealii::Point<dim> > q,
                                    const dealii::FullMatrix<double> W,
                                    const std::vector<unsigned int> global_bc_dofs,
                                    const std::vector<double> normals,
                                    const std::vector<dealii::Tensor<1, dim> > u_gradients, const std::vector<double> JxW);
    // Copy std::vector<Point> to std::vector<double>


    const step27::NumericsParams<dim> & numerics_prms;

    dealii::ConvergenceTable     timing_table;

    dealii::Triangulation<dim>   triangulation;

    dealii::MappingQ<dim>        mapping;

    dealii::FE_Q<dim>            fe; // for DG this is the space Q

    // FE_DGQ<dim> fe_dgq;

    dealii::MGDoFHandler<dim>    mg_dof_handler;


    // For the finite element part of the problem we need
    // an object, which assembles the Laplace operator.
    LaplaceForm<dim> laplace;

    // The boundary elements for the interface
    // need a different quadrature rule which only integrates over surfaces and not over volumes.
    // The matrix entries are assembled by the BEMForm. In this step it is based on collocation.
    BemCollocationForm<dim, FaceIntTraits<dim> > * bem_form;

   // BemFormBase<dim, FaceIntTraits> * bem_form;

    // Similar to the @p SolutionSet we pool all sparsity patterns
    // of the FEM part of the DRS problem in a structure in order
    // to centralize and simplify the memory management.

    // For simplified access we define a list of Ids for the different structures which we will encounter.

    typedef SparsityPatternSet<PatternInfo> PatternSet;

    PatternSet sparsity_patterns;

    typedef SparseMatrixSet<MtxInfo> MatrixSet;

    MatrixSet matrices;

    // The @p blcik_pattern will tell us how to populate the PDE system.
    step27::BlockPatternInfo block_pattern;

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



     const std::vector<SolComps> vN_sol_comps;
     //[] = { Neutral };
  //  static const int n_vN_unit_components = 1;

     const std::vector<SolComps> unit_test_fem_bem_sol_comps; //[] = { Phi, Phi_derivative };
    //static const int n_unit_test_fem_bem_components = 2;



    SolutionSet old_solution;
    SolutionSet new_solution;


    dealii::Vector<double>       system_rhs;


    struct LinearSystem {

        // Before we set up the blocked system of algebraic equations we define
        // some named constants for the block positions of the individual components.
        // In case of a negative index the corresponding component is omitted.
        // This is primarily for testing and development.

        // mapping of physical quantities to rows of blocks


        static const int fem_pot_idx = 0;
        static const int bem_pot_idx = 1;

        unsigned int n_blocks;

        dealii::BlockMatrixArray<double> A;


        dealii::PreconditionIdentity bem_preconditioner;


        LinearSystem() : n_blocks(0) {

            // Automatic adjustment of the block size.
            // This is needed during development.
            if (fem_pot_idx >= 0) n_blocks++;
            if (bem_pot_idx >= 0) n_blocks++;

            A.reinit(n_blocks, n_blocks);

            block_prec.reinit(n_blocks);
        }

        void clear() {
            std::cout << __FUNCTION__ << std::endl;
            A.clear();

            block_prec.clear();
        }

        ~LinearSystem() {

            this->clear();

            std::cout << __FUNCTION__ << " DONE" << std::endl;
        }

        std::vector<unsigned int> block_sizes;


        dealii::BlockVector<double> rhs;


        dealii::BlockVector<double> sol;
        dealii::BlockVector<double> old_sol;


        dealii::BlockTrianglePrecondition<double> block_prec;

        dealii::BlockVector<double> rhs_fem_bem_unit_test;
        dealii::BlockVector<double> sol_fem_bem_unit_test;


        void setup_block_system(const MatrixSet & matrices, const FemBemData<dim> &fem_bem,
                                const PhysicalParams & phys_prm,
                                const BlockPatternInfo& block_pattern,
                                const SolutionSet & sol_set);

        void setup_pc (const GeoMGPC<dim> & mg_pc_Phi,
                       // const GeoMGPC<dim> & mg_pc_Ions,
                       const MatrixSet & matrices);

        void reinit (const unsigned int n_fem_dofs,
                     const unsigned int n_bem_dofs)
        {
            {
                std::vector<unsigned int> tmp;
                std::swap(tmp, block_sizes);
            }


            block_sizes.resize(n_blocks, n_fem_dofs);

            if (bem_pot_idx >= 0)
                block_sizes[bem_pot_idx] = n_bem_dofs;


            rhs.reinit(block_sizes);
            sol.reinit(block_sizes);
            old_sol.reinit(block_sizes);
        }


        void solve();
    };
    // Put fem_bem data before lsys, because lsys connects its internal matrices to the fem_bem SL matrix,
    // so fem_bem must be destroyed after lsys! Otherwise, dealii catches fire.
    FemBemData<dim> fem_bem;

    LinearSystem lsys;


    const unsigned int degree;

    // Although ther different problems can share their sparsity patterns
    // we need separate matrices because the potential has Dirichlet boundary
    // conditions while the ions have not.
    // Therefore, we need two independent data sets and eventually preconditioners.
    MGData<dim> mg_Phi;
    MGData<dim> mg_Ions;

    PhysicalParams phys_data;


    mutable dealii::ConvergenceTable unit_test_table;
};


// @sect3{The <code>FEMBEMPoissonProblem</code> class implementation}

// @sect4{FEMBEMPoissonProblem::FEMBEMPoissonProblem}

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
template <int dim>
FEMBEMPoissonProblem<dim>::FEMBEMPoissonProblem (const step27::NumericsParams<dim> &numerics_prms,
                                                 const std::string &prm_path,
                                                 const QDir &log_dir,
                                                 step27::TestExpr testcase,
                                                 step27::SubModel /*sm*/)
    :
      numerics_prms(numerics_prms),
      triangulation (dealii::Triangulation<dim>::
                     limit_level_difference_at_vertices),
      mapping (numerics_prms.fe_mapping_degree),
      fe (numerics_prms.fe_degree),
      mg_dof_handler (triangulation),
      old_solution(),
      new_solution(),
      degree (numerics_prms.fe_degree),
      fem_bem(testcase, prm_path),
      laplace(mapping, fe, *numerics_prms.fem_q_rule),
      matrices(sparsity_patterns)

    //    ,
    //      bem_form (DRS::BCids::cavity(), // TO DO: pass id of an dielectric interface, not of the particular one
    //                mapping,
    //                fe,
    //                inner_surface_quad)
{
    // All our reference solutions are computed for $\varepsilon = 1$.
    if (testcase != step27::drs)
    {
        this->phys_data.EPSILON_S = 1.;
        this->phys_data.EPSILON_P = 1.;
    }

    bem_form = new BemCollocationForm<dim, FaceIntTraits<dim> >  (DRS::BCids::cavity(), // TO DO: pass id of an dielectric interface, not of the particular one
                                                            mapping,
                                                            fe,
                                                            numerics_prms);


    const SolComps sol_comps [] = { step27::Phi, step27::Phi_derivative, step27::Neutral };


    const_cast<std::vector<SolComps>& >(unit_test_fem_bem_sol_comps).assign(sol_comps, sol_comps + 2);



    std::vector<uint> block_sizes(unit_test_fem_bem_sol_comps.size(), 1); // Dummy sizes.
    old_solution.reinit(unit_test_fem_bem_sol_comps, block_sizes);

    new_solution.reinit(unit_test_fem_bem_sol_comps, block_sizes); //unit_test_fem_bem_sol_comps, block_sizes);

#ifndef nUSE_NEWSETUP
    // ---------------------------- setup_fem_bem_unit_test_block_pattern()

    const double EPS_S = phys_data.EPSILON_S;
    const double EPS_P = phys_data.EPSILON_P;

    block_pattern[IJ(step27::Phi, step27::Phi)].push_back({ MtxInfo::LaplacePhi, 1. });
    block_pattern[IJ(step27::Phi_derivative, step27::Phi_derivative)].push_back({ MtxInfo::SLFullMatrix, +1.*(EPS_S/EPS_P) });

    block_pattern[IJ(step27::Phi, step27::Phi_derivative)].push_back({MtxInfo::FEMBEMUpperOffBlock, -1.*(EPS_P/EPS_S)});
    block_pattern[IJ(step27::Phi_derivative, step27::Phi)].push_back({MtxInfo::BEMFEMLowerOffBlock, +1. });


    // -------------------------- END -------------------
#endif

    // Last but not least, create the coarse grid.
    // At the beginning we set up the coarse mesh.
    DRS::MeshGenerator::DRSCellParams drs_cell_params;

    // Before we set up the triangulation we have to read the parameters
    // which determine its size.
    dealii::ParameterHandler prm_handler;

    DRS::MeshGenerator::DRSCellParams::declare(prm_handler);
    std::cout << "DRS prms taken from : " << numerics_prms.prm_drs_cell_filepath.c_str() << std::endl;
    prm_handler.read_input(numerics_prms.prm_drs_cell_filepath);

    drs_cell_params.get(prm_handler);


    std::ofstream log_drs_cell_prm( (log_dir.absolutePath() + QDir::separator() + QFileInfo(numerics_prms.prm_drs_cell_filepath.c_str()).fileName() + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_drs_cell_prm,
                                  dealii::ParameterHandler::Text);

    DRS::MeshGenerator::drs_cell(triangulation, drs_cell_params, //numerics_prms.prm_drs_cell_filepath,
                                 numerics_prms.n_init_refinements);

}

// @sect4{FEMBEMPoissonProblem::setup_system}

// The following function extends what the
// corresponding one in deal.II's step-16 did, which extended a function from deal.II's step-6.
// For all reactive boundaries, i.e. the electrodes and the dielectric interface,
// this function sorts the indices of the boundary dofs to the back such that the respective indices form contiguous subsets.
// This is necessary for setting up the PDE system as a block of scalar equations.
template <int dim>
void FEMBEMPoissonProblem<dim>::setup_system ()
{
    // As remainder from deal.II's step-16 we keep the output of the
    // degrees of freedom in the multilevel structure.
    dealii::deallog << "Number of degrees of freedom: "
                    << mg_dof_handler.n_dofs();

    dealii::DoFHandler<dim> & dh = static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler);

    // Before we deal with the dof enumeration of the FEM-BEM coupling
    // we cope with the dofs on the electrodes.

    std::set< dealii::types::boundary_id > cathode_boundary_indicator;
    std::vector<bool> cathode_dof_flags;

    std::set< dealii::types::boundary_id > anode_boundary_indicator;
    std::vector<bool> anode_dof_flags;

    cathode_boundary_indicator.insert(DRS::BCids::cathode());
    anode_boundary_indicator.insert(DRS::BCids::anode());

    // In order to extract the DoFs on the selected boundary,
    //  dealii::DoFTools::extract_boundary_dofs needs to know
    // for which component of a multi-component FE problem.
    // Since we cobstruct our multi-component PDE problem from assigning matrices assembled
    // for a scalar, i.e. single-component, problem to different entries in a BlockMatrixArray,
    // Ww have only one component. So the following is trivial but necessary.
    std::vector< bool > component_select(1 /*n_components*/);
    component_select[0] = true;

    cathode_dof_flags.resize(mg_dof_handler.n_dofs());
    dealii::DoFTools::extract_boundary_dofs	(dh,
                                             component_select,
                                             cathode_dof_flags,
                                             cathode_boundary_indicator);
    // Now sort back. At the end ALL boundary dofs will be selected back but their relative order
    // is given by the order in which we tackle the different subboundaries.
    dealii::DoFRenumbering::sort_selected_dofs_back (static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                     cathode_dof_flags);

    anode_dof_flags.resize(mg_dof_handler.n_dofs());
    dealii::DoFTools::extract_boundary_dofs	(static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                             component_select,
                                             anode_dof_flags,
                                             anode_boundary_indicator);

    dealii::DoFRenumbering::sort_selected_dofs_back (static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                     anode_dof_flags);

    // For the FEM-BEM coupling we sort the dofs of the potential
    // on the surface of the cavity to the back.
    // This leads to a dense block in the lower right of the stiffness matrix.

    // Only after the final dof numbering scheme is established
    // we can setup the BEM problem
    bem_form->reinit_bem_data(mg_dof_handler, this->mapping);

    fem_bem.reinit(DRS::BCids::cavity(), mg_dof_handler);
    // After the FEM do the sparsity patterns
    // for the off-diagonal coupling matrices
    {
       // std::cout << __FILE__ << " " << __LINE__ << std::endl;

        // Using operator[] automatically adds an entry to the underlying std::map of the collection
        // of sparsity patterns, if it does not exist, yet.
        this->sparsity_patterns[PatternInfo::FEMBEM].reinit(fem_bem.n_bc_dofs,
                                        fem_bem.n_bc_dofs,
                                        mg_dof_handler.max_couplings_between_dofs());

       //  std::cout << __FILE__ << " " << __LINE__ << std::endl;

        dealii::DoFTools::make_boundary_sparsity_pattern(
                    static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                    bem_form->interface,
                    bem_form->dof_to_boundary_mapping,
                    this->sparsity_patterns(PatternInfo::FEMBEM) );

        // std::cout << __FILE__ << " " << __LINE__ << std::endl;

        this->sparsity_patterns(PatternInfo::FEMBEM).compress();

        //  std::cout << __FILE__ << " " << __LINE__ << std::endl;

        sparsity_patterns.reinit_FEM_BEM(mg_dof_handler, fem_bem);

      //   std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

    //   std::cout << __FILE__ << " " << __LINE__ << std::endl;

    // Only after all the dof sorting is done we can setup the sparsity patterns.
    sparsity_patterns.reinit_FEM (mg_dof_handler);

    for (unsigned int l=0;l<triangulation.n_levels();++l)
        dealii::deallog << "   " << 'L' << l << ": "
                        << mg_dof_handler.n_dofs(l);
    dealii::deallog  << std::endl;

    system_rhs.reinit (mg_dof_handler.n_dofs());

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

    typename dealii::FunctionMap<dim>::type & Phi_dirichlet_boundary_unit_test
            = Phi_dirichlet_boundary;

    Phi_dirichlet_boundary_unit_test[DRS::BCids::anode()]   = fem_bem.u_ref;
    Phi_dirichlet_boundary_unit_test[DRS::BCids::cathode()] = fem_bem.u_ref;

    Phi_dirichlet_boundary_unit_test[DRS::BCids::cylinder_surface()] = fem_bem.u_ref;

    dealii:: VectorTools::interpolate_boundary_values (mapping,
                                                       static_cast<const dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                       Phi_dirichlet_boundary_unit_test,
                                                       constraints/*_fem_bem_unit_test*/);

    constraints.close ();

    hanging_node_constraints.close ();

    hanging_node_constraints.distribute (fem_bem.solution_ref);

    constraints.condense (sparsity_patterns(PatternInfo::FEMwDirBC) ); // .Phi);

    hanging_node_constraints.condense (sparsity_patterns(PatternInfo::FEMwvNBC) ); //.Ions);

    // After all constraints have been applied we can compress the sparsity patterns.
    sparsity_patterns.compress();

    matrices.reinit();

    // The MG preconditioner needs to know about Dirichlet boundary conditions.
    // In later steps, when we also simulate the ion densities which do not have Dirichlet boundary conditions
    // (at least not those of the potential) the differences will become more clear.
    mg_Phi.reinit_Dir_bc(mg_dof_handler, Phi_dirichlet_boundary);

    // Although we do not need the ion-specific MG preconditioner, we keep it for later steps
    // such that we do not have to re-write the assemble_multigrid function.
    mg_Ions.reinit_vN_bc(mg_dof_handler);

    // At the end, print some info about the dof distribution.
    std::cout << "   Number of degrees of freedom: " << mg_dof_handler.n_dofs() << " (by level: ";
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
        std::cout << mg_dof_handler.n_dofs(level)
                  << (level == triangulation.n_levels()-1 ? ")" : ", ");
    std::cout << std::endl;
}

// @sect4{FEMBEMPoissonProblem::initialize_solution}
// Resizes the solution vector to the right size. Overwrites any previous content.
template <int dim>
void FEMBEMPoissonProblem<dim>::initialize_solution()
{
    this->new_solution.reinit(mg_dof_handler.n_dofs());
}

// The error of the BEM part is computed from the definition of a surface L2 norm
// as surface integrals over the squared distance between numerical and reference solution.
template <int dim>
void FEMBEMPoissonProblem<dim>::compute_L2_bem_error ()
{
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    // -------- BEM DATA BEGIN ---------------



    dealii::FEFaceValues<dim> bem_values (mapping,
                                          fe,
                                          *numerics_prms.bem_q_rule_1,
                                          dealii::update_values|dealii::update_normal_vectors|dealii::update_quadrature_points|dealii::update_JxW_values);

    const unsigned char bc_id       = *bem_form->boundary_indicators.begin();

    // Return the normal vectors at the quadrature points.
    // For a face, these are the outward normal vectors to the cell.
    // For a cell of codimension one, the orientation is
    // given by the numbering of vertices.
    const std::vector<dealii::Point<dim> > & normals = bem_values.get_normal_vectors ();

    const std::vector<dealii::Point<dim> > & q  = bem_values.get_quadrature_points ();

    const unsigned int n_bem_quad_points= bem_form->n_q_points;

    const std::vector<double> & JxW     = bem_values.get_JxW_values();


    // For testing the integrals

    std::vector<double>          u_values    (n_bem_quad_points);
    std::vector<double>          u_h_values  (n_bem_quad_points);
    std::vector<double>          dnu_h_values  (n_bem_quad_points);

    std::vector<dealii::Tensor<1, dim> > u_gradients (n_bem_quad_points);


    // -------- BEM DATA END ---------------

    // copy solutions
    fem_bem.u_h = new_solution(Phi); // unit_test_solution_Pot_fem_bem();
    fem_bem.dnu_h = new_solution(Phi_derivative); // unit_test_Phi_derivative();

    fem_bem.dnu_h *= -1;

    typename dealii::MGDoFHandler<dim>::active_cell_iterator
            cell = mg_dof_handler.begin_active(),
            endc = mg_dof_handler.end();

    fem_bem.L2_diff_u_h = 0.;
    fem_bem.L2_diff_dnu_h = 0.;

    fem_bem.L8_diff_u_h = 0.;
    fem_bem.L8_diff_dnu_h = 0.;
    // ------------------------- CELL LOOP BEGIN -----------------------
    for (; cell!=endc; ++cell)
    {
        // ------------- BEM stuff ------------
        if (cell->at_boundary() )
        {
            //  bool print_cell_dofs = true;
            cell->get_dof_indices (local_dof_indices);

            for (unsigned int f = 0; f < dealii::GeometryInfo<dim>::faces_per_cell; f++)
            {
                if(cell->face(f)->boundary_indicator() == bc_id )
                {
                    bem_values.reinit (cell, f);
                    fem_bem.u_ref->value_list (q, u_values);

                    bem_values.get_function_values(fem_bem.u_h, u_h_values);
                    bem_values.get_function_values(fem_bem.dnu_h, dnu_h_values);

                    fem_bem.u_ref->gradient_list(q, u_gradients);

                    for (unsigned int a=0; a < n_bem_quad_points; ++a)
                    {
                        //                        if(a == 0) {
                        //                            std::cout << "u_h_values[a]: " << u_h_values[a] << "\t u_values[a]: "<<  u_values[a] << "\t delta: " << u_h_values[a] - u_values[a] << std::endl;
                        //                        }
                        double diff_u_h_a   = std::fabs(u_h_values[a] - u_values[a]);
                        double diff_dnu_h_a = std::fabs(- dnu_h_values[a] - (normals[a] * u_gradients[a]));

                        fem_bem.L2_diff_u_h   += std::pow(diff_u_h_a, 2.) * JxW[a];

                        fem_bem.L2_diff_dnu_h += std::pow(diff_dnu_h_a, 2.) * JxW[a];

                        fem_bem.L8_diff_u_h = std::max(fem_bem.L8_diff_u_h, diff_u_h_a);
                        fem_bem.L8_diff_dnu_h = std::max(fem_bem.L8_diff_dnu_h, diff_dnu_h_a);
                    }
                }
            }
        }

    } // cell loop END

    fem_bem.L2_diff_u_h = std::sqrt( fem_bem.L2_diff_u_h );
    fem_bem.L2_diff_dnu_h = std::sqrt( fem_bem.L2_diff_dnu_h );
}

// @sect5{Function: setup_block_system}

template <int dim>
void FEMBEMPoissonProblem<dim>::LinearSystem::setup_block_system(const MatrixSet &matrices,
                                                                 const FemBemData<dim> & fem_bem,
                                                                 const PhysicalParams & phys_prm,
                                                                 const BlockPatternInfo& block_pattern,
                                                                 const SolutionSet & sol_set)
{
    typename BlockPatternInfo::const_iterator ij = block_pattern.begin(),
            end_ij = block_pattern.end();

    for (; ij != end_ij; ++ij) {
        IJ key = ij->first;

        uint row = sol_set.block_id(key.first);
        uint col = sol_set.block_id(key.second);

       typename BlockPatternInfo::mapped_type::const_iterator
        m = ij->second.begin(),
                end_m = ij->second.end();
        for (; m != end_m; ++m)
            if (m->m_tag != MtxInfo::SLFullMatrix)// would be more appropriate: MtxInfo::SPInfo::Dense)
                A.enter(matrices(m->m_tag), row, col, m->mtx_coeff);
        else
                A.enter(matrices.fm(m->m_tag), row, col, m->mtx_coeff);
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
//     \end{array}\right)
// \f}
// At the beginning we multiplied the BEM part of the system by $\pi$ in order to save some divisions and to make the
// shorter. Therefore, at this point all matrices are taken as they are, except for the upper off-diagonal block
// for which we have to reverse the sign.
// To stress which block keeps its sign and which one not we explicitly give all coefficients including their signs.


// TO DO: pass vector of PCs for later extension
template <int dim>
void FEMBEMPoissonProblem<dim>::LinearSystem::setup_pc (const GeoMGPC<dim> & mg_pc_Phi,
                                                        const MatrixSet & matrices)
{
    this->block_prec.clear();

    // First, do the diagonal entries.
    if (fem_pot_idx >= 0)
        block_prec.enter(mg_pc_Phi.preconditioner,
                         fem_pot_idx, fem_pot_idx);

    // If we do not add at least the identity matrix to the diagonal of the block preconditioner
    // the solution is going to be complete wrong because we would use a zero matrix as preconditioner for this block!
    if (bem_pot_idx >= 0)
        block_prec.enter( bem_preconditioner, bem_pot_idx, bem_pot_idx);


    // Now, do the off-diagonals
    if (fem_pot_idx >= 0 && bem_pot_idx >= 0)
        block_prec.enter(matrices(MtxInfo::BEMFEMLowerOffBlock), bem_pot_idx, fem_pot_idx, 1 );
}


template <int dim>
void FEMBEMPoissonProblem<dim>::LinearSystem::solve()
{
    dealii::SolverControl solver_control (100000, 1e-9);

    // The problem is not symmetric anymore. Thus we use GMRES
    typedef dealii::SolverGMRES<dealii::BlockVector<double> > SolverType;
    SolverType::AdditionalData a_data(50 /*size of ortho basis*/, false /*right pc?*/ );
    SolverType  solver (solver_control); //, a_data);

    //===================
    //
    //       SOLVE
    //
    //===================
    std::cout << __FUNCTION__ << " start solving" << std::endl;

    solver.solve (A, sol, rhs, block_prec);

    std::cout << "   " << solver_control.last_step()
              << " GMRES iterations needed to obtain convergence. "
              << std::endl;
}

template <int dim>
void FEMBEMPoissonProblem<dim>::block_solve_fix_point()
{
#ifdef TOGGLE_USE_OLD_VERSION
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
#else
    std::cout << __FUNCTION__ << " begin" << std::endl;

    // The last step is the application of the constraints of the FEM problem.
    // Depending on the solution component we need either all constraints,
    // i.e. the Dirichlet data on selected boundaries and the hanging node constraints,
    // or only the hanging node constraints, as it is the case for the BEM component of the solution
    // and solution components which only have von Neumann boundary data.
    // For initialization we assume that no component needs more than the hanging node constraints.
    std::vector<bool> all_constraints(new_solution.n_blocks(), false);
    std::vector<bool> hnc_constraints(new_solution.n_blocks(), true);

    // The potential needs Dirichlet data and not only hnc.
    all_constraints[new_solution.block_id(Phi)] = true;
   // all_constraints[SolutionSet::Phi_unit_id] = true;

    hnc_constraints[new_solution.block_id(Phi)] = false;
   //  hnc_constraints[SolutionSet::Phi_unit_id] = false;

    // Before we solve the linear system we distribute the constraints, so that
    // @p new_solution becomes a valid initial guess for the solution of the linear system.
    new_solution.distribute_constraints(constraints, all_constraints);
    new_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);

    // TO DO: cant we use the BlockVector right away?
    // The last task is to copy the vectors for the individual solution components into
    // a block vector so that we can use deal.II's linear algebra
    // which is built around the BlockMatrixArray class.
    if (LinearSystem::fem_pot_idx >= 0) {
        lsys.sol.block(LinearSystem::fem_pot_idx) = this->new_solution(Phi);
        lsys.rhs.block(LinearSystem::fem_pot_idx) = this->system_rhs;
    }

    // well ...
    if (LinearSystem::bem_pot_idx >= 0)
    {
        // TODO: make part of FemBemData or BemForm
        unsigned int k_fem_index = bem_form->dof_to_boundary_mapping.size() -  fem_bem.n_bc_dofs;
        std::vector<uint> index_range(fem_bem.n_bc_dofs);
        std::generate(index_range.begin(), index_range.end(), [&]{ return k_fem_index++; });

        // bem2fem operations (fem2b2m is analogous):
        fem_bem.newton_potential_rhs.extract_subvector_to(
                    index_range.begin(), index_range.end(),
                    lsys.rhs.block(LinearSystem::bem_pot_idx).begin());
        this->new_solution(Phi_derivative).extract_subvector_to(index_range.begin(), index_range.end(),
                                                                lsys.sol.block(LinearSystem::bem_pot_idx).begin());
    }

    // The MG preconditioners need a fully initialized dof handler.
    GeoMGPC<dim>  mg_pc_Phi (mg_Phi, mg_dof_handler,
                             hanging_node_constraints,
                             this->numerics_prms.n_smoothing_steps);

    // The preconditioners change as the mesh is refined.
    // Therefore, the preconditioners have to be set over and over again.
    lsys.setup_pc(mg_pc_Phi, this->matrices);

    // After all these preparations we can solve the linear algebraic system.
    // XXX    // A few steps before we go to the next mesh.
    // XXX    // for (int nl_iter = 0; nl_iter < this->numerics_prms.n_time_steps; nl_iter++)
    {
        lsys.solve();
        // TO DO: key-based access
        //        double l2_cats_incr = 0;
        //        double l2_neut_incr = 0;
        //        double l2_ans_incr = 0;
        double l2_pot_incr = 0;

        // After solving we have to copy the solution components back
        if (LinearSystem::fem_pot_idx >= 0)
        {
            // In order to monitor the relaxation of the solution during the pseudo-time stepping we compute the
            // euclidean norm of the solution's increment.
            this->old_solution(Phi) -= lsys.sol.block(LinearSystem::fem_pot_idx);
            l2_pot_incr = old_solution(Phi).l2_norm();
            this->old_solution(Phi) = lsys.sol.block(LinearSystem::fem_pot_idx);
        }

        std::cout // << " l2 incr cast : " << l2_cats_incr << ", neut : " << l2_neut_incr
                // << ", ans : " << l2_ans_incr
                << ", pot : " << l2_pot_incr << std::endl;
    }

    lsys.block_prec.clear();

    // After solving copy the solution components back
    // TO DO: no way to do this with pointers, s.t. we don't need all the ifs?

    if (LinearSystem::fem_pot_idx >= 0)
        this->new_solution(Phi) = lsys.sol.block(LinearSystem::fem_pot_idx);

    if (LinearSystem::bem_pot_idx >= 0) {
        fem_bem.solution_num = lsys.sol.block(LinearSystem::bem_pot_idx);

        for (unsigned int K = 0; K < fem_bem.n_bc_dofs; ++K)
        {
            new_solution(Phi_derivative)(bem_form->dof_to_boundary_mapping.size()
                                          -
                                          bem_form->n_bc_dofs + K) = lsys.sol.block(LinearSystem::bem_pot_idx)(K);
        }
    }

    std::cout << __FUNCTION__ << " end" << std::endl;

    new_solution.distribute_constraints(constraints, all_constraints);
    new_solution.distribute_constraints(hanging_node_constraints, hnc_constraints);
#endif
    this->compute_L2_bem_error();
}


// @sect4{Postprocessing}
//
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
template <int dim>
void FEMBEMPoissonProblem<dim>::refine_grid_global ()
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
    soltrans.interpolate(sol_on_old_grid.data(), sol_on_new_grid.data());

    this->new_solution.reinit(mg_dof_handler.n_dofs());
}

template <int dim>
void FEMBEMPoissonProblem<dim>::refine_grid_adaptive ()
{
    dealii::Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    dealii::KellyErrorEstimator<dim>::estimate (static_cast<dealii::DoFHandler<dim>&>(mg_dof_handler),
                                                dealii::QGauss<dim-1>(3),
                                                typename dealii::FunctionMap<dim>::type(),
                                                new_solution(Phi),
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
}

template <int dim>
void FEMBEMPoissonProblem<dim>::output_results (const unsigned int cycle) const
{

    dealii::DataOut<dim> data_out;

    const dealii::DoFHandler<dim> & dh = mg_dof_handler;

    data_out.attach_dof_handler (dh/*mg_dof_handler*/);

    data_out.add_data_vector (this->new_solution(Phi), "potential");
    data_out.add_data_vector (fem_bem.solution_ref, "pot_ref_sol");

    data_out.build_patches(mapping, this->degree);

    std::ostringstream filename;
    filename << "solution-"
             << fem_bem.testcase_name.c_str()
             << "-q" << numerics_prms.fe_degree
             << "-m" << numerics_prms.fe_mapping_degree
             << "-p_bem" << numerics_prms.bem_quad_order_1
             << "-" << cycle
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);


    // error of FEM solution
    if (true)
    {
        // potential
        dealii::Vector<float> difference_per_cell (this->triangulation.n_active_cells());
        dealii::VectorTools::integrate_difference (mapping,
                                                   mg_dof_handler,
                                                   this->new_solution(Phi), //unit_test_solution_Pot_fem_bem(),
                                                   *fem_bem.u_ref,
                                                   difference_per_cell,
                                                   dealii::QGauss<dim>(2*fe.degree+1),
                                                   dealii::VectorTools::L2_norm);
        const double L2_error = difference_per_cell.l2_norm();

        difference_per_cell = 0.;
        dealii::VectorTools::integrate_difference (mapping,
                                                   mg_dof_handler,
                                                   this->new_solution(Phi), //.unit_test_solution_Pot_fem_bem(),
                                                   *fem_bem.u_ref,
                                                   difference_per_cell,
                                                   dealii::QGauss<dim>(2*fe.degree+1),
                                                   dealii::VectorTools::H1_norm);
        const double H1_error = difference_per_cell.l2_norm();

        std::cout << "FEM, dealii::VectorTools::integrate_difference : final Pot ||u - u_h ||_L2 : " << L2_error
                  << " ||u - u_h ||_H1 : " << H1_error << std::endl;

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

        std::cout << "BEM: ||u - u_h ||_L2_surf : " << fem_bem.L2_diff_u_h
                  << ", ||t - t_h ||_L2_surf : " << fem_bem.L2_diff_dnu_h << std::endl;
    }


    // For debugging purposes dump matrices to screen
    if (false) {
        //    printf("HALLO SL_MATRIX:\n");
        //    fem_bem.SL_matrix.print_formatted(std::cout,
        //                                      3, // const unsigned int 	precision = 3,
        //                                      false, // const bool 	scientific = true,
        //                                      0, // const unsigned int 	width = 0,
        //                                      " * ", // const char * 	zero_string = " ",
        //                                      1., //const double 	denominator = 1.,
        //                                      0.); // const double 	threshold = 0.);

        //    printf("HALLO DL_MATRIX:\n");
        //    fem_bem.DL_matrix.print_formatted(std::cout,
        //                                      3, // const unsigned int 	precision = 3,
        //                                      false, // const bool 	scientific = true,
        //                                      0, // const unsigned int 	width = 0,
        //                                      " * ", // const char * 	zero_string = " ",
        //                                      1., //const double 	denominator = 1.,
        //                                      0.); // const double 	threshold = 0.);
    }
}


// @sect4{FEMBEMPoissonProblem::run}

// Like several of the functions above, this
// is almost exactly a copy of of the
// corresponding function in step-6. The only
// difference is the call to
// <code>assemble_multigrid</code> that takes
// care of forming the matrices on every
// level that we need in the multigrid
// method.
template <int dim>
void FEMBEMPoissonProblem<dim>::run ()
{

    // Then, we prepare the finite element space and the initial conditions
    this->mg_dof_handler.distribute_dofs (fe);

    {
        this->old_solution.reinit(this->mg_dof_handler.n_dofs());
        this->new_solution.reinit(mg_dof_handler.n_dofs());


        std::cout << "Cycle " << 0 << ':' << std::endl;
        setup_system ();

        // cycle 0
        lsys.reinit(fem_bem.n_fem_dofs, fem_bem.n_bc_dofs);

        lsys.setup_block_system(this->matrices, this->fem_bem,
                                this->phys_data
                                , block_pattern, new_solution);

        assemble_system_generic (0/*cycle*/, this->numerics_prms.arch);

        assemble_multigrid ();

        block_solve_fix_point ();

        output_results (0/*cycle*/);
    }


    for (unsigned int cycle=1; cycle<this->numerics_prms.n_mesh_refinements; ++cycle)
    {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        refine_grid_global ();

        setup_system ();

        lsys.reinit(fem_bem.n_fem_dofs, fem_bem.n_bc_dofs);

        dealii::Timer timer;
        timer.restart();

        assemble_system_generic (cycle, this->numerics_prms.arch);

        std::cerr << "assemble_system took " << timer.wall_time() << " seconds."
                  << std::endl;

        //compute_bem_assembly_errors();

        assemble_multigrid ();

        block_solve_fix_point ();

        output_results (cycle);
    }


    // Once the computation is done, write convergence table which contains the profiling results.
#ifndef nUSE_FEM

    std::string conv_results_folder("conv_results");
    std::ostringstream filename;
    QDir conv_results;

    conv_results.mkdir(conv_results_folder.c_str());

    filename << conv_results_folder << "/"
             << fem_bem.testcase_name.c_str()
             << "fem_"
             << "bem_conv_"
             << "_q" << numerics_prms.fe_degree
             << "_m" << numerics_prms.fe_mapping_degree
             << "_p_bem" << numerics_prms.bem_quad_order_1
             << ".out" << std::ends;

    std::ofstream out (filename.str().c_str());

    out << "# " <<           ", FE degree : " << numerics_prms.fe_degree
        << ", MAPPING order : " << numerics_prms.fe_mapping_degree
        << ", BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_1 << std::endl;

    // this->unit_test_table.write_text(out);

    out  << std::endl;
    out  << std::endl;


    static const bool eval_red_rate = false;

    if (eval_red_rate) {
        fem_bem.convergence_table
                .evaluate_convergence_rates("L2(Phi)", dealii::ConvergenceTable::reduction_rate_log2);
        //        fem_bem.convergence_table
        //                .evaluate_convergence_rates("L_infty(Phi)", ConvergenceTable::reduction_rate_log2);

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
              << ", BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order_1
                 // <<  ", SG_QUAD_ORDER : " << SQ_ORDER
              << std::endl;

    fem_bem.convergence_table.write_text(std::cout);

    fem_bem.convergence_table.write_text(out);
#else
    std::cout << "init ref : " << numerics_prms.n_init_refinements
              << ", FE degree : " << numerics_prms.fe_degree
              << ", MAPPING order : " << numerics_prms.fe_mapping_degree
              << ", BEM_QUAD_ORDER : " << numerics_prms.bem_quad_order
                 // <<  ", SG_QUAD_ORDER : " << SQ_ORDER
              << std::endl;
#endif // USE_FEM

    std::ofstream timing_out( numerics_prms.timing_file, std::ios_base::app | std::ios_base::out);

    timing_table.write_text(timing_out);

    timing_out.flush();
    timing_out.close();
}


// #define nUSE_LAC_PC
} // END namespace Step16


