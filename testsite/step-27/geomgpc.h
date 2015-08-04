#ifndef GEOMGPC_H
#define GEOMGPC_H

// These, now, are the include necessary for
// the multi-level methods. The first two
// declare classes that allow us to enumerate
// degrees of freedom not only on the finest
// mesh level, but also on intermediate
// levels (that's what the MGDoFHandler class
// does) as well as allow to access this
// information (iterators and accessors over
// these cells).
//
// The rest of the include files deals with
// the mechanics of multigrid as a linear
// operator (solver or preconditioner).
#include <deal.II/lac/precondition.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/lac/iterative_inverse.h>

namespace step27 {

// #define nUSE_LAC_PC

// @sect4{struct: MGData}
template<int dim>
struct MGData {
    // The following four objects are the
    // only additional member variables,
    // compared to step-6. They first three
    // represent the
    // operators that act on individual
    // levels of the multilevel hierarchy,
    // rather than on the finest mesh as do
    // the objects above while the last object
    // stores information about the boundary
    // indices on each level and information
    // about indices lying on a refinement
    // edge between two different refinement
    // levels.
    //
    // To facilitate having objects on each
    // level of a multilevel hierarchy,
    // deal.II has the MGLevelObject class
    // template that provides storage for
    // objects on each level. What we need
    // here are matrices on each level, which
    // implies that we also need sparsity
    // patterns on each level. As outlined in
    // the @ref mg_paper, the operators
    // (matrices) that we need are actually
    // twofold: one on the interior of each
    // level, and one at the interface
    // between each level and that part of
    // the domain where the mesh is
    // coarser. In fact, we will need the
    // latter in two versions: for the
    // direction from coarse to fine mesh and
    // from fine to coarse. Fortunately,
    // however, we here have a self-adjoint
    // problem for which one of these is the
    // transpose of the other, and so we only
    // have to build one; we choose the one
    // from coarse to fine.
    dealii::MGLevelObject<dealii::SparsityPattern>       sparsity_patterns;

    dealii::MGLevelObject<dealii::SparseMatrix<double> > matrices;

    // Some holds for the multigrid constraints.
    dealii::MGConstrainedDoFs constrained_dofs;

    dealii::MGLevelObject<dealii::SparseMatrix<double> > interface_matrices;


    // @sect5{Function: reinit_Dir_bc}
    void reinit_Dir_bc( const dealii::MGDoFHandler<dim> & mg_dof_handler,
                        typename dealii::FunctionMap<dim>::type & dirichlet_boundary)
    {
        // The multigrid constraints have to be
        // initialized. They need to know about
        // the boundary values as well, so we
        // pass the <code>dirichlet_boundary</code>
        // here as well.
        constrained_dofs.clear();
        constrained_dofs.initialize(mg_dof_handler, dirichlet_boundary);

        resize(mg_dof_handler);
    }

    // @sect5{Function: reinit_vN_bc}
    void reinit_vN_bc( const dealii::MGDoFHandler<dim> & mg_dof_handler)
    {
        constrained_dofs.clear();
        constrained_dofs.initialize(mg_dof_handler);

        resize(mg_dof_handler);
    }

private:
    // @sect5{Function: resize}
    void resize(const dealii::MGDoFHandler<dim> & mg_dof_handler)
    {
        const unsigned int n_levels = mg_dof_handler.get_tria().n_levels();
        // Use tmp objects to anihilate any references to the mg sparsity patterns
        // which causes Subscriptor trouble.
       {
            dealii::MGLevelObject<dealii::SparseMatrix<double> > tmp;
            std::swap(tmp, interface_matrices);

            dealii::MGLevelObject<dealii::SparseMatrix<double> > tmp2;
            std::swap(tmp2, matrices);
        }
        // Now for the things that concern the
        // multigrid data structures. First, we
        // resize the multi-level objects to hold
        // matrices and sparsity patterns for every
        // level. The coarse level is zero (this is
        // mandatory right now but may change in a
        // future revision). Note that these
        // functions take a complete, inclusive
        // range here (not a starting index and
        // size), so the finest level is
        // <code>n_levels-1</code>.  We first have
        // to resize the container holding the
        // SparseMatrix classes, since they have to
        // release their SparsityPattern before the
        // can be destroyed upon resizing.
        interface_matrices.resize(0, n_levels-1);
        sparsity_patterns.resize(0, n_levels-1);
        matrices.resize(0, n_levels-1);
        matrices.clear ();

        // Now, we have to provide a matrix on each
        // level. To this end, we first use the
        // MGTools::make_sparsity_pattern function
        // to first generate a preliminary
        // compressed sparsity pattern on each
        // level (see the @ref Sparsity module for
        // more information on this topic) and then
        // copy it over to the one we really
        // want. The next step is to initialize
        // both kinds of level matrices with these
        // sparsity patterns.
        //
        // It may be worth pointing out that the
        // interface matrices only have entries for
        // degrees of freedom that sit at or next
        // to the interface between coarser and
        // finer levels of the mesh. They are
        // therefore even sparser than the matrices
        // on the individual levels of our
        // multigrid hierarchy. If we were more
        // concerned about memory usage (and
        // possibly the speed with which we can
        // multiply with these matrices), we should
        // use separate and different sparsity
        // patterns for these two kinds of
        // matrices.
        for (unsigned int level=0; level<n_levels; ++level)
        {
            dealii::CompressedSparsityPattern csp;
            csp.reinit(mg_dof_handler.n_dofs(level),
                       mg_dof_handler.n_dofs(level));
            dealii::MGTools::make_sparsity_pattern(mg_dof_handler, csp, level);

            sparsity_patterns[level].copy_from (csp);

            interface_matrices[level].reinit(sparsity_patterns[level]);

            matrices[level].reinit(sparsity_patterns[level]);
        }
    }
};

template<int dim>
struct GeoMGPC {

    dealii::FullMatrix<double> coarse_matrix;
    dealii::MGCoarseGridHouseholder<> coarse_grid_solver;

    typedef dealii::SparseMatrix<double> SpMatrix;

    typedef dealii::Vector<double>       Vc;
#ifndef nUSE_LAC_PC
    typedef // PreconditionSOR<SparseMatrix<double> > //
   // IterativeInverseWrapper<Vc>

     dealii::PreconditionChebyshev<SpMatrix, Vc> SmootherType;

    typedef dealii::MGSmootherPrecondition<SpMatrix, SmootherType, Vc> Smoother;

    std::vector<dealii::ReductionControl> mg_red_control;

#else
    typedef dealii::PreconditionSOR<SpMatrix> SmootherType; // TODO: replace by Krylov method

    typedef dealii::MGSmootherRelaxation<SpMatrix, SmootherType, Vc> Smoother;
#endif

    typedef dealii::PreconditionMG<dim, dealii::Vector<double>,
               //  BoundaryDoFTools::FemBem
                 dealii::MGTransferPrebuilt<dealii::Vector<double> > > FemPC;

    dealii::GrowingVectorMemory<>   vector_memory;

    Smoother smoother;

    dealii::MGTransferPrebuilt<dealii::Vector<double> > transfer;

    dealii::mg::Matrix<> interface_up;
    dealii::mg::Matrix<> interface_down;
    dealii::mg::Matrix<> mg_matrix;

     dealii::Multigrid<dealii::Vector<double> > mg;



     typename GeoMGPC<dim>::FemPC preconditioner;

     GeoMGPC(const MGData<dim> & mg_data,
                       const dealii::MGDoFHandler<dim> & mg_dof_handler,
                       const dealii::ConstraintMatrix & hanging_node_constraints,
                       const uint n_smoothing_steps)
        :
          smoother(vector_memory,
#ifndef nUSE_LAC_PC
                   1, // we always use one smoothing step as the iterative solver itself has a step parameter
                   false // TODO: true -> we take advantage of running the smoother on coarser grid more times
#else
                   n_smoothing_steps,
                   false /*use variable smoothing, i.e. double
                              steps for each coarsening step*/
#endif
                   ),
          transfer(hanging_node_constraints,
                   mg_data.constrained_dofs),
          interface_up(mg_data.interface_matrices),
          interface_down(mg_data.interface_matrices),
                mg_matrix(mg_data.matrices),
          mg(mg_dof_handler,
             mg_matrix,
             coarse_grid_solver,
             transfer,
             smoother,
             smoother),
          preconditioner(mg_dof_handler, mg, transfer)
    {
         reinit(mg_data, mg_dof_handler);

             smoother.set_steps(n_smoothing_steps);
    }

     void reinit(const MGData<dim> & mg_data, const dealii::MGDoFHandler<dim> & mg_dof_handler)
     {
         transfer.build_matrices(mg_dof_handler);

         coarse_matrix.copy_from (mg_data.matrices[0]);
         coarse_grid_solver.initialize (coarse_matrix);
         #ifndef nUSE_LAC_PC
         smoother.initialize(mg_data.matrices, typename SmootherType::AdditionalData( 11, 100., false /* non zero starting. According to doc true should be better for smoothing purposes in MG. At least for Poisson this is not the case. */
                                                                                     )
                             );
 #else
           smoother.initialize(mg_data.matrices);
 #endif
     //    smoother.set_steps(n_smoothing_steps);
 #ifndef nUSE_LAC_PC
         mg_red_control.resize( mg_dof_handler.get_tria().n_levels() );
         for (uint l = 0; l < mg_dof_handler.get_tria().n_levels(); l++)
         {
             // With Liesen's polynomial preconditioner we set the number of iterations.
             mg_red_control[l].set_max_steps(20);
           //  smoother.smoothers[l].initialize(mg_data.matrices[l], // typename
            //                                  PreconditionIdentity // SmootherType ::AdditionalData
              //                                () );
             // smoother.smoothers[l].solver.select("bicgstab");
            // smoother.smoothers[l].solver.set_control(mg_red_control[l]);
         }
 #else
         smoother.set_symmetric(true);
 #endif
         mg.set_edge_matrices(interface_down, interface_up);
     }


};


template<typename Vector>
struct IterativeInverseWrapper : private dealii::IterativeInverse<Vector> {

    typedef dealii::IterativeInverse<Vector> Base;

    struct AdditionalData {

        AdditionalData() {}
    };

    IterativeInverseWrapper() : Base() {}

    template<class MATRIX> // , class PRECONDITION>
    void initialize_iter_inv (const MATRIX & A, const AdditionalData & P = AdditionalData() )
    {
        this->Base::initialize(A, dealii::PreconditionIdentity() );
        red_control.set_max_steps(100);
        red_control.set_tolerance(.1);
        // red_control.log_history(true);
        // red_control.set_reduction(.01);
        this->solver.select("bicgstab");
        this->solver.set_control(red_control);
    }

    template<class VECTOR >
    void vmult (VECTOR & dst, const VECTOR & src) const { this->Base::vmult(dst, src); }

    template<class VECTOR >
    void Tvmult (VECTOR & dst, const VECTOR & src) const { AssertThrow (false, dealii::ExcNotImplemented() ); }

    template<class VECTOR >
    void vmult_add (VECTOR & dst, const VECTOR & src) const { AssertThrow (false, dealii::ExcNotImplemented() ); }

    template<class VECTOR >
    void Tvmult_add (VECTOR & dst, const VECTOR & src) const { AssertThrow (false, dealii::ExcNotImplemented() ); }


    void clear() { this->Base::clear(); }

    dealii::ReductionControl red_control;

};


} // END namespace 27

#endif // GEOMGPC_H
