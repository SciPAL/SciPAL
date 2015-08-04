#ifndef BOUNDARYDOFTOOLS_HH
#define BOUNDARYDOFTOOLS_HH

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_q.h>

namespace BoundaryDoFTools {

using namespace dealii;

template <int dim, int spacedim>
void
map_boundary_dofs_to_support_points (const Mapping<dim,spacedim>     &mapping,
                                     const DoFHandler<dim,spacedim>  &dof_handler,
                                     const std::set<unsigned char>   &boundary_indicators,
                                     const std::vector<unsigned int> &dof_to_boundary_mapping,
                                     std::map<unsigned char,
                                     std::vector<Point<spacedim> > >  &support_points)
{
    // We only need the points on a face
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_face;

    // check whether fe has support
    // points
    Assert (dof_handler.get_fe().has_support_points(),
            typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
    //AssertDimension (support_points.size(), dof_handler.n_dofs());

    // now loop over all cells and
    // enquire the support points on
    // each of these. use a dummy
    // quadrature formula where the
    // quadrature points are located at
    // the unit support points to
    // enquire the location of the
    // support points in real space
    //
    // the weights of the quadrature
    // rule are set to invalid values
    // by the used constructor.
    Quadrature<dim-1> q_dummy(dof_handler.get_fe().get_unit_face_support_points());
    FEFaceValues<dim,spacedim> fe_values (mapping,
                                          dof_handler.get_fe(),
                                          q_dummy,
                                          update_quadrature_points);
    {
        support_points.clear();

        std::set<unsigned char>::const_iterator
                b_id = boundary_indicators.begin(),
                endb = boundary_indicators.end();

        for (; b_id != endb; ++b_id)
        {
            std::set<unsigned char> tmp;
            tmp.insert(*b_id);
            support_points[*b_id].resize(dof_handler.n_boundary_dofs(tmp) );
        }
    }

    typename DoFHandler<dim,spacedim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    std::set<unsigned char>::const_iterator
            b_id = boundary_indicators.begin(),
            endb = boundary_indicators.end();

    for (; cell!=endc; ++cell) // Could also be done with a FilteredIterator
    {
        if (cell->at_boundary() )
            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
            {
                b_id = boundary_indicators.begin();

                for (; b_id != endb; ++b_id )
                    if(cell->face(f)->boundary_indicator() == *b_id )
                    {
                        fe_values.reinit (cell, f); //cell);
                        cell->face(f)->get_dof_indices (local_dof_indices);
                        const std::vector<Point<spacedim> > & points
                                = fe_values.get_quadrature_points ();

                        for (unsigned int i=0; i<dofs_per_cell; ++i)
                            support_points[*b_id][dof_to_boundary_mapping[local_dof_indices[i]]] = points[i];
                    }
            }
    }
}

#ifdef uSE_IMPL
//template <typename VECTOR> class Impl {};

//// For FEMBEM prebuilt transfer matrices have to be reimplemented from scratch.
//template <typename VECTOR>
//class FemBemMGTransferPrebuilt : public ::dealii::MGTransferBase<VECTOR>
//{
// Impl<VECTOR> impl;

//public:
//    /**
//       * Constructor without constraint
//       * matrices. Use this constructor
//       * only with discontinuous finite
//       * elements or with no local
//       * refinement.
//       */
//     FemBemMGTransferPrebuilt () : impl() {}
//      /**
//       * Constructor with constraint matrices as well as mg_constrained_dofs.
//       */
// FemBemMGTransferPrebuilt (const ConstraintMatrix& c,
// const MGConstrainedDoFs& mg_c)
////     :
////       impl(c, mg_c)
// //      constraints(&c),
// //      mg_constrained_dofs(&mg_c)
//    {}


//      /**
//       * Actually build the prolongation
//       * matrices for each level.
//       */
// template <int dim, int spacedim>
// void build_matrices (const MGDoFHandler<dim,spacedim> &mg_dof) { // impl.build_matrices(mg_dof);
// }



// virtual void prolongate (const unsigned int to_level,
//                          VECTOR&            dst,
//                          const VECTOR&      src) const  { // impl.prolongate(to_level, dst, src);
// }

// virtual void restrict_and_add (const unsigned int from_level,
//                                VECTOR&            dst,
//                                const VECTOR&      src) const
// { // impl.restrict_and_add(from_level, dst, src);
// }


// template <int dim, class InVector, int spacedim>
// void
// copy_to_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
// MGLevelObject<VECTOR>& dst,
// const InVector& src) const
// {
// // TO DO
//    // this->Base::copy_to_mg(mg_dof, dst, src);
// }


// template <int dim, class OutVector, int spacedim>
// void
// copy_from_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
// OutVector& dst,
// const MGLevelObject<VECTOR> &src) const
// {
// // TO DO
//   //  this->Base::copy_from_mg(mg_dof, dst, src);
// }

// template <int dim, class OutVector, int spacedim>
// void
// copy_from_mg_add (const MGDoFHandler<dim,spacedim>& mg_dof,
// OutVector& dst,
// const MGLevelObject<VECTOR>& src) const
// {
// // TO DO
//     this->Base::copy_from_mg_add(mg_dof, dst, src);
// }



//};


//template <typename VECTOR>
//class Impl
//        :
//        private ::dealii::MGTransferPrebuilt<VECTOR> {

//    typedef ::dealii::MGTransferPrebuilt<VECTOR> Base;

//public:

//    /**
//      * Constructor without constraint
//      * matrices. Use this constructor
//      * only with discontinuous finite
//      * elements or with no local
//      * refinement.
//      */
//    Impl () : Base() {}
//     /**
//      * Constructor with constraint matrices as well as mg_constrained_dofs.
//      */
//Impl (const ConstraintMatrix& c,
//const MGConstrainedDoFs& mg_c)
//    :
//      Base(c, mg_c)
////      constraints(&c),
////      mg_constrained_dofs(&mg_c)
//   {}


//     /**
//      * Actually build the prolongation
//      * matrices for each level.
//      */
//template <int dim, int spacedim>
//void build_matrices (const MGDoFHandler<dim,spacedim> &mg_dof);

//virtual void prolongate (const unsigned int    to_level,
//VECTOR       &dst,
//const VECTOR &src) const  { this->Base::prolongate(to_level, dst, src); }

//virtual void restrict_and_add (const unsigned int    from_level,
//   VECTOR       &dst,
//   const VECTOR &src) const
//{ this->Base::restrict_and_add(from_level, dst, src); }


//     /**
//      * Transfer from a vector on the
//      * global grid to vectors defined
//      * on each of the levels
//      * separately, i.a. an @p MGVector.
//      */
//template <int dim, class InVector, int spacedim>
//void
//copy_to_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
//MGLevelObject<VECTOR>& dst,
//const InVector& src) const
//{
//// TO DO
//    this->Base::copy_to_mg(mg_dof, dst, src);
//}

//     /**
//      * Transfer from multi-level vector to
//      * normal vector.
//      *
//      * Copies data from active
//      * portions of an MGVector into
//      * the respective positions of a
//      * <tt>Vector<number></tt>. In order to
//      * keep the result consistent,
//      * constrained degrees of freedom
//      * are set to zero.
//      */
//template <int dim, class OutVector, int spacedim>
//void
//copy_from_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
//OutVector& dst,
//const MGLevelObject<VECTOR> &src) const
//{
//// TO DO
//    this->Base::copy_from_mg(mg_dof, dst, src);
//}

//     /**
//      * Add a multi-level vector to a
//      * normal vector.
//      *
//      * Works as the previous
//      * function, but probably not for
//      * continuous elements.
//      */
//template <int dim, class OutVector, int spacedim>
//void
//copy_from_mg_add (const MGDoFHandler<dim,spacedim>& mg_dof,
//OutVector& dst,
//const MGLevelObject<VECTOR>& src) const
//{
//// TO DO
//    this->Base::copy_from_mg_add(mg_dof, dst, src);
//}

//     /**
//      * If this object operates on
//      * BlockVector objects, we need
//      * to describe how the individual
//      * vector components are mapped
//      * to the blocks of a vector. For
//      * example, for a Stokes system,
//      * we have dim+1 vector
//      * components for velocity and
//      * pressure, but we may want to
//      * use block vectors with only
//      * two blocks for all velocities
//      * in one block, and the pressure
//      * variables in the other.
//      *
//      * By default, if this function
//      * is not called, block vectors
//      * have as many blocks as the
//      * finite element has vector
//      * components. However, this can
//      * be changed by calling this
//      * function with an array that
//      * describes how vector
//      * components are to be grouped
//      * into blocks. The meaning of
//      * the argument is the same as
//      * the one given to the
//      * DoFTools::count_dofs_per_component
//      * function.
//      */
//void
//set_component_to_block_map (const std::vector<unsigned int> &map)
//{
//// TO DO
//    this->Base::set_component_to_block_map(map);
//}

//     /**
//      * Memory used by this object.
//      */
//std::size_t memory_consumption () const { return this->Base::memory_consumption(); }
//};

//// Copy from deal.II-7.2.0
//// needed for MG preconditioning of BEM matrices.
//template <typename VECTOR>
//template <int dim, int spacedim>
//void Impl<VECTOR>::build_matrices (
//  const MGDoFHandler<dim,spacedim>  &mg_dof)
//{
//  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
//  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;

//  this->Base::sizes.resize(n_levels);
//  for (unsigned int l=0;l<n_levels;++l)
//    this->Base::sizes[l] = mg_dof.n_dofs(l);

//    // reset the size of the array of
//    // matrices. call resize(0) first,
//    // in order to delete all elements
//    // and clear their memory. then
//    // repopulate these arrays
//    //
//    // note that on resize(0), the
//    // shared_ptr class takes care of
//    // deleting the object it points to
//    // by itself
//  this->Base::prolongation_matrices.resize (0);
//  this->Base::prolongation_sparsities.resize (0);

//  for (unsigned int i=0; i<n_levels-1; ++i)
//    {
//      this->Base::prolongation_sparsities.push_back
//        (std_cxx1x::shared_ptr<SparsityPattern> (new SparsityPattern));
//      this->Base::prolongation_matrices.push_back
//        (std_cxx1x::shared_ptr<SparseMatrix<double> > (new SparseMatrix<double>));
//    }

//    // two fields which will store the
//    // indices of the multigrid dofs
//    // for a cell and one of its children
//  std::vector<unsigned int> dof_indices_parent (dofs_per_cell);
//  std::vector<unsigned int> dof_indices_child (dofs_per_cell);

//    // for each level: first build the sparsity
//    // pattern of the matrices and then build the
//    // matrices themselves. note that we only
//    // need to take care of cells on the coarser
//    // level which have children
//  for (unsigned int level=0; level<n_levels-1; ++level)
//    {

//        // reset the dimension of the structure.
//        // note that for the number of entries
//        // per row, the number of parent dofs
//        // coupling to a child dof is
//        // necessary. this, of course, is the
//        // number of degrees of freedom per
//        // cell
//        // increment dofs_per_cell
//        // since a useless diagonal
//        // element will be stored
//      CompressedSimpleSparsityPattern csp (this->Base::sizes[level+1],
//                                           this->Base::sizes[level]);
//      std::vector<unsigned int> entries (dofs_per_cell);
//      for (typename MGDoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
//           cell != mg_dof.end(level); ++cell)
//        if (cell->has_children())
//          {
//            cell->get_mg_dof_indices (dof_indices_parent);

//            Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
//                   ExcNotImplemented());
//            for (unsigned int child=0; child<cell->n_children(); ++child)
//              {
//                  // set an alias to the
//                  // prolongation matrix for
//                  // this child
//                const FullMatrix<double> &prolongation
//                  = mg_dof.get_fe().get_prolongation_matrix (child,
//                                                             cell->refinement_case());

//                Assert (prolongation.n() != 0, ExcNoProlongation());

//                cell->child(child)->get_mg_dof_indices (dof_indices_child);

//                  // now tag the entries in the
//                  // matrix which will be used
//                  // for this pair of parent/child
//                for (unsigned int i=0; i<dofs_per_cell; ++i)
//                  {
//                    entries.resize(0);
//                    for (unsigned int j=0; j<dofs_per_cell; ++j)
//                      if (prolongation(i,j) != 0)
//                        entries.push_back (dof_indices_parent[j]);
//                    csp.add_entries (dof_indices_child[i],
//                                     entries.begin(), entries.end());
//                  }
//              }
//          }

//      this->Base::prolongation_sparsities[level]->copy_from (csp);
//      csp.reinit(0,0);
//      this->Base::prolongation_matrices[level]->reinit (*this->Base::prolongation_sparsities[level]);

//        // now actually build the matrices
//      for (typename MGDoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
//           cell != mg_dof.end(level); ++cell)
//        if (cell->has_children())
//          {
//            cell->get_mg_dof_indices (dof_indices_parent);

//            Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
//                   ExcNotImplemented());
//            for (unsigned int child=0; child<cell->n_children(); ++child)
//              {
//                  // set an alias to the
//                  // prolongation matrix for
//                  // this child
//                const FullMatrix<double> &prolongation
//                  = mg_dof.get_fe().get_prolongation_matrix (child,
//                                                             cell->refinement_case());

//                cell->child(child)->get_mg_dof_indices (dof_indices_child);

//                  // now set the entries in the
//                  // matrix
//                for (unsigned int i=0; i<dofs_per_cell; ++i)
//                  this->Base::prolongation_matrices[level]->set (dof_indices_child[i],
//                                                     dofs_per_cell,
//                                                     &dof_indices_parent[0],
//                                                     &prolongation(i,0),
//                                                     true);
//              }
//          }
//    }


//                                // impose boundary conditions
//                                // but only in the column of
//                                // the prolongation matrix
//  if (this->Base::mg_constrained_dofs != 0)
//  if (this->Base::mg_constrained_dofs->set_boundary_values())
//    {
//      std::vector<unsigned int> constrain_indices;
//      for (int level=n_levels-2; level>=0; --level)
//        {
//          if (this->Base::mg_constrained_dofs->get_boundary_indices()[level].size() == 0)
//            continue;

//                                // need to delete all the columns in the
//                                // matrix that are on the boundary. to achieve
//                                // this, create an array as long as there are
//                                // matrix columns, and find which columns we
//                                // need to filter away.
//          constrain_indices.resize (0);
//          constrain_indices.resize (this->Base::prolongation_matrices[level]->n(), 0);
//          std::set<unsigned int>::const_iterator dof
//            = this->Base::mg_constrained_dofs->get_boundary_indices()[level].begin(),
//            endd = this->Base::mg_constrained_dofs->get_boundary_indices()[level].end();
//          for (; dof != endd; ++dof)
//            constrain_indices[*dof] = 1;

//          const unsigned int n_dofs = this->Base::prolongation_matrices[level]->m();
//          for (unsigned int i=0; i<n_dofs; ++i)
//            {
//              SparseMatrix<double>::iterator
//                start_row = this->Base::prolongation_matrices[level]->begin(i),
//                end_row   = this->Base::prolongation_matrices[level]->end(i);
//              for(; start_row != end_row; ++start_row)
//                {
//                  if (constrain_indices[start_row->column()] == 1)
//                    start_row->value() = 0;
//                }
//            }
//        }
//    }

//                                // to find the indices that describe the
//                                // relation between global dofs and local
//                                // numbering on the individual level, first
//                                // create a temp vector where the ith level
//                                // entry contains the respective global
//                                // entry. this gives a neat way to find those
//                                // indices. in a second step, actually build
//                                // the std::vector<std::pair<uint,uint> > that
//                                // only contains the active dofs on the
//                                // levels.

//  this->Base::copy_indices.resize(n_levels);
//  std::vector<unsigned int> temp_copy_indices;
//  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
//  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
//  for (int level=mg_dof.get_tria().n_levels()-1; level>=0; --level)
//    {
//      this->Base::copy_indices[level].clear();
//      typename MGDoFHandler<dim,spacedim>::active_cell_iterator
//        level_cell = mg_dof.begin_active(level);
//      const typename MGDoFHandler<dim,spacedim>::active_cell_iterator
//        level_end  = mg_dof.end_active(level);

//      temp_copy_indices.resize (0);
//      temp_copy_indices.resize (mg_dof.n_dofs(level), numbers::invalid_unsigned_int);

//        // Compute coarse level right hand side
//        // by restricting from fine level.
//      for (; level_cell!=level_end; ++level_cell)
//        {
//          DoFAccessor<dim, DoFHandler<dim,spacedim> >& global_cell = *level_cell;
//            // get the dof numbers of
//            // this cell for the global
//            // and the level-wise
//            // numbering
//          global_cell.get_dof_indices(global_dof_indices);
//          level_cell->get_mg_dof_indices (level_dof_indices);

//          for (unsigned int i=0; i<dofs_per_cell; ++i)
//          {
//            if(this->Base::mg_constrained_dofs != 0)
//            {
//              if(!this->Base::mg_constrained_dofs->at_refinement_edge(level,level_dof_indices[i]))
//                temp_copy_indices[level_dof_indices[i]] = global_dof_indices[i];
//            }
//            else
//              temp_copy_indices[level_dof_indices[i]] = global_dof_indices[i];
//          }
//        }

//                                // now all the active dofs got a valid entry,
//                                // the other ones have an invalid entry. Count
//                                // the invalid entries and then resize the
//                                // copy_indices object. Then, insert the pairs
//                                // of global index and level index into
//                                // copy_indices.
//      const unsigned int n_active_dofs =
//        std::count_if (temp_copy_indices.begin(), temp_copy_indices.end(),
//                       std::bind2nd(std::not_equal_to<unsigned int>(),
//                                    numbers::invalid_unsigned_int));
//      this->Base::copy_indices[level].resize (n_active_dofs);
//      unsigned int counter = 0;
//      for (unsigned int i=0; i<temp_copy_indices.size(); ++i)
//        if (temp_copy_indices[i] != numbers::invalid_unsigned_int)
//          this->Base::copy_indices[level][counter++] =
//            std::pair<unsigned int, unsigned int> (temp_copy_indices[i], i);
//      Assert (counter == n_active_dofs, ExcInternalError());
//    }
//}
#endif

} // END namespace BoundaryDoFTools

#endif // BOUNDARYDOFTOOLS_HH
