#ifndef DRS_LAC_HH
#define DRS_LAC_HH

namespace step27 {


// Forward declarations.
struct MtxInfo;

template <int dim> class FemBemData;


// @sect3{struct SolComps}
//
// TODO: fix doc
// // For unknown reasons we cannot use a simple enum AND get sizeof(sc)/sizeof(sc[0]) correctly calculated.
// // As a quick fix we define a structure with static integers which we then can use as named constants.
// From a physical point of view this structure summarizes what DRS is about.
enum SolComps {Cations=1, Anions=2, Neutral=3, Phi=4, Phi_derivative=5 };


// @sect3{struct SparsityPatternSet}
//
// Base class for managing a collection of sparsity patterns.
// The main purpose is to do compression and reinitialization by iterating over the elements of a container.
template<typename PatternIds>
struct SparsityPatternSetBase
        :
        public dealii::Subscriptor,
        public std::map<PatternIds,
        dealii::SparsityPattern>
{
    typedef std::map<PatternIds, dealii::SparsityPattern> Base;

    typedef typename Base::iterator iterator;

    typedef typename Base::const_iterator const_iterator;

    SparsityPatternSetBase()
        :
          Base()
    {}


    dealii::SparsityPattern & operator() (PatternIds p_id)
    {
#ifdef PATTERN_ACCESS_DEBUG
        std::cout << "n patterns : " << this->size() << " p_id : " << int(p_id) << std::endl;
#endif
        iterator p = this->find(p_id);
        Assert(p != this->end(), dealii::ExcMessage("Illegal pattern id") );

        return p->second;
    }

// operator[] does not allow constness. Therefore, we define this wrapper which additionally checks
    // whether the requested object exists.
    const dealii::SparsityPattern & operator() (PatternIds p_id) const
    {
        const_iterator p = this->find(p_id);

#ifdef PATTERN_ACCESS_DEBUG
        std::cout << "n patterns : " << this->size() << " p_id : " << int(p_id) << std::endl;
#endif
        Assert(p != this->end(), dealii::ExcMessage("Illegal pattern id") );

        return p->second;
    }

    // Compress all patterns.
    void compress() {

        iterator p = this->begin(), p_end = this->end();

        for (; p != p_end; ++p)
            p->second.compress();
    }

};


// can be moved into FemBemPoissonProblem. Then we can remove the template parameter.
// This class implements the problem-specific setup of the sparsity patterns.
template<typename PatternInfo>
struct SparsityPatternSet
        :
        public SparsityPatternSetBase<typename PatternInfo::Ids>
{
    typedef SparsityPatternSetBase<typename PatternInfo::Ids> Base;

    typedef typename Base::iterator iterator;

    typedef typename Base::const_iterator const_iterator;

    SparsityPatternSet()
        :
          Base()
    #ifdef USE_EL_MASS_MAT
        ,
          cathode_mass((*this)[PatternInfo::FEMCathodeBoundaryMass /*TO DO: additional Id: user_defined?*/]),
          anode_mass((*this)[PatternInfo::FEMAnodeBoundaryMass /*TO DO: additional Id: user_defined?*/]),
          // The following 2 matrices are of size n_boundary_dofs(bc_id) x n_boundary_dofs(bc_id)
          cathode_mass_bc_dofs((*this)[PatternInfo::BEMCathodeBoundaryMass]),
          anode_mass_bc_dofs((*this)[PatternInfo::BEMAnodeBoundaryMass])
    #endif
    {}


#ifdef USE_EL_MASS_MAT
     // In addition to the lookup by pattern type, we add the possilility to address patterns by their physical name.
 //
    dealii::SparsityPattern & cathode_mass;

    dealii::SparsityPattern & anode_mass;

    dealii::SparsityPattern & cathode_mass_bc_dofs;

    dealii::SparsityPattern & anode_mass_bc_dofs;
#endif

    // Reinitialization does not include compression!
    template<int dim>
    void reinit_FEM(const dealii::DoFHandler<dim> & dof_handler)
    {
        // In order to avoid redundant code we group the patterns
        // according to the make_sparsity_pattern function which they need.
        // dealii::SparsityPattern *

        typename PatternInfo::Ids fem_patterns[2] = { PatternInfo::FEMwDirBC, PatternInfo::FEMwvNBC }; // &this->Phi, &this->Ions };

        for (uint s = 0; s < 2; s++)
        {
            dealii::SparsityPattern &si = (*this)[fem_patterns[s]];
            si.reinit (dof_handler.n_dofs(),
                       dof_handler.n_dofs(),
                       dof_handler.max_couplings_between_dofs());
            dealii::DoFTools::make_sparsity_pattern (dof_handler, si);
        }

        // The blocks coupling the FEM with the BEM part
        // and the boundary mass matrices for the electrodes
        // need a tailor-made setup.

#ifdef USE_EL_MASS_MAT
        // kill any old content. This is saver than calling clear().
        {
            std::map<dealii::types::boundary_id, std::map<uint, uint> >  tmp;
            std::swap(tmp, electrode_bc_local_2_global_map);
        }

        build_boundary_pattern(cathode_mass, cathode_mass_bc_dofs, dof_handler, DRS::BCids::cathode() );

        build_boundary_pattern(anode_mass, anode_mass_bc_dofs, dof_handler, DRS::BCids::anode() );
#endif
    }


#ifndef USE_EL_MASS_MAT
    std::map<dealii::types::boundary_id, std::vector<unsigned int> > dof_2_bc_map;

    std::map<dealii::types::boundary_id, std::map<uint, uint> > electrode_bc_local_2_global_map;

    // The problem with the mass matrices on those parts of the boundary
    // which represent the electrodes, is that we need them in two different flavors.
    // We need a representation as sparse matrices with as many rows and columns
    // as there are DoFs on the subboundary. In order to to use the dealii::BlockMatrixArray class
    // we additionally need the representation as sparse matrices with as many rows and columns as there
    // are DoFs in the whole system (i.e. one scalar PDE).
    // The following function assembles these two types of patterns for a given
    // subboundary, identified by @p bc_id.
    template<int dim>
    void build_boundary_pattern(dealii::SparsityPattern & electrode_mass,
                                dealii::SparsityPattern & electrode_mass_bc_dofs,
                                const dealii::DoFHandler<dim> & dof_handler,
                                dealii::types::boundary_id bc_id)
    {
        std::set< dealii::types::boundary_id > electrode_boundary_id;

        electrode_boundary_id.insert( bc_id); // DRS::BCids::cathode());

        typename dealii::FunctionMap<dim>::type electrode_indicator;
        electrode_indicator[bc_id /*DRS::BCids::cathode()*/] = 0;

        std::vector<unsigned int> & electrode_dof_2_bc_map = dof_2_bc_map[bc_id];

        // After all dofs are in place, we have to update
        // the boundary dof lists of the cathode and anode.
        // This can be considered as the inverse operation of @p extract_boundary_dofs.
        dealii::DoFTools::map_dof_to_boundary_indices(dof_handler,
                                                      electrode_boundary_id,
                                                      electrode_dof_2_bc_map);

        electrode_mass_bc_dofs.reinit(dof_handler.n_boundary_dofs(electrode_indicator),
                                      dof_handler.n_boundary_dofs(electrode_indicator),
                                      dof_handler.max_couplings_between_dofs());

        dealii::DoFTools::make_boundary_sparsity_pattern(
                    dof_handler,
                    electrode_indicator,
                    electrode_dof_2_bc_map,
                    electrode_mass_bc_dofs
                    );
        // Before we can use the temporary sparsity pattern we have to compress it.
        electrode_mass_bc_dofs.compress();

        // TO DO: split into two functions } {

        // Now, build the real ones.
        std::vector<unsigned int> row_lengths_c(dof_handler.n_dofs(), 0);

        // For the interior dofs the row lengths are zero.
        unsigned int bc_dof = 0;
        unsigned int n_visited_bc_dofs = 0;

        std::map<uint, uint> & electrode_bc_local_2_global
                = electrode_bc_local_2_global_map[bc_id];

        for (unsigned int r = 0; r < dof_handler.n_dofs(); r++)
            if (electrode_dof_2_bc_map[r] != dealii::DoFHandler<dim>::invalid_dof_index)
            {
                bc_dof = electrode_dof_2_bc_map[r];
                electrode_bc_local_2_global[bc_dof] = r;
                row_lengths_c[r] = electrode_mass_bc_dofs.row_length(bc_dof);
                n_visited_bc_dofs++;
            }

        Assert(n_visited_bc_dofs == dof_handler.n_boundary_dofs(electrode_indicator),
               dealii::ExcMessage("Bug in electrode boundary enumeration"));

        electrode_mass.reinit(dof_handler.n_dofs(),
                              dof_handler.n_dofs(),
                              row_lengths_c);

        // Finally, add the entries
        // what is this for? : std::map<uint, uint>::const_iterator c_end = electrode_bc_local_2_global.end();
        for (unsigned int r = 0; r < dof_handler.n_dofs(); r++)
            if (electrode_dof_2_bc_map[r] != dealii::DoFHandler<dim>::invalid_dof_index)
            {
                bc_dof = electrode_dof_2_bc_map[r];
                dealii::SparsityPattern::const_iterator
                        col =  electrode_mass_bc_dofs.begin(bc_dof),
                        endc = electrode_mass_bc_dofs.end(bc_dof);
                // std::cout << "global row : " << r << ", local bc dof : " << bc_dof << std::endl;

                for (; col != endc; ++col)
                {
                    std::map<uint, uint>::const_iterator c = electrode_bc_local_2_global.find(col->column());
#ifdef nUSE_DEAL_II_8_1
                    Assert(c != c_end,
                           ExcMessage("Illegal bc dof on electrode boundary"));
#endif
                    unsigned int global_col = c->second;
                    electrode_mass.add(r, global_col);
                }
            }
    }
#endif

    // FEM part of FEM-BEM problem
    std::vector<unsigned int> lower_off_block_col_indices;
    std::vector<unsigned int> lower_off_block_row_indices;

    template<int dim, typename FemBemData>
    void reinit_FEM_BEM(const dealii::DoFHandler<dim> & /*dof_handler*/,
                        FemBemData/*<dim> */& fem_bem)
    {

        // The upper right off-diagonal block is just the mass matrix for the boundary dofs.
        // We create a new sparsity pattern for a rectangular matrix and copy the
        // existing boundary sparsity pattern at the right places.

        dealii::SparsityPattern & FEM_BEM_upper_off_block = (*this)[PatternInfo::FEMBEMUpperOffBlock];
        dealii::SparsityPattern & FEM_BEM_lower_off_block = (*this)[PatternInfo::BEMFEMLowerOffBlock];

        dealii::SparsityPattern & BEM_BEM_block = (*this)[PatternInfo::Dense];


        // As a preparatory step we have to get the lengths of the rows of the boundary mass matrix
        const unsigned int n_fem_dofs  = fem_bem.n_fem_dofs;
        const unsigned int n_bulk_dofs = fem_bem.n_bulk_dofs;
        const unsigned int n_bc_dofs = fem_bem.n_bc_dofs;

        // Do the easy part first.
        // To save memory we set the max. number of entries per row to zero.
        // For full matrices we only need the information about the numer of rows and columns.
        BEM_BEM_block.reinit(n_bc_dofs, n_bc_dofs, 0);

        // At this point we cannot use the element access via operator[]
        // because the pattern has to exist and has to be initialized.
        dealii::SparsityPattern & fem_bem_sp = (*this)(PatternInfo::FEMBEM);

        std::vector<unsigned int> row_lengths(n_fem_dofs, 0);
        // For the interior dofs the row lengths are zero.
        for (unsigned int r = n_bulk_dofs; r < n_fem_dofs; r++)
            row_lengths[r] = fem_bem_sp.row_length(r - n_bulk_dofs);

        FEM_BEM_upper_off_block.reinit(n_fem_dofs, // as many rows as the FEM matrix
                                       fem_bem.n_bc_dofs, // as many cols as the BEM matrices.
                                       row_lengths
                                       );

        // Finally, add the entries
        for (unsigned int r = n_bulk_dofs; r < n_fem_dofs; r++)
        {
            dealii::SparsityPattern::const_iterator
                    col =  fem_bem_sp.begin(r - n_bulk_dofs),
                    endc = fem_bem_sp.end(r - n_bulk_dofs);

            for (; col != endc; ++col)
                FEM_BEM_upper_off_block.add(r, col->column());
        }

        std::vector<unsigned int> bc_vector(n_fem_dofs);
        for(uint i = n_bulk_dofs; i < n_fem_dofs; i++)
            bc_vector[i] = n_bc_dofs;

        //  TODO: fem_bem.hypersingular_to_laplace_sparsity.reinit(n_fem_dofs, n_fem_dofs, bc_vector);


        // Next, we go for the lower off-diagonal block
        FEM_BEM_lower_off_block.reinit( fem_bem.n_bc_dofs, // as many rows as the BEM matrices
                                        n_fem_dofs, // as many cols as the FEM matrix
                                        fem_bem.n_bc_dofs // as many entries per row as in the dense double layer matrix
                                        );
        // The pattern itself is rather simple
        lower_off_block_col_indices.resize(fem_bem.n_bc_dofs);
        lower_off_block_row_indices.resize(fem_bem.n_bc_dofs);

        // Set up the column indices for all rows only once since they are identical.
        // They are just the positions in the dense matrix shifted by the number of interior dofs.
        // This is not really smart but currently there is no other option.
        for (unsigned int c = 0; c < fem_bem.n_bc_dofs; c++)
        {
            lower_off_block_row_indices[c] = c; // needed for the scatter operation
            lower_off_block_col_indices[c] = n_bulk_dofs + c;
        }

        // Here, we only loop over as many rows as there are boundary dofs
        for (unsigned int r = 0; r < fem_bem.n_bc_dofs; r++) {
            FEM_BEM_lower_off_block.add_entries(r,
                                                lower_off_block_col_indices.begin(),
                                                lower_off_block_col_indices.end() );

            // TODO:  fem_bem.hypersingular_to_laplace_sparsity.add_entries(r+n_bulk_dofs, lower_off_block_col_indices.begin(), lower_off_block_col_indices.end());



        }
    }
};


// @sect3{struct SparseMatrixSet}
//
// TODO: rename into MatrixSet
template<typename MtxInfo>
struct SparseMatrixSet
        :
        protected std::map<typename MtxInfo::Ids, dealii::SparseMatrix<double> > {

protected:
    typedef std::map<typename MtxInfo::Ids, dealii::FullMatrix<double> /*TODO: replace by SciPAL::FullMatrixAccessor*/ > FullMatrixSetBase;
    FullMatrixSetBase full_matrices;

public:
    typedef std::map<typename MtxInfo::Ids, dealii::SparseMatrix<double> > Base;
    typedef SparsityPatternSet<typename MtxInfo::SPInfo> PatternSet;
    typedef dealii::SparseMatrix<double> Mtx;

    int n_elements() const { return this->size(); }

    SparseMatrixSet(const PatternSet& sps)
        :
          Base(),
          patterns(&sps)
    {}



   typename FullMatrixSetBase::mapped_type & fm(typename MtxInfo::Ids id)
    {
        // Currently, we can only do sparse matrices
        Assert(id == MtxInfo::SLFullMatrix /*SPInfo::Dense*/, dealii::ExcMessage("Illegal matrix id. fm() can only return dense matrices.") );
        return this->full_matrices[id];
    }

    const typename FullMatrixSetBase::mapped_type & fm(typename MtxInfo::Ids id) const
    {
        // Currently, we can only do sparse matrices
        Assert(id == MtxInfo::SLFullMatrix /*SPInfo::Dense*/, dealii::ExcMessage("Illegal matrix id. fm() can only return dense matrices.") );
       typename FullMatrixSetBase::const_iterator m = this->full_matrices.find(id);
          Assert(m != this->full_matrices.end(), dealii::ExcMessage("Illegal matrix id") );

        return m->second;
    }


    const Mtx & operator() (typename MtxInfo::Ids id) const
    {
        typename Base::const_iterator m = this->find(id);
        Assert(m != this->end(), dealii::ExcMessage("Illegal matrix id") );

        return m->second;
    }


    Mtx & operator() (typename MtxInfo::Ids id) {
        typename Base::iterator m = this->find(id);
        Assert(m != this->end(), dealii::ExcMessage("Illegal matrix id") );

        return m->second;
    }

    // Reinitialize a particular matrix.
    void reinit(typename MtxInfo::Ids id)
    {
        typename Base::iterator m = this->find(id);
        Assert(m != this->end(), dealii::ExcMessage("Illegal matrix id") );

        // Currently, we can only do sparse matrices
        Assert(id != MtxInfo::SPInfo::Dense, dealii::ExcNotImplemented() );

        // The matrix exists, hence we can safely access the pattern map.
        m->second.reinit(
                    (*patterns)(mi.pattern_types[id])
                    );

    }

    // Reinitialize the whole set of matrices.
    void reinit()
    {
        typename MtxInfo::Id2Pattern::const_iterator
                m = mi.pattern_types.begin(),
                m_end = mi.pattern_types.end();

        // To avoid unintentional instantiation of non-existing sparsity
        // patterns we use the non-modifying element access operator() of
        // SparsityPatternSet.
        // The dimensions of full matrices are stored in a sparsity pattern as well.
        // However, the sparsity pattern itself is empty. Hence we have to distinguish
        // whether a matrix is sparse or dense.
        for (; m != m_end; ++m)
            if (m->second != MtxInfo::SPInfo::Dense)
                (*this)[m->first].reinit((*patterns)(m->second));
            else
                this->full_matrices[m->first].reinit((*patterns)(m->second).n_rows(), (*patterns)(m->second).n_cols() );


    }

protected:
    dealii::SmartPointer<const PatternSet> patterns;

    // We need the mapping from the matrices to the sparsity pattern types.
    MtxInfo mi;
};



// @sect4{Class: BlockPatternInfo}
//
// All information about the structure of the PDE system is stored in this type,
// i.e. which physical components interact with each other
// (so to say the sparsity pattern of the coefficient matrix)
// and the strength of the interaction (the coefficient of the block).
// In analogy to dealii::BlockMatrixArray we store a
// list of coefficients and its associated matrix types
// for each entry of the interaction matrix.
// Additionally, we store the analogous information about the block in preconditioner
// associated to this matrix block. Preconditioners are identified via the
// MtxInfo tag as well, since preconditioners behave like matrices.
// By default, a matrix is not preconditioned at all.
template<typename MtxInfo>
struct InteractionInfo {

    InteractionInfo(typename MtxInfo::Ids m,
    double c)
        :
          m_tag(m), mtx_coeff(c), pc_tag(MtxInfo::None), pc_coeff(0.)
    {}

    InteractionInfo(typename MtxInfo::Ids m, double c,
                    typename MtxInfo::Ids pc, double pcc)
        :
          m_tag(m), mtx_coeff(c), pc_tag(pc), pc_coeff(pcc)
    {}

    typename MtxInfo::Ids m_tag;
    double mtx_coeff;

    typename MtxInfo::Ids pc_tag;
    double pc_coeff;

};

typedef std::pair<SolComps, SolComps> IJ;
typedef InteractionInfo<MtxInfo> Interaction;

typedef std::map<IJ, std::vector<Interaction> > BlockPatternInfo;


} // END namespace step27

#endif // DRS_LAC_HH
