#ifndef BEM_FORMS_S27_HH
#define BEM_FORMS_S27_HH

#include <step-27/bemforms.h>

namespace step27
{
// @sect5{Function: assemble_sl_dl}
//
// This function assembles the single-layer and double-layer matrix by collocation.
// @param cell an iterator pointing the cell behind the face
// @p param face_id for which its contribution to the boundary method has to be assembled.
// @param DL_matrix a reference to a dealii::FullMatrix of correct dimensions to store the DL entries.
// @param SL_matrix a reference to a dealii::FullMatrix of correct dimensions to store the SL entries.
// @param arch a flag indicating whether compilation on the cpu or gpu via CUDA should be done.
// If its value is @p both the benchmarking comparing speedup and correctness of the CUDA-based assembly is performed.
template<int dim, typename IntegrationTraits, typename CudaDriver>
template<typename cellIterator>
void BemCollocationForm<dim, IntegrationTraits, CudaDriver>::assemble_sl_dl(
        const cellIterator & cell,
        const uint face_id,
        dealii::FullMatrix<double> & DL_matrix,
        dealii::FullMatrix<double> & SL_matrix,
        const Architecture arch)
{
    // Depending on the template parameter for the
    // architecture we invoke one of the non-public member functions
    // for assembling the single- and double-layer matrices.
    switch(arch) {
    case cpu:
        this->assemble_sl_dl_cpu(cell, face_id, DL_matrix, SL_matrix);
        break;
    case cuda:
        this->assemble_sl_dl_gpu(cell, face_id, DL_matrix, SL_matrix);
        break;
    case both: {
        // We run the assembly on cpu and gpu and compare the results.
        // To this end, we need two additional matrices for storing the
        // result of the additional assembly. Since there is not much which could go wrong
        // on the side of the cpu (no crash of the cuda driver, mysterious memory overwrites in the vram, etc ...),
        // we use the result of the cpu assembly as reference. Since it is a plain deal.II implementation most errors
        // on the cpu side should be filtered out by the various assertions inside of deal.II.
        dealii::FullMatrix<double> DL_matrix_reference (DL_matrix);
        dealii::FullMatrix<double> SL_matrix_reference (SL_matrix);

        this->assemble_sl_dl_cpu(cell, face_id, DL_matrix_reference, SL_matrix_reference);
        this->assemble_sl_dl_gpu(cell, face_id, DL_matrix, SL_matrix);

        // After the assembly we can start to check for potential errors.
        // As a first test, we check the max norm of the difference.
        DL_matrix_reference.add(-1, DL_matrix);
        SL_matrix_reference.add(-1, SL_matrix);

        const double dl_max_error = DL_matrix_reference.linfty_norm();
        const double sl_max_error = SL_matrix_reference.linfty_norm();

        // Since we work with double we accept 5e-16 as maximal error.
        // The output is formatted such that one easily recognizes which matrix is possibly wrong.
        //
        const double tol_assembly = 1e-15;
        if (dl_max_error > tol_assembly )
            std::cout << "errors in matrix assembly: DL: " << dl_max_error;

        if (sl_max_error > tol_assembly)
            std::cout << (dl_max_error < tol_assembly ? "                                     " : " ") << ",     SL : " << sl_max_error;

        if (dl_max_error > tol_assembly  || sl_max_error > tol_assembly)
            std::cout << std::endl;

        // ----------
#ifdef BEM_DEBUGffff
    std::cout << "HALLO SL_MATRIX:" << std::endl;
    this->matrices.fm(MtxInfo::SLFullMatrix).print_formatted(std::cout,
                                      12, // const unsigned int 	precision = 3,
                                      false, // const bool 	scientific = true,
                                      0, // const unsigned int 	width = 0,
                                      " * ", // const char * 	zero_string = " ",
                                      1., //const double 	denominator = 1.,
                                      0.); // const double 	threshold = 0.);

     std::cout << "HALLO DL_MATRIX:" << std::endl;
    fem_bem.DL_matrix.print_formatted(std::cout,
                                      12, // const unsigned int 	precision = 3,
                                      false, // const bool 	scientific = true,
                                      0, // const unsigned int 	width = 0,
                                      " * ", // const char * 	zero_string = " ",
                                      1., //const double 	denominator = 1.,
                                      0.); // const double 	threshold = 0.);
#endif
        // -----------
        break;
    }
    default:
        AssertThrow(false, dealii::ExcNotImplemented());
        break;
    }
}










template<int dim, typename IntegrationTraits, typename CudaDriver>
template<typename cellIterator>
void BemCollocationForm<dim, IntegrationTraits, CudaDriver>::assemble_sl_dl_cpu(const cellIterator & /*cell*/,
                                                                    const uint /*face_id*/,
                                                                    dealii::FullMatrix<double> & DL_matrix,
                                                                    dealii::FullMatrix<double> & SL_matrix)
{
   // std::cout << __FUNCTION__ << std::endl;
    // TO DO
    const std::vector<dealii::Point<dim> > &
            x = support_points.begin()->second;

    for (uint i=0; i < n_bem_points; ++i)
        for (uint a=0; a < this->n_q_points; ++a)
        {
            dealii::Point<dim> R = x[i] - this->q[a];
            double G_ia = LaplaceKernel<dim>::single_layer(R);
            double H_ia = this->normals[a] * LaplaceKernel<dim>::double_layer(R);

            for (uint j = 0; j <this->dofs_per_cell; ++j)
            {
                uint J = this->global_bc_dofs[j];
                // loop over quad points only for valid boundary dofs
                if (J != dealii::DoFHandler<dim>::invalid_dof_index)
                {
                    DL_matrix(i, J) += - H_ia * this->W(a, j); // H_ia * w_aj * JxW[a];
                    SL_matrix(i, J) += /*-*/ + G_ia * this->W(a, j); // -G_ia * w_aj * JxW[a];
                }            }
        }
}
} // namespace step27 END



// @sect5{Function: assemble_sl_dl_gpu}
//
// Same as assemble_sl_dl_cpu, but running on the GPU.
#ifndef GPU_PART
template<int dim, typename IntegrationTraits, typename CudaDriver>
template<typename cellIterator>
void step27::BemCollocationForm<dim, IntegrationTraits, CudaDriver>::assemble_sl_dl_gpu(const cellIterator & cell,
                                    const uint face_id,
                                    dealii::FullMatrix<double> & DL_matrix,
                                    dealii::FullMatrix<double> & SL_matrix)
{
}

#else
template<int dim, typename IntegrationTraits, typename CudaDriver>
template<typename cellIterator>
void step27::BemCollocationForm<dim, IntegrationTraits, CudaDriver>::assemble_sl_dl_gpu(const cellIterator & cell,
                                    const uint face_id,
                                    dealii::FullMatrix<double> & DL_matrix,
                                    dealii::FullMatrix<double> & SL_matrix)
{
   //   std::cout << __FUNCTION__ << std::endl;

    // FIXME: move to attribut list of class.
    FullMatrixAccessor<double> local_DL_matrix_k(0,0, this->is_col_major);
    FullMatrixAccessor<double> local_SL_matrix_k(0,0, this->is_col_major);


    // ...then we need to resize the containers
    local_DL_matrix_k.reinit(this->n_bem_points, this->dofs_per_cell);
    local_SL_matrix_k.reinit(this->n_bem_points, this->dofs_per_cell);

    const std::vector<dealii::Point<dim> > &
            x = support_points.begin()->second;

    // For running the assembly on a CUDA device we have to convert
    // the set of support points,
    this->point_to_double_vector(this->x_k, x);

    // ... quadrature points,
    this->point_to_double_vector(this->q_k, this->q);

    // ... and normals at the quadrature points into
    // one-dimensional arrays which contain the coordinates
    // in lexicographical order.
    this->point_to_double_vector(this->n_k, this->normals);

    // FIXME: what doas reinit do?
    this->CudaDriver::reinit(this->n_bem_points, this->dofs_per_cell);

    this->assemble_bem_matrices(local_DL_matrix_k, local_SL_matrix_k,
                              this->x_k, this->q_k, this->n_k, this->W);
    // Write the calculated values into the global matrices
    for (uint i=0; i < this->n_bem_points; ++i) {
        for (uint j = 0; j < this->dofs_per_cell; ++j)
        {
            uint J = this->global_bc_dofs[j];
            // loop over quad points only for valid boundary dofs
            if (J != dealii::DoFHandler<dim>::invalid_dof_index)
            {           // beware of the minus! D:
                DL_matrix(i,J) += - local_DL_matrix_k(i,j); //This actually INVERTS the minus assignment in the cuda kernel! X_X
                SL_matrix(i,J) += - local_SL_matrix_k(i,j); // So IN THE END, SL comes with a + and DL with a -!
            }        }    }
}
#endif

#endif // BEMFORMS_HH
