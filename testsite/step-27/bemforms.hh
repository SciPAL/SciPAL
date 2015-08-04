#ifndef BEM_FORMS_S27_HH
#define BEM_FORMS_S27_HH

#include <step-27/bemforms.h>

namespace step27
{
// @sect5{Function: assemble_sl_dl}
//
// This function assembles SL and DL matrix by collocation.
// @param DL_matrix a reference to a FullMatrix of correct dimensions to store the DL entries.
// @param SL_matrix a reference to a FullMatrix of correct dimensions to store the SL entries.
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
    case both:
     //   break;
    default:
        dealii::ExcNotImplemented();
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

    FullMatrixAccessor<double> local_DL_matrix_k(0,0, this->is_col_major);
    FullMatrixAccessor<double> local_SL_matrix_k(0,0, this->is_col_major);


    // ...then we need to resize the containers
    local_DL_matrix_k.reinit(this->n_bem_points, this->dofs_per_cell);
    local_SL_matrix_k.reinit(this->n_bem_points, this->dofs_per_cell);

    const std::vector<dealii::Point<dim> > &
            x = support_points.begin()->second;

    this->point_to_double_vector(this->x_k, x);

    this->point_to_double_vector(this->q_k, this->q);
    this->point_to_double_vector(this->n_k, this->normals);

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
