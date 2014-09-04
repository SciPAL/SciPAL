/*This file is part of SciPAL.

    SciPAL is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SciPAL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
*/


#ifndef MVTEST_HH
#define MVTEST_HH

#include <step-2/MVTest.h>
#include <step-2/PrecisionName.h>
#include <deal.II/base/convergence_table.h>

// @sect4{Constructor}
//
// @param p : Reference to the full set of runtime parameters.
// @param rt : Reference to the runtimes table.
// @param mv_variant : Selects the implementation of the matrix-vector product.
template <typename T>
step2::MVTest<T>::MVTest(const step2::MVTestUIParams & p,
                         dealii::ConvergenceTable &rt,
                         MVCase mv_variant)
    :
      params(&p),
      results_table(rt)
{
    driver_m = NULL;

    // Depending on the variant we instantiate the driver
    // and choose a suitable header for the column in the results table.
    switch (mv_variant) {

    case Fujimoto_mv:
        driver_m = new FujimotoDriver<T,cublas> (this->params->fj_version);
        col_head = "Fujimoto " + QString::number(this->params->fj_version).toStdString();
        break;
    case cublas_mv:
        driver_m = new CUBlasDriver<T,cublas> ();
        col_head  = "CUBLAS";
        break;
    case atlas_mv:
        driver_m = new step2::CpuBlasDriver<T,blas> ();
        col_head  = "CPU Blas";
        break;
    default: // do nothing
        break;
    }

    if (mv_variant != none)
    {
        assert(driver_m);

        // After checking that we have successfully allocated the driver object
        // we append a string to the column header indicating the precision and number type
        // used in this test.
        col_head += " " + PrecisionName<T>::name();

        // Although it is not an error we use
        std::cerr << "\nTesting " << col_head.c_str() << " mvmult "
                  << std::endl;
        // so that QTCreator highlights this message in red.
    }
}

// @sect4{Destructor}
//
template <typename T>
step2::MVTest<T>::~MVTest()
{
    if (driver_m)
        delete driver_m;
}


// @sect4{Function: run}
//
// This function manages the loop over the matrix sizes. The loops over the repetitions are delegated to the
// different driver classes.
template <typename T>
QString
step2::MVTest<T>::run()
{
    for (size_t i=0; i< this->params->matrix_sizes.size(); i++)
    {
        size_t nr = this->params->matrix_sizes[i].first ;
        size_t nc = this->params->matrix_sizes[i].second;

#ifdef DEBUG
        std::cout  << "Testing MV for " << nr << "x" << nc << " matrix" << std::endl;
#endif
        // For unit testing purposes we reallocate the
        // complete matrix-vector multiplication test for each
        // matrix size to test.
        setup_and_assemble(nr, nc);


        // We copy the reference right-hand side into a local vector
        // and allocate a temporary destination vector.
        std::vector<T> x(this->x_orig.begin(), this->x_orig.end());

        // The reinitialization of @p y is done such that we can use it as matrix in
        // step-3 which tests different matrix-matrix multiplications.
        unsigned int n_elements =  this->y_orig.size();
        std::vector<T> y(n_elements, 0.);

        double elapsed_time
                =
                driver_m->mvmult(y, A, x,
                                 this->params->n_repetitions) / this->params->n_repetitions;

        // After the test we check whether @p y @p == @p y_orig.
        // As error measure we use the max-norm of the difference divided by the reference solution.
        // This yields an error roughly independent of the matrix size.

        dealii::Vector<T> diff ( n_elements );
        for (unsigned int i = 0; i < n_elements; i++)
            diff(i) = (y[i] - y_orig(i)) / y_orig(i);

        double linfty_error = diff.linfty_norm();

        // If the error appears to be too large, we print a warning message.
        // For sufficiently small vectors we dump them onto the screen.
        if (linfty_error> (sizeof(T)<8 ? 1e-5 : 1e-14))
        {
            std::cerr << nr << "x" << nc << " matrix : "
                      << "|| (y - y_orig)/y_orig||_infty = " << linfty_error
                      << " MVTest probably failed!"
                      << std::endl;

            // The maximum length of a vector which can still be reasonably displayed
            // in a terminal is something like 20. Since this is mostly needed when
            // debugging this magic number should not harm anything.
            if (y.size() < 20) {
                std::cerr << y_orig << std::endl;

                // We do not need the @p diff vector anymore, so we use it
                // for storing the result of the matrix-vector product
                // such that it can be put into the out stream.
                std::copy(y.begin(), y.end(), diff.begin());
                std::cerr << diff << std::endl;
            }
        }

        // Finally, add the measured runtimes to the table collecting the results
        // and set the precision of the column's contents.
        this->results_table.add_value(col_head.c_str(), elapsed_time);

        this->results_table.set_precision(col_head.c_str(), 12);
    }

    // We return the header of the results column which simplifies the further processing of the results.
    return col_head.c_str();
}


// @sect4{Function: setup_and_assemble}
//
// The second purpose of the MVTest class is the creation of reasonable test data.
// This is done in this function. This is the place where one could use the deal.II examples
// as blackbox matrix generators.
// For the time being we define the entries of the
// matrix $A \in \mathbb{R}^{m \times n}$,
// the source vector $x\in \mathbb{R}^{n}$ and
// destination $y\in \mathbb{R}^{m}$ as
// \f{eqnarray}
// A_{ik} & = & i + k\,, \\
// x_k & = & \frac{1}{k}\,, \\
// y_i& = & \sum_{k=1}^n A_{ik} x_k  = n + i\sum_{k=1}^n\frac{1}{k} \quad i \in
// \lbrace 1, \ldots , m\rbrace\,.
// \f}
// Having different assembly routines is another possible extension of this program.
template <typename T>
void
step2::MVTest<T>::setup_and_assemble(unsigned int nr, unsigned int nc)
{
    // First we have to setup our linear system, i.e. to allocate enough memory.
    this->A.reinit(nr, nc);
    this->x_orig.reinit(nc);
    this->y_orig.reinit(nr);

    // Then we assemble. It would be nice to have that in parallel as well.
    y_orig = 0.;
    int tmp = 1;
    for (unsigned int r = 0; r < nr; ++r)
        for (unsigned int c = 0; c < nc; ++c)
        {
            x_orig(c)  = (c+1);
            A(r,c)     = tmp; tmp++;//r+1 + 1./(c+1);
            y_orig(r) +=  A(r,c)*x_orig(c);
        }
}
#endif
