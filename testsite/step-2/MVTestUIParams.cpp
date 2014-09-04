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


#include <step-2/MVTestUIParams.h>

// @sect4{Function: declare}
//
// Declare the accepted parameters in the parameter-handler.
void
step2::TestUIParamsBase::declare(dealii::ParameterHandler & prm)
{

    // This parameterhandler-class provides the functionality to divide the parameter-file into sections
    // and to define allowed intervals for parameters with numerical values.

    prm.enter_subsection("Dimensions of test problems.");

    prm.declare_entry ("log2(min n cols)",
                       "4",
                       dealii::Patterns::Double(2,15),
                       "Binary logarithm of minimal number of columns of "
                       "upper triangular matrix R. "
                       "Allowed range : [2,15]");

    prm.declare_entry ("log2(max n cols)",
                       "12",
                       dealii::Patterns::Double(2,16.01),
                       "Binary logarithm of minimal number of columns of "
                       "upper triangular matrix R. "
                       "Allowed range : [3,16]");

    prm.declare_entry ("n rows growth rate",
                       "1.4",
                       dealii::Patterns::Double(1.1,20),
                       "In each test instance the number of rows of R "
                       "is increased by this factor. "
                       "Allowed range : [1.1,10]");

    prm.declare_entry("n repetitions",
                      "1",
                      dealii::Patterns::Integer(1),
                      "Repeat the test for a given matrix size this many times.");

    prm.leave_subsection();
}



// @sect4{Function: get}
//
// Read parameters into the object
void
step2::TestUIParamsBase::get(dealii::ParameterHandler & prm)
{

    prm.enter_subsection("Dimensions of test problems.");

    min_n_cols = pow(2., prm.get_double ("log2(min n cols)") );

    max_n_cols = pow(2., prm.get_double ("log2(max n cols)") );

    n_rows_growth_rate = prm.get_double("n rows growth rate");

    n_repetitions = prm.get_integer("n repetitions");

    prm.leave_subsection();
}


// @sect4{Function: declare}
//
void
step2::MVTestUIParams::declare(dealii::ParameterHandler &prm)
{
    // First, declare the general simulation parameters, i.e.
    // to specify the different precisions for which to run
    // the tests and working directory. The working directory specifies
    // the toplevel directory in which the program should
    // run and store its results independent of where it has been started.
    SimUIParams::declare(prm);
    DeviceParams::declare(prm);
    TestUIParamsBase::declare(prm);

    prm.enter_subsection("Dimensions of test problems.");

    prm.declare_entry("n random trials",
                      "0",
                      dealii::Patterns::Integer(0),
                      "Perform the test for this many randomly chosen matrix sizes."
                      "If the value given is 0, then a list of matrix sizes is used "
                      "which is generated deterministically from the growth strategy.");
    prm.leave_subsection();
}


// @sect4{Function: get}
//
void
step2::MVTestUIParams::get(dealii::ParameterHandler &prm)
{
    this->SimUIParams::get(prm);

    this->DeviceParams::get(prm);

    this->TestUIParamsBase::get(prm);

    prm.enter_subsection("Dimensions of test problems.");

    n_random_trials = prm.get_integer("n random trials");

    prm.leave_subsection();

    if (n_random_trials>0)
        create_random_matrix_sizes();
    else
        create_regular_matrix_sizes();
}


// @sect4{Function: create_random_matrix_sizes}
//
void
step2::MVTestUIParams::create_random_matrix_sizes()
{
    std::cerr << __FUNCTION__ << std::endl;

    srand ( time(NULL) );

    int maxMatrixDimension = this->max_n_cols;

    int random_trials = this->n_random_trials;

    this->matrix_sizes.resize(random_trials);

    for (size_t i=0; i < random_trials; i++)
    {
        size_t nrows= rand() % maxMatrixDimension + 1;
        size_t ncols= rand() % maxMatrixDimension + 1;

        int mod = 16;
        nrows=nrows + ( mod - nrows % mod );
        ncols=ncols + ( mod - ncols % mod );

        this->matrix_sizes[i] = matrixSize(nrows, ncols);
    }
}


// @sect4{Function: create_regular_matrix_sizes}
//
void step2::MVTestUIParams::create_regular_matrix_sizes()
{
    std::cerr << __FUNCTION__ << std::endl;

    float max_nc = this->max_n_cols;


    for (int nc = this->min_n_cols; nc < max_nc;
         nc=std::max(nc+1, int(nc*this->n_rows_growth_rate) )
         )
//        for (int nr = nc; nr <= max_nc;
//             nr=std::max(nr+1, int(nr*this->n_rows_growth_rate) )
//             )

        { int nr = max_nc;
            this->matrix_sizes.push_back( step2::matrixSize(nr,nc));
            // For rectangular matrices the performance of
            // the CPUBlas changes when @p nr and @p nc get swapped.
            // Therefore, we add both possibilities to the list of our
            // test sizes.
            // if (nc != nr)
             //   this->matrix_sizes.push_back( step2::matrixSize(nc,nr));
        }
}

