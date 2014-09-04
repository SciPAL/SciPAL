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


#ifndef MVTESTUIPARAMS_H
#define MVTESTUIPARAMS_H

#include <step-2/SimUIParams.h>
#include <step-2/DeviceParams.h>
#include <step-2/MVMultDriverInterface.h>

namespace step2 {

// @sect3{Class: TestUIParamsBase}
//
// All factorization methods have to know the dimensions of the input
// matrix and the destination where the results should be stored.
// Therefore these details are collected in a separate class
// which allows us to re-use them in other projects.
struct TestUIParamsBase {

    TestUIParamsBase() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    // problem size. A minimal and maximal column number is given.
    int  min_n_cols;
    int  max_n_cols;

    // growth-rate of rows until column number is reached.
    double n_rows_growth_rate;

    // For time measurements it is convenient to repeat the test several times and average over the results.
    int n_repetitions;


private:
    // Dummy implementation to avoid the automatic creation by the compiler.
    // Whenever an object should be copied this results in an error at compile-time.
    TestUIParamsBase(const TestUIParamsBase & ) {}

    TestUIParamsBase & operator = (const TestUIParamsBase & /*other*/)
    {
        return *this;
    }
};


// @sect3{Class: MVTestUIParams}
//
// All run-time parameters are assembled in this structure
// so that they can be passed around as a whole.
// We compose the different subsections by public inheritance.
struct MVTestUIParams
        :
        public SimUIParams,
        public DeviceParams,
        public TestUIParamsBase
{
    MVTestUIParams()
        :
          SimUIParams(),
          DeviceParams(),
          TestUIParamsBase()
    {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);


    int n_random_trials;

    std::vector<step2::matrixSize> matrix_sizes;

protected:
    void create_random_matrix_sizes();
    void create_regular_matrix_sizes();

private:
    MVTestUIParams (const MVTestUIParams & /*other*/)
        :
          SimUIParams(),
          DeviceParams(),
          TestUIParamsBase()
    {}

    MVTestUIParams & operator = (const MVTestUIParams & /*other*/)
    {
        return *this;
    }
};

}
#endif // MVTESTUIPARAMS_H
