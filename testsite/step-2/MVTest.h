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


#ifndef MVTEST_H
#define MVTEST_H

#include <QDir>
#include <step-2/Fujimoto_driver_step-2.h>
#include <step-2/cuda_driver_step-2.h>
#include <step-2/CpuBlas_driver_step-2.h>
#include <step-2/MVTestUIParams.h>

#include <deal.II/base/convergence_table.h>

namespace step2 {

// @sect3{Class: MVTest}
//
// This class drives the tests of the different tuning parameters
// for Fujimoto's matrix-vector product or the test of the CUBLAS
// or ATLAS (CPU) matrix-vector product.
// The sole template parameter of this class is the number type,
// i.e. real float and double
// (or the complex counterparts in a future version).
template <typename T>
class MVTest {

public:
    typedef ::FullMatrixAccessor<T> FullMatrixAccessor;

    MVTest(const MVTestUIParams & p,
           dealii::ConvergenceTable & results_table,
           MVCase variant=atlas_mv);

    ~MVTest();

    virtual QString run();


protected:
    typename step2::MVMultDriverInterface<T> * driver_m;

    virtual void setup_and_assemble(unsigned int nr,
                                    unsigned int nc);

    const MVTestUIParams * params;

    FullMatrixAccessor A;

    dealii::Vector<T> x_orig, y_orig;

    std::string col_head;

    dealii::ConvergenceTable & results_table;
};
}
#endif // MVTEST_H
