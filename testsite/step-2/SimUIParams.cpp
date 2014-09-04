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


#include <step-2/SimUIParams.h>

// @sect4{Function: declare}
//
// This function declares the parameters
// such that the @p Parameterhandler can try to retrieve them.
void
step2::SimUIParams::declare(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    prm.declare_entry("Run directory", "./test_me",
                      dealii::Patterns::Anything(),
                      "Specify a directory results of "
                      "the test are to be stored. This can be either "
                      "an absolute path or path relative to the directory "
                      "where the program has been started. The default is "
                      "subdir called test_me-<date> where <date> will be replaced "
                      "by the date at which the program has been started. "
                      "this simplifies keeping the projects directory clean "
                      "");

    prm.leave_subsection();

    // Now. declare the parameters of the simulation.
    prm.enter_subsection("Global parameters");

    prm.declare_entry("Run CPU-BLAS vs CUBLAS float", "true",
                      dealii::Patterns::Bool(),
                      "CPU-BLAS and CUBLAS float.");

    prm.declare_entry("Run CPU-BLAS vs CUBLAS double", "true",
                      dealii::Patterns::Bool(),
                      "CPU-BLAS and CUBLAS double.");

    prm.declare_entry("Run CPU-BLAS vs Fujimoto float", "true",
                      dealii::Patterns::Bool(),
                      "CPU-BLAS and Fujimoto float.");

    prm.declare_entry("Run CPU-BLAS vs Fujimoto double", "false",
                      dealii::Patterns::Bool(),
                      "CPU-BLAS and Fujimoto double.");

    prm.declare_entry("Run Fujimoto vs CUBLAS float", "true",
                      dealii::Patterns::Bool(),
                      "Fujimoto vs CUBLAS float.");

    prm.declare_entry("Run Fujimoto vs CUBLAS double", "false",
                      dealii::Patterns::Bool(),
                      "Fujimoto vs CUBLAS double.");

    prm.leave_subsection();


    prm.enter_subsection("Fujimoto parameters");

    prm.declare_entry("Fujimoto variant", "generic",
                      dealii::Patterns::Selection("original|generic|no bitshift|double optimization"),
                      "In case of float the <generic> version is a verbatim copy "
                      "of the original version. In case of double the reading of the matrix "
                      "induces 50% idle threads. <no bitshift> is identical to the generic version "
                      " except that index arithmetic is done in a more human readable fashion. "
                      "<double optimization> provides optimized access to the matrix entries "
                      "in the case of double precision.");

    prm.leave_subsection();
}


// @sect4{Function: get}
//
// This function is responsible for reading the parameters from the prm-file.
void
step2::SimUIParams::get(dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Simulation basics");

    this->run_dir = (prm.get("Run directory").c_str() );

    prm.leave_subsection();

    prm.enter_subsection("Global parameters");

    run_cpublas_vs_cublas_float = prm.get_bool("Run CPU-BLAS vs CUBLAS float");
    run_cpublas_vs_cublas_double = prm.get_bool("Run CPU-BLAS vs CUBLAS double");
    run_cpublas_vs_Fujimoto_float = prm.get_bool("Run CPU-BLAS vs Fujimoto float");
    run_cpublas_vs_Fujimoto_double = prm.get_bool("Run CPU-BLAS vs Fujimoto double");
    run_Fujimoto_vs_cublas_float = prm.get_bool("Run Fujimoto vs CUBLAS float");
    run_Fujimoto_vs_cublas_double = prm.get_bool("Run Fujimoto vs CUBLAS double");

    prm.leave_subsection();


    prm.enter_subsection("Fujimoto parameters");

    std::string v = prm.get("Fujimoto variant");


    if (v == "original")
        this->fj_version = 0;
    else if (v == "generic")
        this->fj_version = 1;
    else if (v == "no bitshift")
        this->fj_version = 2;
    else if ("double optimization")
        this->fj_version = 3;

    prm.leave_subsection();

    // Depending on the selections made in the parameter file
    // the two groups get populated.
    if (run_cpublas_vs_Fujimoto_float || run_Fujimoto_vs_cublas_float){
        float_tests.push_back(Fujimoto_mv);
        float_vs.insert(std::pair<MVCase, int>(Fujimoto_mv, float_tests.size() - 1));
    }

    if (run_cpublas_vs_Fujimoto_double || run_Fujimoto_vs_cublas_double){
        double_tests.push_back(Fujimoto_mv);
        double_vs.insert(std::pair<MVCase, int>(Fujimoto_mv, double_tests.size() - 1));
    }

    if (run_cpublas_vs_cublas_float || run_Fujimoto_vs_cublas_float){
        float_tests.push_back(cublas_mv);
        float_vs.insert(std::pair<MVCase, int>(cublas_mv, float_tests.size() - 1));
    }

    if (run_cpublas_vs_cublas_double || run_Fujimoto_vs_cublas_double){
        double_tests.push_back(cublas_mv);
        double_vs.insert(std::pair<MVCase, int>(cublas_mv, double_tests.size() - 1));
    }

    if (run_cpublas_vs_cublas_float || run_cpublas_vs_Fujimoto_float){
        float_tests.push_back(atlas_mv);
        float_vs.insert(std::pair<MVCase, int>(atlas_mv, float_tests.size() - 1));
    }

    if (run_cpublas_vs_cublas_double || run_cpublas_vs_Fujimoto_double){
        double_tests.push_back(atlas_mv);
        double_vs.insert(std::pair<MVCase, int>(atlas_mv, double_tests.size() - 1));
    }

}
