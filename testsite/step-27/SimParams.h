// @sect3{File: SimParams.h}
#ifndef SIMPARAMS_H
#define SIMPARAMS_H

#include <QDir>
#include <QTime>


// deal.II:
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>

// For the Architecture enum:
#include <step-27/Architecture.h>
#include <step-27/coupling_type.h>

// Go into the project's namespace again:
namespace step27 {


enum TestExpr { sum, prod, dipole, sum_prod, sum_prod_dipole, drs };

enum SubModel { fem_bem_unit_test, pure_von_Neumann, full_DRS };

// @sect4{Class: SimParams}
//
// Class to set and get the parameters of the project, using the very handy
// deal.II ParameterHandler class.
    struct SimParams {

        // For simplicity we set up the map in the header file, such that everything is in one place.
        SimParams() {

            testcase_str2enum["sum"] = sum;
            testcase_str2enum["prod"] = prod;
            testcase_str2enum["sum_prod"] = sum_prod;
            testcase_str2enum["dipole"] = dipole;
            testcase_str2enum["sum_prod_dipole"] = sum_prod_dipole;
            testcase_str2enum["drs"] = drs;

            sub_model_str2enum["FEM_BEM_unit_test"] = fem_bem_unit_test;
            sub_model_str2enum["pure_von_Neumann"] = pure_von_Neumann;
            sub_model_str2enum["full_DRS"] = full_DRS;
        }


        static void declare(dealii::ParameterHandler & prm);

        void get(dealii::ParameterHandler & prm);

        // Directory where simulation output should be stored.
        QDir run_dir;

        // Directory where the parameter settings of a run should be logged.
        // This will always be @p run_dir + "/log".
        QDir prm_log_dir;

        // Directory and file where the output of the timing should be stored.
        QDir timing_output_dir;



        // The following flag allows to run either one of the convergence tests or
        // the DRS simulation. The latter is going to be effective only in a later step.
        TestExpr fem_bem_test_case;

        SubModel sub_model;

        Architecture arch;

    private:
// Private copy constructor and assignement operator, to prevent
// anything from using them.
        SimParams (const SimParams & /*other*/) {}

        SimParams & operator = (const SimParams & /*other*/) { return *this; }

        // Deal.II's parameter handler returns strings but we need enums.
        std::map<std::string, TestExpr> testcase_str2enum;

        std::map<std::string, SubModel> sub_model_str2enum;
    };

}


#endif // SIMPARAMS_H
