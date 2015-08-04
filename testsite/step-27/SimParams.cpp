// @sect3{File: SimParams.cpp}

#include <step-27/SimParams.h>

#include <deal.II/base/exceptions.h>

#include <QDebug>
// @sect4{Class: SimParams}
// @sect5{Function: declare}
//
// Declare the general parameters of the project. Currently, this is only the directory for the results.
// In later steps, where we add the CUDA parallelization, there will be additional parameters,
// for instance for selecting the type of parallelization.
void
step27::SimParams::declare(dealii::ParameterHandler & prm)
{
    prm.declare_entry("Run directory",
                      "./test_me-",
                      dealii::Patterns::Anything(),
                      "Specify a directory where the results of "
                      "the simulation should be stored. This can be either "
                      "an absolute path or path relative to the directory "
                      "where the program has been started. The default is a"
                      "subdir called test_me-<date> where <date> will be replaced "
                      "by the date at which the program has been started. "
                      "This simplifies keeping the projects directory clean.");

    prm.declare_entry("Use time tag",
                      "true",
                      dealii::Patterns::Bool(),
                      "If set to true, the results of the run will be stored in a subdirectory "
                      "of a given run directory which named after the date and time at the "
                      "start of the program. This allows to group runs without "
                      "overwriting previous results");

    prm.declare_entry("FEM-BEM reference solution", "sum_prod",
                      dealii::Patterns::Selection("sum|prod|sum_prod|dipole|sum_prod_dipole"),
                      "These test cases can be used for convergence tests. Currently available are:\n"
                      "sum: f(x,y,z) = 2*x+y+z\n "
                      "prod: f(x,y,z) = 0.01*x*y*z\n "
                      "sum_prod: f(x,y,z) = sum + prod\n "
                      "dipole: f(x,y,z) = a/|r-r0| - a/|r-r1|, "
                      "where a is the strength of the dipole and r0,r1 are points on the z axis, "
                      "symmetric w.r.t. to the origin. "
                      "a and the lengths of r0, r1 can be tuned in protein_dipole_data.prm\n "
                      "sum_prod_dipole: f(x,y,z) = sum + prod + dipole\n "
                      "drs: do not use a test case, but actually run the drs simulation. Not available in step-27."
                      ".");

    prm.declare_entry("Sub model", "FEM_BEM_unit_test",
                      dealii::Patterns::Selection("FEM_BEM_unit_test|pure_von_Neumann|full_DRS"),
                      "To validate that the full DRS simulation (choose full_DRS) works it is useful to reduce "
                      "the DRS problem to the Poisson equation for the potential which either has to fulfill "
                      "the FEM-BEM problem at the cavity interface and "
                      "Dirichlet boundary conditions elsewhere (FEM_BEM_unit_test)  "
                      "or a Poisson equation which is only subject to "
                      "von Neumann boundary data (pure_von_Neumann). In the latter "
                      "case the solution is not unique anymore. To fix that, the problem is made pseudo time-dependent "
                      "by adding a time derivative which makes the global constant visible again. "
                      "This parameter allows to test the different aspects of the overall simulation independently.");

    prm.declare_entry("Architecture", "cpu",
                      dealii::Patterns::Selection("cpu|cuda|both"),
                      "Where to run the assembly of the BEM part. \n"
                      "cpu: do it on the CPU\n"
                      "cuda: use your GPU\n"
                      "both: use both (for speedup tests, nor for production runs).");
}


// @sect5{Function: get}
//
// The get function is just a mirror of the declare function.
void
step27::SimParams::get(dealii::ParameterHandler & prm)
{
    std::string tmp_run_dir = prm.get("Run directory");

    std::string testcase = prm.get("FEM-BEM reference solution");

      std::map<std::string, TestExpr>::const_iterator selection = testcase_str2enum.find(testcase);
      AssertThrow(selection != testcase_str2enum.end(),
             dealii::ExcMessage("Illegal test case selected. Check your parameter file!"));

    fem_bem_test_case = selection->second;


    std::string sm = prm.get("Sub model");

      std::map<std::string, SubModel>::const_iterator sm_selection = sub_model_str2enum.find(sm);
      AssertThrow(sm_selection != sub_model_str2enum.end(),
             dealii::ExcMessage("Illegal submodel selected. Check your parameter file!"));

    sub_model = sm_selection->second;


    bool use_time_tag = prm.get_bool("Use time tag");

    // Depending on the user's input the path is extended.
    // We can safely use a slash as QT will convert it into the separator appropriate for the operating system.
    if (use_time_tag)
        tmp_run_dir += (
                QString("/"
                        +
                        QDateTime::currentDateTime().toString("ddd-yyyy-MM-dd/hh_mm_ss")
                        ).remove(".")).toStdString();

    // Get the run directory's absolute path and create it, if it does not exist.
    this->run_dir.setPath(tmp_run_dir.c_str() );
    this->run_dir.makeAbsolute();

    if(!this->run_dir.exists())
        this->run_dir.mkpath(".");

    // Do the same for the log directory such that it is a subdirectory of the run directory.
    this->prm_log_dir.setPath(this->run_dir.absolutePath() + QDir::separator() + "log");
    this->prm_log_dir.makeAbsolute();

    if(!this->prm_log_dir.exists())
        this->prm_log_dir.mkpath(".");


    // Set up the subdirectory for the output of the timing results.
    this->timing_output_dir.setPath(this->run_dir.absolutePath()
                                  + QDir::separator()
                                  + "timing");

    if(!this->timing_output_dir.exists())
        this->timing_output_dir.mkpath(".");

    // Keep for later:
#ifndef nUSE_MULTIARCH
    // get the Architecture:
    std::string arch_string = prm.get("Architecture").c_str();

    std::cout << "chosen architecture : " << arch_string.c_str() << std::endl;

    if(arch_string == "cuda"){
        this->arch = step27::cuda;
    } else if (arch_string == "both") {
        this->arch = step27::both;
    } else if (arch_string == "cpu") {
        this->arch = step27::cpu;
    } else {
        // Throw a fit if an invalid architecture has been provided
#ifdef QT_NO_DEBUG
        AssertThrow(false,
                    dealii::StandardExceptions::ExcMessage("Unknown Architecture string. Must be \"cuda\", \"cpu\", or \"both\".") );
#else
        Assert(false,
               dealii::StandardExceptions::ExcMessage("Unknown Architecture string. Must be \"cuda\", \"cpu\", or \"both\".") );
#endif
    }
#endif

}
