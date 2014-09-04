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


#ifndef STEP_2_HH
#define STEP_2_HH
#include <cstdio>
// We walk through the program top-down starting with the file which contains the main() function.
//
// We need several headers. First of all
// STL's vector class.
#include <vector>

// To store the measured runtimes in a table we use
// deal.II's
#include <deal.II/base/convergence_table.h>

// The basis for all CUDA-related programs is provided by the CUDA runtime library.
#include <cuda_runtime.h>
// Our SciPAL library provides a little helper class which offers a more comfortable interface
// to the technical data about the available GPUs.
#include <../SciPAL/include/base/GPUInfo.h>

// Finally, one of the greatest aids in all the less-scientific issues in a program is QT's string class
#include <QString>

// Last but not least, we need
// the headers for the matrix-vector benchmark test
#include <step-2/MVTest.h>
#include <step-2/MVTest.hh>

// and for the user interface.
#include <step-2/MVTestUIParams.h>

// We encapsulate the whole project into a dedicated namespace
// in order to be able to re-use
// parts of this program in others.
namespace step2 {

// @sect3{Class: MyFancySimulation}
//
// To make this test facility extendible, we implement
// a class for a simple user interface. Its primary tasks are:
// - Management of run-time parameters by a simple text-based parameter file.
// - Setting device parameters according to the user's parameters.
// - Preprocessing and output of results.
// The result of this program is a comparison of the performance of different
// implementations of the matrix-vector product for dense matrices and vectors.
// The runtimes are collected in a central table @p results_table. It is an instance
// of the dealii::ConvergenceTable class which has such nice member functions
// like @p write_tex() which writes the contents to a tex file such that it can be
// directly included in the tex source of an article.
// <br />
// The non-public part of this class is declared only as protected and not as private
// in order to be able to reuse this class in other projects doing similar benchmark tests.
// This is also the reason why the two main functions of this class are declared as <i>virtual</i>.
class MyFancySimulation {

public:
   MyFancySimulation(int argc, char *argv[]);

   virtual void run();

protected:
    virtual void save_results();


    dealii::ConvergenceTable results_table;

    MVTestUIParams params;

    SciPAL::GPUInfo gpu_info;

    QString launch_dir;

    QString prm_dir;

    QString prm_log_dir;
};

}

// @sect4{Constructor: MyFancySimulation}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
step2::MyFancySimulation::MyFancySimulation(int argc, char *argv[])
{
    // Before setting up the simulation we
    // figure out how many GPUs are available.
    cudaGetDeviceCount(&gpu_info.n_CUDA_devices);
    std::cout
            << "N available CUDA devices : "
            << gpu_info.n_CUDA_devices << std::endl;

    gpu_info.get();

    // Then, we declare and read parameters from a file. Basically, the parameter
    // file must have the same name as the binary. The extension has
    // to be ".prm". What has been read will be dumped into a log file.
    // We could use the parameter file to read a device ID and once the parameters
    // are read we could reset the @p gpu_info object to the new device.
    dealii::ParameterHandler prm_handler;

    MVTestUIParams::declare(prm_handler);

    // Before reading the parameters we set up the directories needed
    // during the execution of the program.
    // As point of reference we use directory where the program has been started.
     QDir cwd = QDir::current();
     cwd.makeAbsolute();

    // This is also the current working directory.
   this->launch_dir = cwd.absolutePath();

    // By default, the parameter file has the same name as the binary
    // and is supposed to be in a subdirectory prm of that directory
    // where the program has been started.
    std::string prm_filename;

    if (argc == 1)
    {
        QFileInfo tmp(argv[0]);
        this->prm_dir = tmp.absolutePath() + "/prm";
        prm_filename  = tmp.fileName().toStdString();
        prm_filename += ".prm";
    }
    else
    {
        // Whatever gets passed as first command line argument is considered as path
        // to a parameter file.
        std::cout << "Given parameter file : " << argv[1] << std::endl;

        // We convert the sequence of characters into something more meaningful.
        QFileInfo tmp(argv[1]);

        // Before we proceed, let us figure out whether the given parameter file exists.
        // Note: If the file is a symlink that points to a non existing file,
        // false is returned as well.
        if(!tmp.exists())
        {
            std::cerr << "The following parameter file does not exist:\n"
                      << argv[1] << std::endl;

            qFatal("Cannot proceed without proper path to parameter file");
        }

        // Next, we subdivide the given filename into its path and name
        // so that the corresponding subdirectories can be created.
        this->prm_dir = tmp.absolutePath();

        prm_filename = tmp.fileName().toStdString();

        std::cout << "xx Parameter file path : "
                  << tmp.absolutePath().toStdString().c_str()
                  << std::endl;
    }

    // For control purposes it is useful to notify the user of which parameter file has been actually used.
    std::cout << "Parameter file : " << prm_filename  << std::endl;

    // Before the parameter file can be read, we have to make sure that
    // its directory exists. In case of the default parameter file
    // the directory will be created.
    cwd.mkpath(this->prm_dir);

    std::cout << "step-2: prm path : " << this->prm_dir.toStdString().c_str()  << std::endl;

    QDir::setCurrent(this->prm_dir);

    prm_handler.read_input (prm_filename);

    this->params.get(prm_handler);

    // Having read the parameters we can create the toplevel run directory which was given in the parameter file.
    QDir::setCurrent(this->launch_dir);
    cwd.mkpath(this->params.run_dir);
    QDir::setCurrent(this->params.run_dir);
    cwd = QDir::current();
    this->params.run_dir = cwd.absolutePath();

    std::cout << "Entering run dir : " << this->params.run_dir.toStdString().c_str()  << std::endl;



    // Create the log directory and write what has been actually read
    // into a log file. Basically, this is just another parameter file
    // and thus can be used again as input to another run after stripping the .log suffix.
    this->prm_log_dir = this->params.run_dir + "/log";
    cwd.mkpath(this->prm_log_dir);



    QDir::setCurrent(this->prm_log_dir);

    std::ofstream log_out_text( (prm_filename +".log" ).c_str() );
    prm_handler.print_parameters (log_out_text,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    assert(QDir::setCurrent(this->params.run_dir ) );

    // We know all the matrix sizes we are going to test and therefore
    // initialize @p results_table with the number of rows and columns of the test matrices.
    for(unsigned int i = 0; i < this->params.matrix_sizes.size(); i++)
    {
        this->results_table.add_value("rows",  (int) this->params.matrix_sizes[i].first);
        this->results_table.add_value("columns", (int) this->params.matrix_sizes[i].second);
    }
}

// @sect4{Function: run}
//
// Running the performance test is done here.
void step2::MyFancySimulation::run()
{
    // The test runs are subdivided into groups for
    // single and double precision.
    // The list of matrix-vector product implementations
    // which are to be executed are stored as vector.
    // To keep the code short we introduce some references.
    const std::vector<MVCase> & float_tests  = this->params.float_tests;
    const std::vector<MVCase> & double_tests = this->params.double_tests;

    // For the comparisons of the runtimes we need the reverse, i.e. a mapping
    // from the name of the implementation to its index in the list of tests.
    const std::map<MVCase, int> & float_vs  = this->params.float_vs;
    const std::map<MVCase, int> & double_vs = this->params.double_vs;

    // In the benchmark tests we measure the runtimes for different matrix sizes.
    // The plot command for this raw data is assembled in
    QString gpplots_runtimes;

    // The interesting thing is which matrix-vector product is faster. To quantify this
    // we will compute the speedup of the CUDA implementations over ATLAS.
    // Basically, gnuplot can do this for us. The necessary plot commands are
    // stored in a separate string.
    QString gpplots_speedup;

    // The string @p gnuplot will contain all commands we want to pass to gnuplot
    // and merges the plot style settings with the plot commands.
    // For a full documentation of gnuplot's capabilities
    // have look at its <a href="http://www.gnuplot.info">online doc</a>.
    QString gnuplot
            =
            //We use the postscript terminal for publication-ready results.
            "set term postscript landscape enhanced color solid "
            " linewidth 2.0 \"Helvetica\" 20\n"
            "set xlabel \"matrix entries\"\n"
            "set ylabel \"execution time\"\n"
            "set logscale xy\n"
            "set grid\n"
            "set output \"runtimes.ps\"\n"
            "set key inside left Left box lw 0.5\n";

    gnuplot += "plot ";

    // Before we can complete the plot command we have to run the tests.
    // Float comes first.
    std::vector<MVCase>::const_iterator
            t = float_tests.begin(),
            end_t = float_tests.end();

    int col = 3;

    for (; t != end_t; ++t)
    {
        MVTest<float> mv_test(this->params, results_table, *t);
        if(col > 3) gpplots_runtimes += ",";
        // @p gpplots_runtimes contains all plot commands for the pure runtimes.
        // The relative path to the data file reflects the fact that we will call gnuplot
        // from a subdirectory of the run directory.
        gpplots_runtimes += QString("\"../MVTest_results.out\" using ($1*$2):")
                // Each run of @p mv_test appends an additional column with results to <i>MVTest_results.out</i>.
                // By stepping through the columns we can plot the respective runtime
                + QString::number(col++)
                + QString(" title \"") + mv_test.run() + QString("\" w p");
    }

    // Then come the double tests.
    t = double_tests.begin(), end_t = double_tests.end();

    for (; t != end_t; ++t)
    {
        MVTest<double> mv_test(this->params, results_table, *t);
        if(col > 3) gpplots_runtimes += ",";
        // After the tests have finished we can compose the plotting command for the double results
        gpplots_runtimes += QString("\"../MVTest_results.out\" using ($1*$2):")
                +
                QString::number(col++)
                +
                QString(" title \"") + mv_test.run() + QString("\" w p");
    }

    // Add the runtime plot commands to the global command string.
    gnuplot += gpplots_runtimes;

    // Then, set up titles for the speedup graphs.
    gnuplot +=
            "\nset ylabel \"speedup\"\n"
            "set output \"speedup.ps\"\n"
            "set key inside left Left\n";

    // For the speedups a linear scaling of the y-axis is more suitable
    gnuplot += "unset logscale y\n";

    // The order of the execution of the tests determines the order
    // in which the runtimes appear in @p results_table.
    // To compose the command string for plotting the speedups we use our lookup tables @p float_vs and @p double_vs
    // to find out which columns we have to use for
    // the calculation of the speedup.
    // The column index is associated with the kind of test (like @p atlas_mv). Using the test as key
    // the map returns the respective column to which we still have to add
    // a 3 since the map starts with index 0 but the first 2 columns
    // in our outputfile are reserved for the rows and columns of the matrix.
    gnuplot += "plot";

    int offset = 3 + this->params.float_vs.size();

    if(this->params.run_cpublas_vs_cublas_float)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(float_vs.at(atlas_mv) + 3).arg(float_vs.at(cublas_mv) + 3)
                + QString(" title \"CPU BLAS vs CUBLAS (float)\" w p");

    if(this->params.run_cpublas_vs_cublas_double)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(double_vs.at(atlas_mv) + offset).arg(double_vs.at(cublas_mv) + offset)
                + QString(" title \"CPU BLAS vs CUBLAS (double)\" w p");

    if(this->params.run_cpublas_vs_Fujimoto_float)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(float_vs.at(atlas_mv) + 3).arg(float_vs.at(Fujimoto_mv) + 3)
                + QString(" title \"CPU BLAS vs Fujimoto (float)\" w p");

    if(this->params.run_cpublas_vs_Fujimoto_double)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(double_vs.at(atlas_mv) + offset).arg(double_vs.at(Fujimoto_mv) + offset)
                + QString(" title \"CPU BLAS vs Fujimoto (double)\" w p");

    if(this->params.run_Fujimoto_vs_cublas_float)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(float_vs.at(Fujimoto_mv) + 3).arg(float_vs.at(cublas_mv) + 3)
                + QString(" title \"Fujimoto vs CUBLAS (float)\" w p");

    if(this->params.run_Fujimoto_vs_cublas_double)
        gpplots_speedup += ((gpplots_speedup.isEmpty()) ? QString("") : QString(", "))
                + QString("\"../MVTest_results.out\" using ($1*$2):($%1 / $%2)")
                .arg(double_vs.at(Fujimoto_mv) + offset).arg(double_vs.at(cublas_mv) + offset)
                + QString(" title \"Fujimoto vs CUBLAS (double)\" w p");

    // When the plotting commands for the speedup are complete we can add them to the global command string.
    gnuplot += gpplots_speedup;

    // Besides plotting data gnuplot can also be employed as interface to the console. We use this to
    // convert the images from postscript to pdf and to remove old ps-files.
    gnuplot += "\n!ps2pdf runtimes.ps runtimes.pdf";
    gnuplot += "\n!rm runtimes.ps";
    gnuplot += "\n!ps2pdf speedup.ps speedup.pdf";
    gnuplot += "\n!rm speedup.ps";


    std::cout << "Results are processed and saved.\n";

    // Finally, save the runtimes to disc.
    this->save_results();

    // Then, we create the subdirectory for the plots from which gnuplot is going to be called.
     QDir cwd = QDir::current();
     QString plot_dir = this->params.run_dir + "/plot";
             cwd.mkpath("./plot");

     QDir::setCurrent(plot_dir);


    //We write the plotting commands from the string @p gnuplot to a file <i>plot.gp</i>,
    QFile plotscript("plot.gp");

    bool success = plotscript.open(QIODevice::WriteOnly);
    if (!success)
        std::cerr << "Opening gnuplot file failed!" << std::endl;

    plotscript.write(gnuplot.toStdString().c_str());

    plotscript.close();

    if (! plotscript.exists() )
        std::cerr << "Writing gnuplot file failed!" << std::endl;


    // open a pipe to gnuplot and load the plot file we have just created.
    FILE *gp = popen("gnuplot -persist", "w");
    fprintf(gp, "load \"plot.gp\"\n");
    fflush(gp);
    // After closing the pipe to gnuplot everything should have been plotted.
    pclose(gp);


    QDir::setCurrent(this->params.run_dir);

    std::cout << "Finished." << std::endl;
}

// @sect4{Function: save_results}
//
// The final step in the benchmarking is to
// save the content of @p results_table to a text file.
// It should end up in the run directory given in the parameter file.
// The header of the results file contains the information
// about the sampling of the matrix sizes which allows to recreate
// the input data at a later time solely from the results file.
void step2::MyFancySimulation::save_results()
{
    std::string filename("MVTest_results");

    filename += ".out";

    std::ofstream out(filename.c_str());

    out << "# max n_cols  : " << this->params.max_n_cols << std::endl
        << "# min n_cols  : " << this->params.min_n_cols << std::endl
        << "# growth rate : " << this->params.n_rows_growth_rate << std::endl
        << std::endl;

    results_table.write_text(out);
}

#endif
