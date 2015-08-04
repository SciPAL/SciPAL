#ifndef SIMULATION_HANDLER_H
#define SIMULATION_HANDLER_H

// SciPAL includes
#include <base/GPUInfo.h>

// The parameter class for your simulation.
#include <step-27/SimParams.h>

namespace step27 {

// @sect4{Class: SimulationManager}
//
// The class used to encapsulate the simulation.
// It takes care of reading the parameter files, setting up the
// Step16::FEMBEMPoissonProblem used for the actual dielectric relaxation
// spectroscopy simulation.

template<typename DRSType, typename SimParamsType>
class SimulationManager {

public:
    SimulationManager(int argc, char *argv[]);

    int run();

private:

    // DRS is always a 3D problem.
    static const int dim = 3;

    std::string prm_drs_numerics_filepath;
    std::string prm_drs_cell_filepath;

    step27::NumericsParams<dim> numerics_prms;

    QFile timing_output_file;

protected:
    SimParamsType params;

    // A simulation has to know where to find its parameters.
    QString  prm_path;
};
}

// @sect5{Constructor: SimulationHandler}
//
// The constructor is responsible for reading parameters
// and initializing the device, i.e. the selected graphics card.
// @param argc : The number of command line arguments. This is always $\ge 1$,
// as by default the zeroth argument is the name of program itself.
// @param argv : Pointer to the array of command line arguments.
// @param g : Reference to the object containing the GPU info from the system.
template<typename DRSType, typename SimParamsType>
step27::SimulationManager<DRSType, SimParamsType>::SimulationManager(int argc,
                                             char *argv[])
{
    // Declare and read parameters from a file. Basically, the parameter
    // file must have the same name as the binary. The extension has
    // to be ".prm". What has been read will be dumped into a log file.
    dealii::ParameterHandler prm_handler;

    SimParamsType::declare(prm_handler);

    // Get the current working directory ...
    QDir cwd = QDir::current();

    // i.e. where the program has been started.
    QDir launch_dir = cwd;

    std::cout << "CWD : " << launch_dir.absolutePath().toStdString().c_str() << std::endl;
    
    // By default, the master parameter file has the same name as the binary
    // and is supposed to be in a subdirectory prm of that directory,
    // where the program has been started.
    // As it may contain paths to the run directory and where results are to be stored
    // it has to be processed before all other parameter files.
    std::string master_prm_filename;
    if (argc == 1)
    {
        this->prm_path = (launch_dir.absolutePath().toStdString() + "/prm/").c_str();

        QFileInfo tmp(argv[0]);
        master_prm_filename  = tmp.baseName().toStdString();

        this->prm_drs_cell_filepath     = (this->prm_path + master_prm_filename.c_str() + "_drs_cell.prm").toStdString();
        this->prm_drs_numerics_filepath = (this->prm_path + master_prm_filename.c_str() + "_drs_numerics.prm").toStdString();

        master_prm_filename = this->prm_path.toStdString() + master_prm_filename + ".prm";

        cwd.setPath("./prm");
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

            qFatal("Cannot proceed without proper path to paramter file");
        }

        // Next, we subdivide the given filename into its path and filename
        // so that the corresponding subdirectories can be created.
        this->prm_path = tmp.absolutePath();
        cwd.setPath(prm_path);
        cwd.makeAbsolute();
        master_prm_filename = tmp.fileName().toStdString();
        this->prm_drs_cell_filepath = tmp.absolutePath().toStdString() + "/"
                + tmp.baseName().toStdString()
                + "_drs_cell."
                + tmp.completeSuffix().toStdString();
        this->prm_drs_numerics_filepath = tmp.absolutePath().toStdString() + "/"
                + tmp.baseName().toStdString()
                + "_drs_numerics."
                + tmp.completeSuffix().toStdString();

        std::cout << "Parameter file path : "
                  << tmp.absolutePath().toStdString().c_str()
                  << std::endl;
    }

    std::cout << "Parameter file : " << master_prm_filename  << std::endl;

    std::cout << "DRSCell Parameter file  : " << this->prm_drs_cell_filepath  << std::endl;

    std::cout << "Numerics Parameter file : " << this->prm_drs_numerics_filepath  << std::endl;

    // Before the parameter file can be read, we have to make sure that
    // its directory exists. In case of the default parameter file
    // the directory will be created.
    if (!cwd.exists() )
        launch_dir.mkpath( cwd.absolutePath() );

    QDir::setCurrent(cwd.absolutePath());

    prm_handler.read_input (master_prm_filename);

    QDir::setCurrent(launch_dir.absolutePath());

    this->params.get(prm_handler);

    // Create toplevel run directory.
    cwd.setPath(this->params.run_dir.absolutePath());

    // The following lets a directory make its own path.
    if (!cwd.exists())
        cwd.mkpath( "." );


    // After the run directory we create the log directory.
    if (!this->params.prm_log_dir.exists())
        this->params.prm_log_dir.mkpath(".");

    // Now, change to the run directory
    QDir::setCurrent(cwd.absolutePath());


    // and write what has been actually read
    // into log file. Basically, this is just another parameter file
    // and thus could be used again as input to another run after stripping the .log suffix.
    std::ofstream log_master_prm( (this->params.prm_log_dir.absolutePath() + QDir::separator()
                                   + QString(master_prm_filename.c_str()) + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_master_prm,
                                  dealii::ParameterHandler::Text);

    // At this point the toplevel run dir must exist.
    // Thus, we can change to it without any further sanity test.
    QDir::setCurrent(this->params.run_dir.absolutePath());

    // Now we set the timing output file
    QFileInfo tmp(master_prm_filename.c_str());

    this->timing_output_file.setFileName( params.timing_output_dir.absolutePath()
                                        + QDir::separator()
                                        + tmp.baseName() + "_timing.dat");

    ////////////////////////////////////
    // dealii::ParameterHandler prm_handler;

    step27::NumericsParams<dim>::declare(prm_handler);
    prm_handler.read_input(prm_drs_numerics_filepath);

    numerics_prms.get(prm_handler);

    numerics_prms.timing_file = this->timing_output_file.fileName().toStdString();

    numerics_prms.prm_drs_cell_filepath = this->prm_drs_cell_filepath;

    numerics_prms.arch = this->params.arch;

    std::ofstream log_numerics_prm( (this->params.prm_log_dir.absolutePath() + QDir::separator() + QFileInfo(prm_drs_numerics_filepath.c_str()).fileName() + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_numerics_prm,
                                  dealii::ParameterHandler::Text);

}

// @sect5{Function: run}
//
// The run function. This function sets up the FEMBEMPoissonProblem and executes it.
// At the current stage it is not much more than a delegator for calling the FEMBEMPoissonProblem's @p run()
// function. It smain purpose is to handle the reading of the parameters needed for setting up the FEM BEM problem.
template<typename DRSType, typename SimParamsType>
int step27::SimulationManager<DRSType, SimParamsType>::run()
{
    /*
    dealii::ParameterHandler prm_handler;

    step27::NumericsParams<dim>::declare(prm_handler);
    prm_handler.read_input(prm_drs_numerics_filepath);

    numerics_prms.get(prm_handler);

    numerics_prms.timing_file = this->timing_output_file.fileName().toStdString();

    numerics_prms.prm_drs_cell_filepath = this->prm_drs_cell_filepath;

    numerics_prms.arch = this->params.arch;
*/
    try {
        dealii::deallog.depth_console (0);
        // Set up the @p Step16::FEMBEMPoissonProblem ...
        DRSType laplace_problem(numerics_prms,
                                this->prm_path.toStdString(),
                                this->params.prm_log_dir,
                                this->params.fem_bem_test_case,
                                this->params.sub_model);
        // ...and run it
        laplace_problem.run ();

        // and catch execptions
    } catch (std::exception &exc) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;

        return 1;
    } catch (...) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    return 0;
}
#endif // SIMULATION_HANDLER_H
