/**
 * @page step_42                 The step-42 tutorial program
@htmlonly
<table class="tutorial" width="100%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
      </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ClassGPUInfo">Class: GPUInfo</a>
      <ul>
        <li><a href="#Functionget">Function: get</a>
      </ul>
        <li><a href="#ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#ConstructorMyFancySimulation">Constructor: MyFancySimulation</a>
        <li><a href="#Functionrun">Function: run</a>
      </ul>
        <li><a href="#Functionmain">Function: main</a>
        <li><a href="#ClassCUDADriver">Class: CUDADriver</a>
      <ul>
        <li><a href="#ConstructorCUDADriver">Constructor: CUDADriver</a>
        <li><a href="#Functiongemm_tests">Function: gemm_tests</a>
        <li><a href="#Functiongemv_tests">Function: gemv_tests</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
      </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
    <ul>
        <li><a href="#plain-ClassGPUInfo">Class: GPUInfo</a>
      <ul>
        <li><a href="#plain-Functionget">Function: get</a>
      </ul>
        <li><a href="#plain-ClassMyFancySimulation">Class: MyFancySimulation</a>
      <ul>
        <li><a href="#plain-ConstructorMyFancySimulation">Constructor: MyFancySimulation</a>
        <li><a href="#plain-Functionrun">Function: run</a>
      </ul>
        <li><a href="#plain-Functionmain">Function: main</a>
        <li><a href="#plain-ClassCUDADriver">Class: CUDADriver</a>
      <ul>
        <li><a href="#plain-ConstructorCUDADriver">Constructor: CUDADriver</a>
        <li><a href="#plain-Functiongemm_tests">Function: gemm_tests</a>
        <li><a href="#plain-Functiongemv_tests">Function: gemv_tests</a>
      </ul>
      </ul>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Introduction"></a><h1>Introduction</h1>

This tutorial program illustrates the capabilities of SciPal's DSEL.
We show the evaluation of differnt expressions corresponding to the different (cu)BLAS levels.
The tests are implemented in cuda_driver_step-42.hh.
The perallelization strategy can be choosen by changing the template parameter BW in cuda_driver_step-42.h.  

 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * @endcode
 * 
 * This is standard C++.
 * 
 * @code
 * #include <iostream>
 * #include <vector>
 * 
 * @endcode
 * 
 * Driver class for GPU part of th program.
 * 
 * @code
 * #include <cuda_driver_step-42.h>
 * #include <cuda_driver_step-42.hh>
 * 
 * @endcode
 * 
 * The parameter class for your simulation.
 * 
 * @code
 * #include <SimParams.h>
 * 
 * 
 * @endcode
 * 
 * deal.II includes
 * 
 * @code
 * #include <deal.II/lac/vector.h>
 * 
 * 
 * namespace step42 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassGPUInfo"></a> 
 * <h3>Class: GPUInfo</h3>
 * 
 * 
 * This an auxiliary structure which collects data related to the GPUs
 * available in the computer on which this program gets executed.
 * For details about the cuda<something> functions have a look at
 * the CUDA reference manual.
 * 
 * @code
 * struct GPUInfo {
 * 
 *     int n_CUDA_devices;
 * 
 *     int current_device_id;
 * 
 *     cudaDeviceProp prop;
 * 
 * 
 *     GPUInfo(int DevNo)
 *         :
 *           current_device_id(DevNo)
 *     {
 * 
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Functionget"></a> 
 * <h4>Function: get</h4>
 * 
 * 
 * this function must be called to retrieve the
 * GPU-related information.
 * 
 * @code
 *     void get()
 *     {
 *         static const int KB = 1024;
 *         static const int MB = KB*KB;
 * 
 * @endcode
 * 
 * Retrieve information about the currently selected GPU.
 * 
 * @code
 *         std::cout << "current device ID : " << this->current_device_id << std::endl;
 * 
 *         cudaGetDeviceProperties(&prop, this->current_device_id);
 * 
 *         printf("Currently used GPU: %s \n",prop.name);
 *         printf("Compute Capability: %d.%d \n",prop.major,prop.minor);
 *         printf("ClockRate: %uMHz \n",prop.clockRate/1000);
 *         printf("Warpsize: %d \n",prop.warpSize);
 *         printf("Number of Multiprocessors: %d \n",prop.multiProcessorCount);
 * 
 *         printf("Shared Memory: %luKB\n",prop.sharedMemPerBlock/KB);
 *         printf("Constant Memory: %luKB \n",prop.totalConstMem/KB);
 *         printf("Global Memory: %luMB \n",prop.totalGlobalMem/MB);
 *         printf("the device %s concurrently copy memory between host and device while executing a kernel\n",
 *                (prop.deviceOverlap? "can": "cannot"));
 *     }
 * 
 * @endcode
 * 
 * To keep the compiler from automatically generating
 * a copy constructor and an assignment operator we provide
 * dummy implementations and declare them as private.
 * In case one of them is needed the compiler will complain at compile-time
 * and one can think about whether they are really needed, i.e. one has to
 * review one's software design.
 * 
 * @code
 * private:
 *     GPUInfo (const GPUInfo & / *other* /) {}
 * 
 *     GPUInfo & operator = (const GPUInfo & / *other* /) { return *this; }
 * 
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ClassMyFancySimulation"></a> 
 * <h3>Class: MyFancySimulation</h3>
 * 
 * 
 * To make this test facility extendible, we implement
 * a class for a simple user interface. Its primary tasks are
 * - management of run-time parameters by a simple text-based parameter file
 * - setting device parameters according to the user's parameters
 * - preprocessing and output of results
 * 
 * @code
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     MyFancySimulation(int argc, char *argv[], GPUInfo &g);
 * 
 *     void run();
 * 
 * private:
 *     GPUInfo & gpuinfo;
 * 
 * protected:
 *     SimParams params;
 * };
 * 
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorMyFancySimulation"></a> 
 * <h4>Constructor: MyFancySimulation</h4>
 * 
 * 
 * The constructor is responsible for reading parameters
 * and initializing the device, i.e. the selected graphics card.
 * @param argc : The number of command line arguments. This is always $\ge 1$, as by default the zeroth argument is the name of program itself.
 * @param argv : Pointer to the array of command line arguments.
 * @param g : Reference to the object containing the GPU info from the system.
 * 
 * @code
 * step42::MyFancySimulation::MyFancySimulation(int argc,
 *                                                    char *argv[],
 *                                                    step42::GPUInfo &g)
 *     : gpuinfo(g)
 * {
 * @endcode
 * 
 * Before setting up the simulation we
 * figure out how many GPUs are available
 * 
 * @code
 *     cudaGetDeviceCount(&gpuinfo.n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices : "
 *             << gpuinfo.n_CUDA_devices << std::endl;
 * 
 * @endcode
 * 
 * Declare and read parameters from a file. Basically, the parameter
 * file must have the same name as the binary. The extension has
 * to be ".prm". What has been read will be dumped into a log file.
 * 
 * @code
 *     dealii::ParameterHandler prm_handler;
 * 
 *    SimParams::declare(prm_handler);
 * 
 * @endcode
 * 
 * Get the current working directory ...
 * 
 * @code
 *     QDir cwd = QDir::current();
 * 
 * @endcode
 * 
 * i.e. where the program has been started.
 * 
 * @code
 *     QDir launch_dir = cwd;
 * 
 * @endcode
 * 
 * By default, the parameter file has the same name as the binary
 * and is supposed to be in a subdirectory prm of that directory,
 * where the program has been started.
 * 
 * @code
 *     std::string prm_filename;
 *     if (argc == 1)
 *     {
 *         prm_filename  = argv[0];
 *         prm_filename += ".prm";
 * 
 *         cwd.setPath("./prm");
 *     }
 *     else
 *     {
 * @endcode
 * 
 * Whatever gets passed as first command line argument is considered as path
 * to a parameter file.
 * 
 * @code
 *         std::cout << "Given parameter file : " << argv[1] << std::endl;
 * 
 * @endcode
 * 
 * We convert the sequence of characters into something more meaningful.
 * 
 * @code
 *         QFileInfo tmp(argv[1]);
 * 
 * @endcode
 * 
 * Before we proceed, let us figure out whether the given parameter file exists.
 * Note: If the file is a symlink that points to a non existing file,
 * false is returned as well.
 * 
 * @code
 *         if(!tmp.exists())
 *         {
 *             std::cerr << "The following parameter file does not exist:\n"
 *                       << argv[1] << std::endl;
 * 
 *             qFatal("Cannot proceed without proper path to paramter file");
 *         }
 * 
 * @endcode
 * 
 * Next, we subdivide the given filename into its path and filename
 * so that the corresponding subdirectories can be created.
 * 
 * @code
 *         QString prm_path = tmp.absolutePath();
 *         cwd.setPath(prm_path);
 *         cwd.makeAbsolute();
 *         prm_filename = tmp.fileName().toStdString();
 * 
 *         std::cout << "Parameter file path : "
 *                   << tmp.absolutePath().toStdString().c_str()
 *                   << std::endl;
 *     }
 * 
 *     std::cout << "Parameter file : " << prm_filename  << std::endl;
 * 
 * @endcode
 * 
 * Before the parameter file can be read, we have to make sure that
 * its directory exists. In case of the default parameter file
 * the directory will be created.
 * 
 * @code
 *     if (!cwd.exists() )
 *         launch_dir.mkpath( cwd.absolutePath() );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_handler.read_input (prm_filename);
 * 
 *     QDir::setCurrent(launch_dir.absolutePath());
 * 
 *     this->params.get(prm_handler);
 * 
 * @endcode
 * 
 * Create toplevel run directory
 * 
 * @code
 *     cwd.setPath(this->params.run_dir.absolutePath());
 * 
 * @endcode
 * 
 * The following lets a directory make its own path.
 * 
 * @code
 *     if (!cwd.exists())
 *         cwd.mkpath( "." );
 * 
 * @endcode
 * 
 * Now, change to the run dir
 * 
 * @code
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     cwd.setPath("./log");
 *     cwd.makeAbsolute();
 *     if (!cwd.exists())
 *         cwd.mkpath(".");
 * 
 * @endcode
 * 
 * Create the log directory and write what has been actually read
 * into log file. Basically, this is just another parameter file
 * and can thus be used again as input to another run after stripping the .log suffix.
 * 
 * @code
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(("./" + QString(prm_filename.c_str()).split("/").last()).toStdString().c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 * @endcode
 * 
 * At this point the toplevel run dir must exist.
 * Thus, we can change to it without any further sanity test.
 * 
 * @code
 *     QDir::setCurrent(this->params.run_dir.absolutePath());
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionrun"></a> 
 * <h4>Function: run</h4>
 * 
 * 
 * Actual call to run function.
 * 
 * @code
 * void step42::MyFancySimulation::run()
 * {
 * 
 * @endcode
 * 
 * instantiate an object of the driver class
 * 
 * @code
 *     step42::CUDADriver testcase;
 * 
 * @endcode
 * 
 * ... and run the computation on the GPU (or other dedicated parallel hardware).
 * 
 * @code
 *     testcase.gemm_tests();
 * 
 *     testcase.gemv_tests();
 * 
 *     testcase.complex_tests();
 * 
 *     testcase.feature_demonstration();
 * 
 *     std::cout << "Done." << std::endl;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functionmain"></a> 
 * <h3>Function: main</h3>
 * 
 * 
 * As usual, the main function is pretty boring.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *     using namespace step42;
 * 
 * @endcode
 * 
 * At the beginning we figure out
 * how many CUDA devices are available.
 * 
 * @code
 *     int n_CUDA_devices;
 * 
 *     cudaGetDeviceCount(&n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices : "
 *             <<  n_CUDA_devices
 *             << std::endl;
 * 
 * @endcode
 * 
 * This command is used to set the GPU on which
 * we want to run our computations.
 * For a list of GPUs execute nvidia-smi.
 * 
 * @code
 *     int DevNo = 0;
 *     cudaSetDevice(DevNo);
 *     GPUInfo gpu_info(DevNo);
 * 
 * @endcode
 * 
 * Before we can instantiate the simulation we have to get the GPU info.
 * 
 * @code
 *     gpu_info.get();
 * 
 *     MyFancySimulation machma(argc, argv, gpu_info);
 * 
 *     machma.run();
 * 
 * }
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_42_H
 * #define CUDADriver_STEP_42_H
 * 
 * #include <lac/release/blas_wrapper.hh>
 * #include <lac/release/cublas_wrapper.hh>
 * #include <base/CudaComplex.h>
 * 
 * @endcode
 * 
 * We encapsulate each project
 * into a dedicated namespace
 * in order to be able to re-use
 * parts of a test program in others.
 * 
 * @code
 * namespace step42 {
 * 
 * @endcode
 * 
 * 
 * <a name="ClassCUDADriver"></a> 
 * <h3>Class: CUDADriver</h3>
 * 
 * 
 * This class manages the communication between host and device.
 * In particular the issue of memory transfers from and to the device.
 * The dummy implementation given here is supposed to give an
 * impression how this management could be done.
 * For worked out examples have a look at the other steps from
 * previous lab courses.
 * The documentation of the member functions is kept together
 * with their definitions.
 * 
 * @code
 * class CUDADriver {
 * 
 *     typedef double Number;
 *     typedef SciPAL::CudaComplex<Number> cplxNumber;
 * 
 *     typedef cublas BW;
 * 
 * public:
 * 
 *     CUDADriver();
 * 
 *     ~CUDADriver() { BW::Shutdown(); }
 * 
 *     void gemm_tests();
 * 
 *     void gemv_tests();
 * 
 *     void complex_tests();
 * 
 *     void feature_demonstration();
 * 
 * private:
 * 
 * 
 * };
 * 
 * } // namespace step42 END
 * 
 * #endif // CUDADriver_STEP_42_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_42_HH
 * #define CUDA_DRIVER_STEP_42_HH
 * 
 * 
 * @endcode
 * 
 * The declaration of the interface to the CUDA-backend
 * is contained in the following header.
 * 
 * @code
 * #include <cuda_driver_step-42.h>
 * #include <cuda_kernel_wrapper_step-42.cu.h>
 * 
 * 
 * 
 * #include <lac/development/cublas_Matrix.h>
 * 
 * 
 * #include <lac/development/cublas_Vector.h>
 * #include <lac/blas++.h>
 * #include <base/CudaComplex.h>
 * 
 * @endcode
 * 
 * We have to include
 * 
 * @code
 * #include <cuda_runtime_api.h>
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstructorCUDADriver"></a> 
 * <h4>Constructor: CUDADriver</h4>
 * 
 * 
 * The constructor of the driver class allocates memory
 * on the GPU and copies data from host to device.
 * Furthermore, it keeps a pointer to the original host data for the
 * case that it has to be modified by the GPU results.
 * @param v_h : Pointer to a linear array in host-side memory that is to be copied to the GPU.
 * @param n : Number of entries of @p v_h.
 * 
 * @code
 * step42::CUDADriver::CUDADriver() {
 * 
 *     BW::Init();
 * }
 * 
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functiongemm_tests"></a> 
 * <h4>Function: gemm_tests</h4>
 * 
 * 
 * This function tests the various special cases contained
 * in BLAS' gemm function.
 * 
 * @code
 * void step42::CUDADriver::gemm_tests()
 * {
 * 
 * @endcode
 * 
 * some dummy vectors
 * 
 * @code
 *     const unsigned int n_rows = 4;
 *     const unsigned int n_cols = 3;
 *     const unsigned int n_elements = n_rows * n_cols;
 * 
 *     std::vector<Number>
 *             a(n_elements, 1.),
 *             b(n_elements, 2.),
 *             c(n_rows * n_rows, 1.23);
 * 
 *     for (unsigned int i = 0; i < b.size(); i++ )
 *         b[i] = i+1;
 * 
 * 
 *     SciPAL::Matrix<Number, BW>
 *             A(n_rows, n_cols, a),
 *             B(n_cols, n_rows, b),
 *             C(n_rows, n_rows, c);
 * 
 *      Number alpha = 1.1;
 *      Number beta = 2.;
 * 
 *      std::cout << "A : " << std::endl;
 *      A.print();
 * 
 *      std::cout << "B : " << std::endl;
 *      B.print();
 * 
 *      std::cout << "C : " << std::endl;
 *      C.print();
 * 
 * 
 * 
 *      std::cout << " ============ C = " << alpha << " * A ======" << std::endl;
 *      C = alpha * A; // * B + beta * C;
 *        std::cout << "C : " << std::endl;
 *      C.print();
 * 
 *      std::cout << " ============ C = A * B ======" << std::endl;
 *      C = A * B;   std::cout << "C : " << std::endl; C.print();
 * 
 *      std::cout << " ============ C = B * A ======" << std::endl;
 *      C = B * A;   std::cout << "C : " << std::endl; C.print();
 * 
 *      std::cout << " ============ C = alpha * A * B ======" << std::endl;
 *      C = alpha * A * B; std::cout << "C : " << std::endl; C.print();
 * 
 * 
 * @endcode
 * 
 * sMMaM test
 * 
 * @code
 *      std::cout << " ============ C = " << alpha << " * A * B + "<<"C======" << std::endl;
 *      c.clear();
 *      c.resize(n_rows * n_rows, 1.);
 *      SciPAL::Matrix<Number, BW>
 *              D(n_rows, n_rows, c);
 * 
 * 
 * 
 *      C = D;
 *      std::cout << "C : " << std::endl; C.print();
 *      C = alpha * A * B + // beta *
 *              C; std::cout << "C : " << std::endl; C.print();
 * 
 * 
 * @endcode
 * 
 * gemm test
 * 
 * @code
 *      std::cout << " ============ C = " << alpha << " * A * B + "  << beta  << " * C======" << std::endl;
 *      C = D;
 *      std::cout << "C : " << std::endl; C.print();
 * 
 *      C = alpha * A * B + beta * C;
 *      std::cout << "C : " << std::endl; C.print();
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Functiongemv_tests"></a> 
 * <h4>Function: gemv_tests</h4>
 * 
 * 
 * This function tests the various special cases contained
 * in BLAS' gemv function.
 * 
 * 
 * Currently, it is rather a test for the vector arithmetic.
 * test vector expressions
 * 
 * @code
 * void step42::CUDADriver::gemv_tests()
 * {
 * #ifndef nUSE_ARRAY_EXPRESSIONS
 *      const unsigned int n_rows = 4;
 *      const unsigned int n_cols = 4;
 *      const unsigned int n_elements = n_rows * n_cols;
 * 
 *      Number alpha = 1.1;
 *      Number beta = 2.;
 * 
 *      std::vector<Number>
 *                 a(n_elements, 1.),
 *                 b(n_elements, 2.);
 * 
 * 
 *      for (unsigned int i = 0; i < a.size(); i++ )
 *          a[i] = i+1;
 * 
 *      SciPAL::Vector<Number, BW> vA, vB(n_elements), vC;
 *        vA = a;
 * @endcode
 * 
 * This sets all elements of vB to 2.3, note: vector needs to be initialized.
 * 
 * @code
 *        vB = SciPAL::Literal<Number>(2.3);
 *        vC = a;
 * 
 * 
 *        std::cout << "vA : " << std::endl;
 *        vA.print();
 * 
 *        std::cout << "vB : " << std::endl;
 *        vB.print();
 * 
 *        std::cout << " ============ vC = " << alpha << " * vA ======" << std::endl;
 *        vC = alpha * vA;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vC = " << alpha << " * vA + vB ======" << std::endl;
 *        vC = alpha * vA + vB;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vA = sin(vC) ======" << std::endl;
 *        const unsigned int n_sin_elements = n_elements;
 *        std::vector<Number> d(n_sin_elements);
 *        for(uint i = 0; i < d.size(); i++)
 *            d[i] = i* 2.* M_PI / d.size();
 * 
 *        SciPAL::Vector<Number, BW> vD; vD = d; //(n_sin_elements, 1, d);
 *       vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
 *          std::cout << "sin(vD) : " << std::endl;
 *        vD.print();
 * @endcode
 * 
 * vC = alpha * sin(vA) + vB;
 * 

 * 
 * 

 * 
 * After the element-wise sine of a vector we do the same for a matrix.
 * 
 * @code
 *        SciPAL::Matrix<Number, BW>
 *                A(n_rows, n_cols, d);
 *        A = sin(A);
 *        A.print();
 * 
 * 
 *        std::cout << " ============ linear combination test ======" << std::endl;
 *        vC = 2.0 * vA;
 *        std::cout << "vC = 2.0 * vA" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + vB;
 *        std::cout << "vC = vA + vB" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + 2.0 * vB;
 *        std::cout << "vC = vA + 2.0 * vB" << std::endl;
 *        vC.print();
 * 
 *        vC = 2.0 * vA + 3.0 * vB;
 *        std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * vC = 2.0 * vA + 3.0 * vB + 4.0 * vD;
 * std::cout << "vC = 2.0 * vA + 3.0 * vB + 4.0 * vD" << std::endl;
 * vC.print();
 * 

 * 
 * combined expr test
 * 
 * @code
 *        vC =sin(2.0 * vA + 3.0 * vB);
 *        std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise sqrt
 * 
 * @code
 *        vC = sqrt(vC);
 *        std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise *
 * 
 * @code
 *        vC = vA && vB;
 *        std::cout << "vC = vA .* vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise /
 * 
 * @code
 *        vC = vA || vB;
 *        std::cout << "vC = vA ./ vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * combined expr test
 * 
 * @code
 *        vC = (vA + vB) || (vA - vB);
 *        std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * combined expr test
 * vC =abs(cos(2.0 * vA + 3.0 * vB));
 * std::cout << "vC = abs(cos(2.0 * vA + 3.0 * vB)" << std::endl;
 * vC.print();
 * 

 * 
 * 
 * @code
 * #endif
 * }
 * 
 * void step42::CUDADriver::complex_tests()
 * {
 *     std::cout<<"Entering tests for complex number array exptessions."<<std::endl;
 * #ifndef nUSE_ARRAY_EXPRESSIONS
 *      const unsigned int n_rows = 4;
 *      const unsigned int n_cols = 4;
 *      const unsigned int n_elements = n_rows * n_cols;
 * 
 *      SciPAL::CudaComplex<Number> alpha(1.1);
 *      SciPAL::CudaComplex<Number> beta (2.);
 * 
 *      std::vector<std::complex<Number> >
 *                 a(n_elements, 1.),
 *                 b(n_elements, std::complex<Number>(2., 2.0));
 * 
 * 
 *      for (unsigned int i = 0; i < a.size(); i++ )
 *          a[i] = std::complex<Number>(i+1, (i+1)/2.); //generate some inputs
 * 
 *        SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vA, vB, vC;
 *        vA = a;
 *        vB = b;
 *        vC = a;
 * @endcode
 * 
 * vA(n_elements, 1, a),
 * vB(n_elements, 1, b),
 * vC(n_elements, 1, a);
 * 

 * 
 * 
 * @code
 *        std::cout << "vA : " << std::endl;
 *        vA.print();
 * 
 *        std::cout << "vB : " << std::endl;
 *        vB.print();
 * 
 *        std::cout << " ============ vC = " << alpha.real() << " * vA ======" << std::endl;
 *        vC = alpha * vA;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vC = " << alpha.real() << " * vA + vB ======" << std::endl;
 *        vC = alpha * vA + vB;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vA = sin(vC) ======" << std::endl;
 *        const unsigned int n_sin_elements = n_elements;
 *        std::vector<std::complex<Number> > d(n_sin_elements);
 *        for(uint i = 0; i < d.size(); i++)
 *            d[i] = std::complex<Number>(i* 2.* M_PI / d.size(), i* 4.* M_PI / d.size()) ;
 * 
 *        SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vD; vD = d; //(n_sin_elements, 1, d);
 *       vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
 *          std::cout << "sin(vD) : " << std::endl;
 *        vD.print();
 * @endcode
 * 
 * vC = alpha * sin(vA) + vB;
 * 

 * 
 * 
 * @code
 *        std::cout << " ============ Matrix A = sin(A) ======" << std::endl;
 * @endcode
 * 
 * After the element-wise sine of a vector we do the same for a matrix.
 * 
 * @code
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                A(n_rows, n_cols, d);
 *        A = sin(A);
 *        A.print();
 * 
 *        std::cout << " ============ Matrix B = sqrt(A) ======" << std::endl;
 * @endcode
 * 
 * After the element-wise sine of a vector we do the same for a matrix.
 * 
 * @code
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                B(A);
 *        B = sqrt(A);
 *        B.print();
 * 
 *        std::cout << " ============ Matrix A = .exp(B) ======" << std::endl;
 * @endcode
 * 
 * After the element-wise sine of a vector we do the same for a matrix.
 * SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 * A(B);
 * 
 * @code
 *        A = exp(B);
 *        A.print();
 * 
 *        std::cout << " ============ Matrix A = B*C ======" << std::endl;
 * @endcode
 * 
 * After the element-wise sine of a vector we do the same for a matrix.
 * 
 * @code
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                C(A);
 *        A = B*C;
 *        A.print();
 * 
 * 
 *        std::cout << " ============ linear combination test ======" << std::endl;
 * 
 * @endcode
 * 
 * This does not work : vC = 2.0 * vA; because of mismatching type.
 * Thus we hav to use numbers wrapped in CudaComplexes
 * 
 * @code
 *        vC = beta * vA; ;
 *        std::cout << "vC = 2.0 * vA" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + vB;
 *        std::cout << "vC = vA + vB" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + beta * vB;
 *        std::cout << "vC = vA + 2.0 * vB" << std::endl;
 *        vC.print();
 * 
 *        vC = cplxNumber(2.0) * vA + cplxNumber(3.0) * vB;
 *        std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * combined expr test
 * 
 * @code
 *        vC =sin(cplxNumber(2.0) * vA + cplxNumber(3.0) * vB);
 *        std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise sqrt
 * 
 * @code
 *        vC = sqrt(vC);
 *        std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise *
 * 
 * @code
 *        vC = vA && vB;
 *        std::cout << "vC = vA .* vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * test pointwise /
 * 
 * @code
 *        vC = vA || vB;
 *        std::cout << "vC = vA ./ vB" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * combined expr test
 * 
 * @code
 *        vC = (vA + vB) || (vA - vB);
 *        std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
 *        vC.print();
 * 
 * @endcode
 * 
 * combined expr test
 * vC =abs(cos(2.0 * vA + 3.0 * vB));
 * std::cout << "vC = abs(cos(2.0 * vA + 3.0 * vB)" << std::endl;
 * vC.print();
 * 

 * 
 * 
 * @code
 * #endif
 * }
 * 
 * 
 * void step42::CUDADriver::feature_demonstration()
 * {
 *     Number * bla = new Number[3];
 *     Number * bla2 = new Number[3];
 *     Number * bla3 = new Number[3];
 *     std::vector<Number*> h_testV(5);
 *     h_testV[0] = bla; std::cout<<"ptr1 " << bla << std::endl;
 *     h_testV[1] = bla2;std::cout<<"ptr2 " << bla2 << std::endl;
 *     h_testV[2] = bla3;std::cout<<"ptr3 " << bla3 << std::endl;
 * 
 *     SciPAL::Vector<Number*, cublas> d_testV(5);
 * 
 *  d_testV = h_testV;
 *  d_testV.print();
 * 
 * 
 * 
 * 
 * }
 * 
 * #endif // CUDA_DRIVER_STEP_42_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * 
 * 
 * 
 * @endcode
<a name="Results"></a><h1>Results</h1><br>

The programm generates the following output:<br>
<br>
N available CUDA devices : 4<br>
current device ID : 0<br>
Currently used GPU: GeForce GTX TITAN <br>
Compute Capability: 3.5 <br>
ClockRate: 875MHz <br>
Warpsize: 32 <br>
Number of Multiprocessors: 14 <br>
Shared Memory: 48KB<br>
Constant Memory: 64KB <br>
Global Memory: 6143MB <br>
the device can concurrently copy memory between host and device while executing a kernel<br>
N available CUDA devices : 4<br>
Parameter file : /home/jhagemann/Documents/CUDA/Praktikum_2014/testsite/build-step-42-Desktop/step-42.prm<br>
cublas init succeeded<br>
A : <br>
Matrix dims : 4 3<br>
         1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000 ;<br>
B : <br>
Matrix dims : 3 4<br>
         1.0000          4.0000          7.0000         10.0000 ;<br>
         2.0000          5.0000          8.0000         11.0000 ;<br>
         3.0000          6.0000          9.0000         12.0000 ;<br>
C : <br>
Matrix dims : 4 4<br>
         1.2300          1.2300          1.2300          1.2300 ;<br>
         1.2300          1.2300          1.2300          1.2300 ;<br>
         1.2300          1.2300          1.2300          1.2300 ;<br>
         1.2300          1.2300          1.2300          1.2300 ;<br>
 ============ C = 1.1000 * A ======<br>
C : <br>
Matrix dims : 4 3<br>
         1.1000          1.1000          1.1000 ;<br>
         1.1000          1.1000          1.1000 ;<br>
         1.1000          1.1000          1.1000 ;<br>
         1.1000          1.1000          1.1000 ;<br>
 ============ C = A * B ======<br>
C : <br>
Matrix dims : 4 4<br>
         6.0000         15.0000         24.0000         33.0000 ;<br>
         6.0000         15.0000         24.0000         33.0000 ;<br>
         6.0000         15.0000         24.0000         33.0000 ;<br>
         6.0000         15.0000         24.0000         33.0000 ;<br>
 ============ C = B * A ======<br>
C : <br>
Matrix dims : 3 3<br>
        22.0000         22.0000         22.0000 ;<br>
        26.0000         26.0000         26.0000 ;<br>
        30.0000         30.0000         30.0000 ;<br>
 ============ C = alpha * A * B ======<br>
C : <br>
Matrix dims : 4 4<br>
         6.6000         16.5000         26.4000         36.3000 ;<br>
         6.6000         16.5000         26.4000         36.3000 ;<br>
         6.6000         16.5000         26.4000         36.3000 ;<br>
         6.6000         16.5000         26.4000         36.3000 ;<br>
 ============ C = 1.1000 * A * B + C======<br>
C : <br>
Matrix dims : 4 4<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
C : <br>
Matrix dims : 4 4<br>
         7.6000         17.5000         27.4000         37.3000 ;<br>
         7.6000         17.5000         27.4000         37.3000 ;<br>
         7.6000         17.5000         27.4000         37.3000 ;<br>
         7.6000         17.5000         27.4000         37.3000 ;<br>
 ============ C = 1.1000 * A * B + 2.0000 * C======<br>
C : <br>
Matrix dims : 4 4<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
         1.0000          1.0000          1.0000          1.0000 ;<br>
C : <br>
Matrix dims : 4 4<br>
         8.6000         18.5000         28.4000         38.3000 ;<br>
         8.6000         18.5000         28.4000         38.3000 ;<br>
         8.6000         18.5000         28.4000         38.3000 ;<br>
         8.6000         18.5000         28.4000         38.3000 ;<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
vA : <br>
1.0000<br>
2.0000<br>
3.0000<br>
4.0000<br>
5.0000<br>
6.0000<br>
7.0000<br>
8.0000<br>
9.0000<br>
10.0000<br>
11.0000<br>
12.0000<br>
13.0000<br>
14.0000<br>
15.0000<br>
16.0000<br>
vB : <br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
2.3000<br>
 ============ vC = 1.1000 * vA ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC : <br>
1.1000<br>
2.2000<br>
3.3000<br>
4.4000<br>
5.5000<br>
6.6000<br>
7.7000<br>
8.8000<br>
9.9000<br>
11.0000<br>
12.1000<br>
13.2000<br>
14.3000<br>
15.4000<br>
16.5000<br>
17.6000<br>
 ============ vC = 1.1000 * vA + vB ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC : <br>
3.4000<br>
4.5000<br>
5.6000<br>
6.7000<br>
7.8000<br>
8.9000<br>
10.0000<br>
11.1000<br>
12.2000<br>
13.3000<br>
14.4000<br>
15.5000<br>
16.6000<br>
17.7000<br>
18.8000<br>
19.9000<br>
 ============ vA = sin(vC) ======<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
sin(vD) : <br>
0.0000<br>
0.3827<br>
0.7071<br>
0.9239<br>
1.0000<br>
0.9239<br>
0.7071<br>
0.3827<br>
0.0000<br>
-0.3827<br>
-0.7071<br>
-0.9239<br>
-1.0000<br>
-0.9239<br>
-0.7071<br>
-0.3827<br>
Matrix dims : 4 4<br>
         0.0000          1.0000          0.0000         -1.0000 ;<br>
         0.3827          0.9239         -0.3827         -0.9239 ;<br>
         0.7071          0.7071         -0.7071         -0.7071 ;<br>
         0.9239          0.3827         -0.9239         -0.3827 ;<br>
 ============ linear combination test ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = 2.0 * vA<br>
2.0000<br>
4.0000<br>
6.0000<br>
8.0000<br>
10.0000<br>
12.0000<br>
14.0000<br>
16.0000<br>
18.0000<br>
20.0000<br>
22.0000<br>
24.0000<br>
26.0000<br>
28.0000<br>
30.0000<br>
32.0000<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA + vB<br>
3.3000<br>
4.3000<br>
5.3000<br>
6.3000<br>
7.3000<br>
8.3000<br>
9.3000<br>
10.3000<br>
11.3000<br>
12.3000<br>
13.3000<br>
14.3000<br>
15.3000<br>
16.3000<br>
17.3000<br>
18.3000<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA + 2.0 * vB<br>
5.6000<br>
6.6000<br>
7.6000<br>
8.6000<br>
9.6000<br>
10.6000<br>
11.6000<br>
12.6000<br>
13.6000<br>
14.6000<br>
15.6000<br>
16.6000<br>
17.6000<br>
18.6000<br>
19.6000<br>
20.6000<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = 2.0 * vA + 3.0 * vB<br>
8.9000<br>
10.9000<br>
12.9000<br>
14.9000<br>
16.9000<br>
18.9000<br>
20.9000<br>
22.9000<br>
24.9000<br>
26.9000<br>
28.9000<br>
30.9000<br>
32.9000<br>
34.9000<br>
36.9000<br>
38.9000<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
vC = sin(2.0 * vA + 3.0 * vB)<br>
0.5010<br>
-0.9954<br>
0.3275<br>
0.7229<br>
-0.9291<br>
0.0504<br>
0.8872<br>
-0.7888<br>
-0.2306<br>
0.9808<br>
-0.5856<br>
-0.4933<br>
0.9962<br>
-0.3358<br>
-0.7167<br>
0.9324<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
sqrt(sin(2.0 * vA + 3.0 * vB))<br>
0.7078<br>
-nan<br>
0.5723<br>
0.8502<br>
-nan<br>
0.2245<br>
0.9419<br>
-nan<br>
-nan<br>
0.9903<br>
-nan<br>
-nan<br>
0.9981<br>
-nan<br>
-nan<br>
0.9656<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA .* vB<br>
2.3000<br>
4.6000<br>
6.9000<br>
9.2000<br>
11.5000<br>
13.8000<br>
16.1000<br>
18.4000<br>
20.7000<br>
23.0000<br>
25.3000<br>
27.6000<br>
29.9000<br>
32.2000<br>
34.5000<br>
36.8000<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA ./ vB<br>
0.4348<br>
0.8696<br>
1.3043<br>
1.7391<br>
2.1739<br>
2.6087<br>
3.0435<br>
3.4783<br>
3.9130<br>
4.3478<br>
4.7826<br>
5.2174<br>
5.6522<br>
6.0870<br>
6.5217<br>
6.9565<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = (vA + vB) || (vA - vB)<br>
-2.5385<br>
-14.3333<br>
7.5714<br>
3.7059<br>
2.7037<br>
2.2432<br>
1.9787<br>
1.8070<br>
1.6866<br>
1.5974<br>
1.5287<br>
1.4742<br>
1.4299<br>
1.3932<br>
1.3622<br>
1.3358<br>
Entering tests for complex number array exptessions.<br>
vA : <br>
(1.0000,0.5000)<br>
(2.0000,1.0000)<br>
(3.0000,1.5000)<br>
(4.0000,2.0000)<br>
(5.0000,2.5000)<br>
(6.0000,3.0000)<br>
(7.0000,3.5000)<br>
(8.0000,4.0000)<br>
(9.0000,4.5000)<br>
(10.0000,5.0000)<br>
(11.0000,5.5000)<br>
(12.0000,6.0000)<br>
(13.0000,6.5000)<br>
(14.0000,7.0000)<br>
(15.0000,7.5000)<br>
(16.0000,8.0000)<br>
vB : <br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
(2.0000,2.0000)<br>
 ============ vC = 1.1000 * vA ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC : <br>
(1.1000,0.5500)<br>
(2.2000,1.1000)<br>
(3.3000,1.6500)<br>
(4.4000,2.2000)<br>
(5.5000,2.7500)<br>
(6.6000,3.3000)<br>
(7.7000,3.8500)<br>
(8.8000,4.4000)<br>
(9.9000,4.9500)<br>
(11.0000,5.5000)<br>
(12.1000,6.0500)<br>
(13.2000,6.6000)<br>
(14.3000,7.1500)<br>
(15.4000,7.7000)<br>
(16.5000,8.2500)<br>
(17.6000,8.8000)<br>
 ============ vC = 1.1000 * vA + vB ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC : <br>
(3.1000,2.5500)<br>
(4.2000,3.1000)<br>
(5.3000,3.6500)<br>
(6.4000,4.2000)<br>
(7.5000,4.7500)<br>
(8.6000,5.3000)<br>
(9.7000,5.8500)<br>
(10.8000,6.4000)<br>
(11.9000,6.9500)<br>
(13.0000,7.5000)<br>
(14.1000,8.0500)<br>
(15.2000,8.6000)<br>
(16.3000,9.1500)<br>
(17.4000,9.7000)<br>
(18.5000,10.2500)<br>
(19.6000,10.8000)<br>
 ============ vA = sin(vC) ======<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
sin(vD) : <br>
(0.0000,0.0000)<br>
(0.5069,0.8025)<br>
(1.7743,1.6273)<br>
(4.9176,2.0007)<br>
(11.5920,0.0000)<br>
(23.4544,-9.7076)<br>
(39.3600,-39.3536)<br>
(46.7171,-112.7812)<br>
(0.0000,-267.7449)<br>
(-224.7278,-542.5401)<br>
(-910.7432,-910.7430)<br>
(-2609.8788,-1081.0471)<br>
(-6195.8239,-0.0000)<br>
(-12554.7625,5200.3529)<br>
(-21075.2262,21075.2262)<br>
(-25016.1799,60394.4008)<br>
 ============ Matrix A = sin(A) ======<br>
Matrix dims : 4 4<br>
              (0.0000,0.0000)               (11.5920,0.0000)               (0.0000,-267.7449)               (-6195.8239,-0.0000) ;<br>
              (0.5069,0.8025)               (23.4544,-9.7076)               (-224.7278,-542.5401)               (-12554.7625,5200.3529) ;<br>
              (1.7743,1.6273)               (39.3600,-39.3536)               (-910.7432,-910.7430)               (-21075.2262,21075.2262) ;<br>
              (4.9176,2.0007)               (46.7171,-112.7812)               (-2609.8788,-1081.0471)               (-25016.1799,60394.4008) ;<br>
 ============ Matrix B = sqrt(A) ======<br>
Matrix dims : 4 4<br>
              (0.0000,0.0000)               (3.4047,0.0000)               (11.5703,-11.5703)               (0.0000,-78.7136) ;<br>
              (0.8533,0.4703)               (4.9416,-0.9822)               (13.4632,-20.1491)               (22.7422,114.3327) ;<br>
              (1.4460,0.5627)               (6.8927,-2.8547)               (13.7339,-33.1567)               (66.0668,159.4994) ;<br>
              (2.2613,0.4424)               (9.1867,-6.1383)               (10.3690,-52.1286)               (142.0462,212.5872) ;<br>
 ============ Matrix A = .exp(B) ======<br>
Matrix dims : 4 4<br>
              (1.0000,0.0000)               (30.1051,0.0000)               (57574.4027,88891.4742)               (-0.9849,0.1729) ;<br>
              (2.0925,1.0636)               (77.7179,-116.4362)               (188398.4991,-677321.8222)               (2478670308.2967,7110591049.5679) ;<br>
              (3.5914,2.2651)               (-944.8110,-278.7214)               (-155874.1540,-908387.6504)               (-36969559720391065091893624832.0000,32546478504973833507591684096.0000) ;<br>
              (8.6715,4.1076)               (9664.0784,1410.2392)               (-9181.9538,-30505.7201)               (24739089814809277710757969316927550983550924828268148469792768.0000,-42255371958319387588871286501581214085436953737457640304279552.0000) ;<br>
 ============ Matrix A = B*C ======<br>
Matrix dims : 4 4<br>
              (398.2076,-694.2871)               (97112.8830,-753383.5362)               (-14073630.3226,-10290157.0855)               (-3326070599327342616883763880003935756409361661731012954684391424.0000,-1947301738778217435486587252179419213821722464553705172119846912.0000) ;<br>
              (-166.1929,1046.6511)               (40505.0867,1151639.2786)               (-16849711.7775,-14261783.7370)               (5393792506046452154689315654067589162203322579174980530995200000.0000,1867507568047586354750250687532817671749905479758043067776499712.0000) ;<br>
              (61.0726,1568.4189)               (391571.8626,1661075.7553)               (-28602606.7067,-15832878.4997)               (8374137670507090950721132455600367553100777241987657454835466240.0000,1154191925654070651593339440679955036526394559791946308182343680.0000) ;<br>
              (541.8601,2260.5506)               (1048688.0803,2299606.4475)               (-46124382.0113,-14731102.9816)               (12497043688746812991450891416460595282122078443892479714889039872.0000,-743002158757707251081603709153251540147342286698506061975912448.0000) ;<br>
 ============ linear combination test ======<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = 2.0 * vA<br>
(2.0000,1.0000)<br>
(4.0000,2.0000)<br>
(6.0000,3.0000)<br>
(8.0000,4.0000)<br>
(10.0000,5.0000)<br>
(12.0000,6.0000)<br>
(14.0000,7.0000)<br>
(16.0000,8.0000)<br>
(18.0000,9.0000)<br>
(20.0000,10.0000)<br>
(22.0000,11.0000)<br>
(24.0000,12.0000)<br>
(26.0000,13.0000)<br>
(28.0000,14.0000)<br>
(30.0000,15.0000)<br>
(32.0000,16.0000)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA + vB<br>
(3.0000,2.5000)<br>
(4.0000,3.0000)<br>
(5.0000,3.5000)<br>
(6.0000,4.0000)<br>
(7.0000,4.5000)<br>
(8.0000,5.0000)<br>
(9.0000,5.5000)<br>
(10.0000,6.0000)<br>
(11.0000,6.5000)<br>
(12.0000,7.0000)<br>
(13.0000,7.5000)<br>
(14.0000,8.0000)<br>
(15.0000,8.5000)<br>
(16.0000,9.0000)<br>
(17.0000,9.5000)<br>
(18.0000,10.0000)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA + 2.0 * vB<br>
(5.0000,4.5000)<br>
(6.0000,5.0000)<br>
(7.0000,5.5000)<br>
(8.0000,6.0000)<br>
(9.0000,6.5000)<br>
(10.0000,7.0000)<br>
(11.0000,7.5000)<br>
(12.0000,8.0000)<br>
(13.0000,8.5000)<br>
(14.0000,9.0000)<br>
(15.0000,9.5000)<br>
(16.0000,10.0000)<br>
(17.0000,10.5000)<br>
(18.0000,11.0000)<br>
(19.0000,11.5000)<br>
(20.0000,12.0000)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = 2.0 * vA + 3.0 * vB<br>
(8.0000,7.0000)<br>
(10.0000,8.0000)<br>
(12.0000,9.0000)<br>
(14.0000,10.0000)<br>
(16.0000,11.0000)<br>
(18.0000,12.0000)<br>
(20.0000,13.0000)<br>
(22.0000,14.0000)<br>
(24.0000,15.0000)<br>
(26.0000,16.0000)<br>
(28.0000,17.0000)<br>
(30.0000,18.0000)<br>
(32.0000,19.0000)<br>
(34.0000,20.0000)<br>
(36.0000,21.0000)<br>
(38.0000,22.0000)<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
vC = sin(2.0 * vA + 3.0 * vB)<br>
(542.4820,-79.7800)<br>
(-810.8521,-1250.6183)<br>
(-2173.9477,3418.9097)<br>
(10909.7895,1505.9188)<br>
(-8618.9820,-28669.5197)<br>
(-61113.3864,53734.8541)<br>
(201949.6025,90270.4846)<br>
(-5322.3112,-601278.5869)<br>
(-1480175.6988,693324.2720)<br>
(3388089.3347,2874298.2980)<br>
(3271858.2586,-11625849.6106)<br>
(-32437062.9723,5064072.7195)<br>
(49209951.4403,74447052.4502)<br>
(128346252.3997,-205848381.5924)<br>
(-653986778.5243,-84380263.6616)<br>
(531227762.5025,1711927887.7734)<br>
line :1002, Vector<T,BW>  operator=<br>
<br>
sqrt(sin(2.0 * vA + 3.0 * vB))<br>
(23.3538,-1.7081)<br>
(18.4340,-33.9215)<br>
(30.6398,55.7920)<br>
(104.6972,7.1918)<br>
(103.2426,-138.8453)<br>
(100.6579,266.9183)<br>
(459.9762,98.1252)<br>
(545.8845,-550.7379)<br>
(277.7886,1247.9352)<br>
(1978.7805,726.2802)<br>
(2770.3190,-2098.2872)<br>
(443.2390,5712.5759)<br>
(8320.1894,4473.8797)<br>
(13618.5324,-7557.6566)<br>
(1646.3731,-25626.1063)<br>
(34085.8049,25112.0356)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA .* vB<br>
(1.0000,3.0000)<br>
(2.0000,6.0000)<br>
(3.0000,9.0000)<br>
(4.0000,12.0000)<br>
(5.0000,15.0000)<br>
(6.0000,18.0000)<br>
(7.0000,21.0000)<br>
(8.0000,24.0000)<br>
(9.0000,27.0000)<br>
(10.0000,30.0000)<br>
(11.0000,33.0000)<br>
(12.0000,36.0000)<br>
(13.0000,39.0000)<br>
(14.0000,42.0000)<br>
(15.0000,45.0000)<br>
(16.0000,48.0000)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = vA ./ vB<br>
(0.3750,-0.1250)<br>
(0.7500,-0.2500)<br>
(1.1250,-0.3750)<br>
(1.5000,-0.5000)<br>
(1.8750,-0.6250)<br>
(2.2500,-0.7500)<br>
(2.6250,-0.8750)<br>
(3.0000,-1.0000)<br>
(3.3750,-1.1250)<br>
(3.7500,-1.2500)<br>
(4.1250,-1.3750)<br>
(4.5000,-1.5000)<br>
(4.8750,-1.6250)<br>
(5.2500,-1.7500)<br>
(5.6250,-1.8750)<br>
(6.0000,-2.0000)<br>
line :1015, Vector<T,BW>  operator=<br>
<br>
vC = (vA + vB) || (vA - vB)<br>
(-2.0769,0.6154)<br>
(-3.0000,4.0000)<br>
(2.6000,4.8000)<br>
(3.0000,2.0000)<br>
(2.5135,1.0811)<br>
(2.1765,0.7059)<br>
(1.9541,0.5138)<br>
(1.8000,0.4000)<br>
(1.6878,0.3258)<br>
(1.6027,0.2740)<br>
(1.5362,0.2359)<br>
(1.4828,0.2069)<br>
(1.4389,0.1841)<br>
(1.4024,0.1657)<br>
(1.3714,0.1506)<br>
(1.3448,0.1379)<br>
ptr1 0x15f3470<br>
ptr2 0x15f3490<br>
ptr3 0x4ce8700<br>
0x15f3470<br>
0x15f3490<br>
0x4ce8700<br>
0<br>
0<br>
Done.<br>
cublas shutdown succeeded<br>
 * <a name="PlainProg"></a>
 * <h1> The plain program</h1>
 * 
 * (If you are looking at a locally installed CUDA HPC Praktikum version, then the
 * program can be found at <i>
 *  .. /.. /testsite / /step-42 /step-cu.cc
 * </i>. Otherwise, this is only
 * the path on some remote server.)
 @code

 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #include <iostream>
 * #include <vector>
 * 
 * #include <cuda_driver_step-42.h>
 * #include <cuda_driver_step-42.hh>
 * 
 * #include <SimParams.h>
 * 
 * 
 * #include <deal.II/lac/vector.h>
 * 
 * 
 * namespace step42 {
 * 
@endcode
 <a name="plain-ClassGPUInfo"></a>
@code
 * struct GPUInfo {
 * 
 *     int n_CUDA_devices;
 * 
 *     int current_device_id;
 * 
 *     cudaDeviceProp prop;
 * 
 * 
 *     GPUInfo(int DevNo)
 *         :
 *           current_device_id(DevNo)
 *     {
 * 
 *     }
 * 
@endcode
 <a name="plain-Functionget"></a>
@code
 *     void get()
 *     {
 *         static const int KB = 1024;
 *         static const int MB = KB*KB;
 * 
 *         std::cout << "current device ID : " << this->current_device_id << std::endl;
 * 
 *         cudaGetDeviceProperties(&prop, this->current_device_id);
 * 
 *         printf("Currently used GPU: %s \n",prop.name);
 *         printf("Compute Capability: %d.%d \n",prop.major,prop.minor);
 *         printf("ClockRate: %uMHz \n",prop.clockRate/1000);
 *         printf("Warpsize: %d \n",prop.warpSize);
 *         printf("Number of Multiprocessors: %d \n",prop.multiProcessorCount);
 * 
 *         printf("Shared Memory: %luKB\n",prop.sharedMemPerBlock/KB);
 *         printf("Constant Memory: %luKB \n",prop.totalConstMem/KB);
 *         printf("Global Memory: %luMB \n",prop.totalGlobalMem/MB);
 *         printf("the device %s concurrently copy memory between host and device while executing a kernel\n",
 *                (prop.deviceOverlap? "can": "cannot"));
 *     }
 * 
 * private:
 *     GPUInfo (const GPUInfo & / *other* /) {}
 * 
 *     GPUInfo & operator = (const GPUInfo & / *other* /) { return *this; }
 * 
 * };
 * 
 * 
@endcode
 <a name="plain-ClassMyFancySimulation"></a>
@code
 * class MyFancySimulation {
 * 
 * public:
 * 
 *     MyFancySimulation(int argc, char *argv[], GPUInfo &g);
 * 
 *     void run();
 * 
 * private:
 *     GPUInfo & gpuinfo;
 * 
 * protected:
 *     SimParams params;
 * };
 * 
 * }
 * 
 * 
@endcode
 <a name="plain-ConstructorMyFancySimulation"></a>
@code
 * step42::MyFancySimulation::MyFancySimulation(int argc,
 *                                                    char *argv[],
 *                                                    step42::GPUInfo &g)
 *     : gpuinfo(g)
 * {
 *     cudaGetDeviceCount(&gpuinfo.n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices : "
 *             << gpuinfo.n_CUDA_devices << std::endl;
 * 
 *     dealii::ParameterHandler prm_handler;
 * 
 *    SimParams::declare(prm_handler);
 * 
 *     QDir cwd = QDir::current();
 * 
 *     QDir launch_dir = cwd;
 * 
 *     std::string prm_filename;
 *     if (argc == 1)
 *     {
 *         prm_filename  = argv[0];
 *         prm_filename += ".prm";
 * 
 *         cwd.setPath("./prm");
 *     }
 *     else
 *     {
 *         std::cout << "Given parameter file : " << argv[1] << std::endl;
 * 
 *         QFileInfo tmp(argv[1]);
 * 
 *         if(!tmp.exists())
 *         {
 *             std::cerr << "The following parameter file does not exist:\n"
 *                       << argv[1] << std::endl;
 * 
 *             qFatal("Cannot proceed without proper path to paramter file");
 *         }
 * 
 *         QString prm_path = tmp.absolutePath();
 *         cwd.setPath(prm_path);
 *         cwd.makeAbsolute();
 *         prm_filename = tmp.fileName().toStdString();
 * 
 *         std::cout << "Parameter file path : "
 *                   << tmp.absolutePath().toStdString().c_str()
 *                   << std::endl;
 *     }
 * 
 *     std::cout << "Parameter file : " << prm_filename  << std::endl;
 * 
 *     if (!cwd.exists() )
 *         launch_dir.mkpath( cwd.absolutePath() );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_handler.read_input (prm_filename);
 * 
 *     QDir::setCurrent(launch_dir.absolutePath());
 * 
 *     this->params.get(prm_handler);
 * 
 *     cwd.setPath(this->params.run_dir.absolutePath());
 * 
 *     if (!cwd.exists())
 *         cwd.mkpath( "." );
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     cwd.setPath("./log");
 *     cwd.makeAbsolute();
 *     if (!cwd.exists())
 *         cwd.mkpath(".");
 * 
 *     QDir::setCurrent(cwd.absolutePath());
 * 
 *     prm_filename += ".log";
 *     std::ofstream log_out_text(("./" + QString(prm_filename.c_str()).split("/").last()).toStdString().c_str());
 *     prm_handler.print_parameters (log_out_text,
 *                                   dealii::ParameterHandler::Text);
 * 
 *     QDir::setCurrent(this->params.run_dir.absolutePath());
 * }
 * 
 * 
 * 
@endcode
 <a name="plain-Functionrun"></a>
@code
 * void step42::MyFancySimulation::run()
 * {
 * 
 *     step42::CUDADriver testcase;
 * 
 *     testcase.gemm_tests();
 * 
 *     testcase.gemv_tests();
 * 
 *     testcase.complex_tests();
 * 
 *     testcase.feature_demonstration();
 * 
 *     std::cout << "Done." << std::endl;
 * }
 * 
 * 
 * 
@endcode
 <a name="plain-Functionmain"></a>
@code
 * int main(int argc, char *argv[])
 * {
 *     using namespace step42;
 * 
 *     int n_CUDA_devices;
 * 
 *     cudaGetDeviceCount(&n_CUDA_devices);
 *     std::cout
 *             << "N available CUDA devices : "
 *             <<  n_CUDA_devices
 *             << std::endl;
 * 
 *     int DevNo = 0;
 *     cudaSetDevice(DevNo);
 *     GPUInfo gpu_info(DevNo);
 * 
 *     gpu_info.get();
 * 
 *     MyFancySimulation machma(argc, argv, gpu_info);
 * 
 *     machma.run();
 * 
 * }
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDADriver_STEP_42_H
 * #define CUDADriver_STEP_42_H
 * 
 * #include <lac/release/blas_wrapper.hh>
 * #include <lac/release/cublas_wrapper.hh>
 * #include <base/CudaComplex.h>
 * 
 * namespace step42 {
 * 
@endcode
 <a name="plain-ClassCUDADriver"></a>
@code
 * class CUDADriver {
 * 
 *     typedef double Number;
 *     typedef SciPAL::CudaComplex<Number> cplxNumber;
 * 
 *     typedef cublas BW;
 * 
 * public:
 * 
 *     CUDADriver();
 * 
 *     ~CUDADriver() { BW::Shutdown(); }
 * 
 *     void gemm_tests();
 * 
 *     void gemv_tests();
 * 
 *     void complex_tests();
 * 
 *     void feature_demonstration();
 * 
 * private:
 * 
 * 
 * };
 * 
 * } // namespace step42 END
 * 
 * #endif // CUDADriver_STEP_42_H
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * #ifndef CUDA_DRIVER_STEP_42_HH
 * #define CUDA_DRIVER_STEP_42_HH
 * 
 * 
 * #include <cuda_driver_step-42.h>
 * #include <cuda_kernel_wrapper_step-42.cu.h>
 * 
 * 
 * 
 * #include <lac/development/cublas_Matrix.h>
 * 
 * 
 * #include <lac/development/cublas_Vector.h>
 * #include <lac/blas++.h>
 * #include <base/CudaComplex.h>
 * 
 * #include <cuda_runtime_api.h>
 * 
 * 
@endcode
 <a name="plain-ConstructorCUDADriver"></a>
@code
 * step42::CUDADriver::CUDADriver() {
 * 
 *     BW::Init();
 * }
 * 
 * 
 * 
 * 
@endcode
 <a name="plain-Functiongemm_tests"></a>
@code
 * void step42::CUDADriver::gemm_tests()
 * {
 * 
 *     const unsigned int n_rows = 4;
 *     const unsigned int n_cols = 3;
 *     const unsigned int n_elements = n_rows * n_cols;
 * 
 *     std::vector<Number>
 *             a(n_elements, 1.),
 *             b(n_elements, 2.),
 *             c(n_rows * n_rows, 1.23);
 * 
 *     for (unsigned int i = 0; i < b.size(); i++ )
 *         b[i] = i+1;
 * 
 * 
 *     SciPAL::Matrix<Number, BW>
 *             A(n_rows, n_cols, a),
 *             B(n_cols, n_rows, b),
 *             C(n_rows, n_rows, c);
 * 
 *      Number alpha = 1.1;
 *      Number beta = 2.;
 * 
 *      std::cout << "A : " << std::endl;
 *      A.print();
 * 
 *      std::cout << "B : " << std::endl;
 *      B.print();
 * 
 *      std::cout << "C : " << std::endl;
 *      C.print();
 * 
 * 
 * 
 *      std::cout << " ============ C = " << alpha << " * A ======" << std::endl;
 *      C = alpha * A; // * B + beta * C;
 *        std::cout << "C : " << std::endl;
 *      C.print();
 * 
 *      std::cout << " ============ C = A * B ======" << std::endl;
 *      C = A * B;   std::cout << "C : " << std::endl; C.print();
 * 
 *      std::cout << " ============ C = B * A ======" << std::endl;
 *      C = B * A;   std::cout << "C : " << std::endl; C.print();
 * 
 *      std::cout << " ============ C = alpha * A * B ======" << std::endl;
 *      C = alpha * A * B; std::cout << "C : " << std::endl; C.print();
 * 
 * 
 *      std::cout << " ============ C = " << alpha << " * A * B + "<<"C======" << std::endl;
 *      c.clear();
 *      c.resize(n_rows * n_rows, 1.);
 *      SciPAL::Matrix<Number, BW>
 *              D(n_rows, n_rows, c);
 * 
 * 
 * 
 *      C = D;
 *      std::cout << "C : " << std::endl; C.print();
 *      C = alpha * A * B + // beta *
 *              C; std::cout << "C : " << std::endl; C.print();
 * 
 * 
 *      std::cout << " ============ C = " << alpha << " * A * B + "  << beta  << " * C======" << std::endl;
 *      C = D;
 *      std::cout << "C : " << std::endl; C.print();
 * 
 *      C = alpha * A * B + beta * C;
 *      std::cout << "C : " << std::endl; C.print();
 * }
 * 
 * 
@endcode
 <a name="plain-Functiongemv_tests"></a>
@code
 * void step42::CUDADriver::gemv_tests()
 * {
 * #ifndef nUSE_ARRAY_EXPRESSIONS
 *      const unsigned int n_rows = 4;
 *      const unsigned int n_cols = 4;
 *      const unsigned int n_elements = n_rows * n_cols;
 * 
 *      Number alpha = 1.1;
 *      Number beta = 2.;
 * 
 *      std::vector<Number>
 *                 a(n_elements, 1.),
 *                 b(n_elements, 2.);
 * 
 * 
 *      for (unsigned int i = 0; i < a.size(); i++ )
 *          a[i] = i+1;
 * 
 *      SciPAL::Vector<Number, BW> vA, vB(n_elements), vC;
 *        vA = a;
 *        vB = SciPAL::Literal<Number>(2.3);
 *        vC = a;
 * 
 * 
 *        std::cout << "vA : " << std::endl;
 *        vA.print();
 * 
 *        std::cout << "vB : " << std::endl;
 *        vB.print();
 * 
 *        std::cout << " ============ vC = " << alpha << " * vA ======" << std::endl;
 *        vC = alpha * vA;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vC = " << alpha << " * vA + vB ======" << std::endl;
 *        vC = alpha * vA + vB;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vA = sin(vC) ======" << std::endl;
 *        const unsigned int n_sin_elements = n_elements;
 *        std::vector<Number> d(n_sin_elements);
 *        for(uint i = 0; i < d.size(); i++)
 *            d[i] = i* 2.* M_PI / d.size();
 * 
 *        SciPAL::Vector<Number, BW> vD; vD = d; //(n_sin_elements, 1, d);
 *       vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
 *          std::cout << "sin(vD) : " << std::endl;
 *        vD.print();
 * 
 * 
 *        SciPAL::Matrix<Number, BW>
 *                A(n_rows, n_cols, d);
 *        A = sin(A);
 *        A.print();
 * 
 * 
 *        std::cout << " ============ linear combination test ======" << std::endl;
 *        vC = 2.0 * vA;
 *        std::cout << "vC = 2.0 * vA" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + vB;
 *        std::cout << "vC = vA + vB" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + 2.0 * vB;
 *        std::cout << "vC = vA + 2.0 * vB" << std::endl;
 *        vC.print();
 * 
 *        vC = 2.0 * vA + 3.0 * vB;
 *        std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC =sin(2.0 * vA + 3.0 * vB);
 *        std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
 *        vC.print();
 * 
 *        vC = sqrt(vC);
 *        std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
 *        vC.print();
 * 
 *        vC = vA && vB;
 *        std::cout << "vC = vA .* vB" << std::endl;
 *        vC.print();
 * 
 *        vC = vA || vB;
 *        std::cout << "vC = vA ./ vB" << std::endl;
 *        vC.print();
 * 
 *        vC = (vA + vB) || (vA - vB);
 *        std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
 *        vC.print();
 * 
 * 
 * #endif
 * }
 * 
 * void step42::CUDADriver::complex_tests()
 * {
 *     std::cout<<"Entering tests for complex number array exptessions."<<std::endl;
 * #ifndef nUSE_ARRAY_EXPRESSIONS
 *      const unsigned int n_rows = 4;
 *      const unsigned int n_cols = 4;
 *      const unsigned int n_elements = n_rows * n_cols;
 * 
 *      SciPAL::CudaComplex<Number> alpha(1.1);
 *      SciPAL::CudaComplex<Number> beta (2.);
 * 
 *      std::vector<std::complex<Number> >
 *                 a(n_elements, 1.),
 *                 b(n_elements, std::complex<Number>(2., 2.0));
 * 
 * 
 *      for (unsigned int i = 0; i < a.size(); i++ )
 *          a[i] = std::complex<Number>(i+1, (i+1)/2.); //generate some inputs
 * 
 *        SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vA, vB, vC;
 *        vA = a;
 *        vB = b;
 *        vC = a;
 * 
 *        std::cout << "vA : " << std::endl;
 *        vA.print();
 * 
 *        std::cout << "vB : " << std::endl;
 *        vB.print();
 * 
 *        std::cout << " ============ vC = " << alpha.real() << " * vA ======" << std::endl;
 *        vC = alpha * vA;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vC = " << alpha.real() << " * vA + vB ======" << std::endl;
 *        vC = alpha * vA + vB;
 *          std::cout << "vC : " << std::endl;
 *        vC.print();
 * 
 *        std::cout << " ============ vA = sin(vC) ======" << std::endl;
 *        const unsigned int n_sin_elements = n_elements;
 *        std::vector<std::complex<Number> > d(n_sin_elements);
 *        for(uint i = 0; i < d.size(); i++)
 *            d[i] = std::complex<Number>(i* 2.* M_PI / d.size(), i* 4.* M_PI / d.size()) ;
 * 
 *        SciPAL::Vector<SciPAL::CudaComplex<Number>, BW> vD; vD = d; //(n_sin_elements, 1, d);
 *       vD = sin(vD); // For this to work the device-side apply() function has to be explicitly specialized.
 *          std::cout << "sin(vD) : " << std::endl;
 *        vD.print();
 * 
 *        std::cout << " ============ Matrix A = sin(A) ======" << std::endl;
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                A(n_rows, n_cols, d);
 *        A = sin(A);
 *        A.print();
 * 
 *        std::cout << " ============ Matrix B = sqrt(A) ======" << std::endl;
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                B(A);
 *        B = sqrt(A);
 *        B.print();
 * 
 *        std::cout << " ============ Matrix A = .exp(B) ======" << std::endl;
 *        A = exp(B);
 *        A.print();
 * 
 *        std::cout << " ============ Matrix A = B*C ======" << std::endl;
 *        SciPAL::Matrix<SciPAL::CudaComplex<Number>, BW>
 *                C(A);
 *        A = B*C;
 *        A.print();
 * 
 * 
 *        std::cout << " ============ linear combination test ======" << std::endl;
 * 
 *        vC = beta * vA; ;
 *        std::cout << "vC = 2.0 * vA" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + vB;
 *        std::cout << "vC = vA + vB" << std::endl;
 *        vC.print();
 * 
 * 
 *        vC = vA + beta * vB;
 *        std::cout << "vC = vA + 2.0 * vB" << std::endl;
 *        vC.print();
 * 
 *        vC = cplxNumber(2.0) * vA + cplxNumber(3.0) * vB;
 *        std::cout << "vC = 2.0 * vA + 3.0 * vB" << std::endl;
 *        vC.print();
 * 
 *        vC =sin(cplxNumber(2.0) * vA + cplxNumber(3.0) * vB);
 *        std::cout << "vC = sin(2.0 * vA + 3.0 * vB)" << std::endl;
 *        vC.print();
 * 
 *        vC = sqrt(vC);
 *        std::cout << "sqrt(sin(2.0 * vA + 3.0 * vB))" << std::endl;
 *        vC.print();
 * 
 *        vC = vA && vB;
 *        std::cout << "vC = vA .* vB" << std::endl;
 *        vC.print();
 * 
 *        vC = vA || vB;
 *        std::cout << "vC = vA ./ vB" << std::endl;
 *        vC.print();
 * 
 *        vC = (vA + vB) || (vA - vB);
 *        std::cout << "vC = (vA + vB) || (vA - vB)" << std::endl;
 *        vC.print();
 * 
 * 
 * #endif
 * }
 * 
 * 
 * void step42::CUDADriver::feature_demonstration()
 * {
 *     Number * bla = new Number[3];
 *     Number * bla2 = new Number[3];
 *     Number * bla3 = new Number[3];
 *     std::vector<Number*> h_testV(5);
 *     h_testV[0] = bla; std::cout<<"ptr1 " << bla << std::endl;
 *     h_testV[1] = bla2;std::cout<<"ptr2 " << bla2 << std::endl;
 *     h_testV[2] = bla3;std::cout<<"ptr3 " << bla3 << std::endl;
 * 
 *     SciPAL::Vector<Number*, cublas> d_testV(5);
 * 
 *  d_testV = h_testV;
 *  d_testV.print();
 * 
 * 
 * 
 * 
 * }
 * 
 * #endif // CUDA_DRIVER_STEP_42_HH
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * 
 * 
 * 
 * / *This file is part of SciPAL.
 * 
 *     SciPAL is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU Lesser General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     SciPAL is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU Lesser General Public License for more details.
 * 
 *     You should have received a copy of the GNU Lesser General Public License
 *     along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
 * * /
 * 
 * 
 * 
 * 
 * 
 @endcode
 */
