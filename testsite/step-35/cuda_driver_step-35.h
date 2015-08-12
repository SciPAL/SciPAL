//@sect3{File: cuda_driver_step-35.h}
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

Copyright  Lutz Künneke, Jan Lebert 2014
*/

#ifndef CUDADriver_STEP_35_H
#define CUDADriver_STEP_35_H

//std
#include<iostream>

//CUDA
#include <cuda_runtime.h>

//CUDA cufft
#include <cufft.h>

//SciPAL
#include <base/PrecisionTraits.h>
#include <lac/cublas_wrapper.hh>
#include <lac/cublas_Vector.h>
#include <lac/VectorCustomOperations.h>
#include <numerics/FFTWrapper.h>

//deal.II
#include <deal.II/lac/vector.h>

//Our stuff
//#include "ImageDecomposition.h"
#include "cuda_driver_step-35.hh"
#include "cuda_kernel_wrapper_step-35.cu.h"
#include "patch.h"
#include "cuda_helper.h"
#include "preprocessor_directives.h"
#include "ADMMParams.h"

#include <csignal>

namespace step35 {
// We encapsulate each project into a dedicated namespace
// in order to be able to re-use parts of a test program in others.

//Forward declarations
template<typename Mpatch,typename T,ParallelArch c> class CUDADriver;
template<typename Mpatch,typename T,ParallelArch c> struct constructStruct;
template<typename Mpatch,typename T,ParallelArch c> struct runStruct;

//@sect4{Struct: constructStruct}
//template specialized struct for the constructor and destructor of
//CUDADriver class based on which Mpatch we're dealing with

//Specialization for cset_dyadic
template<typename T,ParallelArch c>
struct constructStruct<cset_dyadic,T,c> {
    //@sect5{Function: constructFun}
    //dummy function needed by specialization, don't do anything
    //@param pt pointer to the friended CUDADriver instance
    // FIXME: in C++ we use references, not pointers.
    void constructFun(CUDADriver<cset_dyadic, T, c> *pt) {
        // FIXME: remove and comment out name of function argument. Then the compiler stops complaining as well.
        pt = pt; //supress compiler warning
    }

    //@sect5{Function: destructFun}
    //dummy function needed by specialization, don't do anything
    //@param pt pointer to the friended CUDADriver instance
     // FIXME: in C++ we use references, not pointers.
    void destructFun(CUDADriver<cset_dyadic, T, c> *pt) {
          // FIXME: remove and comment out name of function argument. Then the compiler stops complaining as well.
        pt = pt; //supress compiler warning
    }
};


//Specialization for cluster<cset_small,T>
template<typename T,ParallelArch c>
struct constructStruct<cluster<cset_small,T>,T,c> {
    //@sect5{Function: constructFun}
    //Callen when CUDADriver is constructed, start worker threads which do the real work
    //@param pt pointer to the friended CUDADriver instance
     // FIXME: in C++ we use references, not pointers.
    void constructFun(CUDADriver<cluster<cset_small,T>,T,c> *pt) {
        if (pt->inf->ext_depth > 1) {//TODO 3d
            std::cerr << "Exact dykstra with 3d images is not yet implemented!" << std::endl;
            std::abort();
        }
        // FIXME: C++ -> std::cout
        printf("CUDA init...\n");

        //set up worker pool
        // FIXME: pointer-free solution. Read Scott Meyer's effective C++ series.
        pt->thread_handles = new std::thread[pt->stream_count];

        //start up threads
        for (int i = 0; i < pt->stream_count; i++) {
            pt->thread_handles[i] = std::thread(&CUDADriver<cluster<cset_small,T>,T,c>::stream_handler, pt, i);
        }
    }

    //@sect5{Function: destructFun}
    //Called when CUDADriver is destructed, kill the stuff we've started in constructFun()
    //@param pt pointer to the friended CUDADriver instance
    void destructFun(CUDADriver<cluster<cset_small,T>,T,c> *pt) {
        pt->queue->add(pt->cend);
        //waits for all threads to end
        for (int i = 0; i < pt->stream_count; i++) {
            if (pt->thread_handles[i].joinable()) {
                pt->thread_handles[i].join();
            }
        }
        delete[] pt->thread_handles;
    }
};

//@sect4{Struct: runStruct}
//template specialized struct for the iterate function of
//CUDADriver class based on which Mpatch we're dealing with
template<typename T,ParallelArch c>
//Specialization for cluster<cset_small, T>
struct runStruct<cluster<cset_small, T>, T, c > {
    //@sect5{Function: run_func}
    //Run function of CUDADriver for cluster<cset_small, T>,
    //adds cluster<cset_small, T> items to queue
    //@param pt pointer to the friended CUDADriver instance
    void run_func(CUDADriver<cluster<cset_small,T>,T, c> *pt) {
        //First element of linked list
        cluster<cset_small,T>* cpoint;
        cpoint=pt->croot;

        //Reset all $q$ values to zero
        while ( cpoint->next != NULL ) {
            cpoint->reset();
            cpoint=cpoint->next;
        }

        //The number of Iterations is limited, this is to be replaced by a convergence check in the future
        for (int rep=0; rep<3; rep++) {
            //Start at first element of linked list
            cpoint=pt->croot;
            //Go through linked list to it's end
            while ( cpoint != NULL ) {
                //Add element to queue
                pt->queue->blocking_add(cpoint);
                //Go to next element
                cpoint = cpoint->next;
            }
        }
    }
};

template<typename T,ParallelArch c>
//Specialization for cset_dyadic
struct runStruct<cset_dyadic,T,c > {
    //@sect5{Function: run_func}
    //@brief Run function of CUDADriver for cset_dyadic, only calls nostream_handler()
    //@param pt pointer to the friended CUDADriver instance
    void run_func(CUDADriver<cset_dyadic,T, c> *pt) {
        pt->nostream_handler2();
    }
};

// @sect4{Class: CUDADriver}
//
// Main class, sets up an instance of queue and info class. A number of threads is started which
// each handles a cuda stream and processes items obtained by the queue. The main thread adds items
// to the queue until everything is processed.
template<typename Mpatch,typename Mdouble, ParallelArch c>
class CUDADriver {
    friend  struct runStruct<Mpatch,Mdouble,c>;
  public:
    //complex type: double2 for double, float2 for float
    typedef typename PrecisionTraits<Mdouble, c>::ComplexType complex;

    //Thread-safe queue used to distribute work on worker threads
    wqueue<Mpatch, Mdouble, c> *queue;

    //Info class instance which contains general information about the image which is currently processed
    ImageInfo<Mdouble, complex, c> *inf;

    //Pointer to root and current element of patch cluster list
    Mpatch* croot,*cpoint;

    //Element used to signal all worker threads to shut down
    Mpatch* cend;

    //Temporary fields that need not to be known outside of the driver
    // FIXME: raw pointers are bad practice, source for constant trouble (e.g. memory leaks) and anyway error-prone.
    complex *fm1,*fm1_h;
    //Mdouble *tmp_h; // ,*tmp2_d, *tmp_d,*lag1,*lag2,*tmp_haar2,*tmp_lagr,tmp_haar
    dealii::Vector<Mdouble> tmp_h;

    // FIXME: For the device side arrays use this:
    SciPAL::Vector<Mdouble, cublas> tmp_d, tmp2_d,lag1,lag2,tmp_haar2,tmp_lagr,tmp_haar;


    cufftHandle *plan_fft,*iplan_fft;

    //Number of CUDA streams and worker threads to create
    int stream_count;
    //Worker thread handles
    std::thread* thread_handles;

    //Contains dimensions of rectangles
    std::vector<int> dimRect;

#ifdef TIME_KERNELS
    // Total GPU time for dykstra kernels in ms
    float kernel_time;
    // Output file where kernel execution time is written at shutdown
    std::ofstream out_time;
#endif

    // @sect5{Constructor: CUDADriver}
    // @param mcroot start object of the linked list
    // @param im image
    // @param fpsf psf function
    // @param cs_h array of $c_s$
    // @param nx width of the image
    // @param ny height of the image
    // @param n maximum frame size
    // @param gamma regularization parameter w.r.t. x
    // @param sigma estimated gauss variance
    // @param newreg regularisation Method
    // @param dim dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    // FIXME: cs_h must some a reference to some decent type not a raw pointer
    CUDADriver(Mpatch* mcroot, std::vector<Mdouble> &input_image, std::vector<Mdouble> &fpsf, Mdouble* cs_h,
               // FIXME: use ImageStack class developed for SOFI and ask that for nx, ny, nz, if needed
               const int nx, const int ny, const int nz, const int n, Mdouble gamma, int sigma,
               regularisationType regType, const int dim)
        :
          // FIXED: use RAII
          inf(new ImageInfo<Mdouble, complex, c>(input_image, fpsf, cs_h, nx, ny, nz, n,
                                            gamma, sigma, regType, dim)),
    // FIXME: more vector instantiations go here
          tmp_d(inf->ext_num_pix),tmp2_d(inf->ext_num_pix),lag1(inf->ext_num_pix),lag2(inf->ext_num_pix),
          tmp_haar2(inf->ext_num_pix),tmp_haar(inf->ext_num_pix),tmp_lagr(inf->ext_num_pix),
          tmp_h(inf->ext_num_pix),
          dimRect(11)
    {
        getLastCudaError("CUDA in error state while driver init\n");
        //Number of CUDA streams (and thus std::threads) to use, 5 seems to be
        //reasonable as we don't want worker threads to fight for work load and
        //the PCIe bus is the bottelneck to the GPU
        stream_count = 5;

        //Root element of linked list
        croot=mcroot;

        //Wait for async memcpys
        cudaDeviceSynchronize();

        //Calculate reasonable maximum queue length
        unsigned int queue_len=(unsigned int) (ny/(n))*(nx/(n));
        //Set up queue
        queue=new wqueue<Mpatch,Mdouble, c>(queue_len);

        //The info class holds all real variables

        //Assign dimensions of rectangles for decomposition
        for (int i = 0; i < 11; i++)
            dimRect[i] = 1 << i;

        //Cufft setup
        int fm1_size=inf->ext_depth*inf->ext_width*(inf->ext_height/2+1)*sizeof(complex);
        fm1_h=new complex[inf->ext_depth*inf->ext_width*(inf->ext_height/2+1)];
        //Create cufft plans
        plan_fft=new cufftHandle();
        iplan_fft=new cufftHandle();
#ifdef DOUBLE_PRECISION
        cufftPlan3d(plan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_D2Z);
        cufftPlan3d(iplan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_Z2D);
#else
        cufftPlan3d(plan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_R2C);
        cufftPlan3d(iplan_fft, inf->ext_depth, inf->ext_width, inf->ext_height,  CUFFT_C2R);
#endif

        //Allocate space on device
        // FIXME: device mem alloc has to be done in the CTor call of the SciPAL::Vector objects and this has to happen in the initializer list of the CTor of this class as is good an common practice in C++ to guarantee a defined state of a class attribute right from the start.
        checkCudaErrors(cudaMalloc((void **)&fm1, fm1_size));
        // checkCudaErrors(cudaMalloc((void **)&tmp_d, inf->n_bytes_per_frame));
        // FIXME: why is there host allocation when device arrays are alloc'd?
        //tmp_h=new Mdouble[inf->ext_num_pix];
        //checkCudaErrors(cudaMalloc((void **)&tmp2_d, inf->n_bytes_per_frame));
        //checkCudaErrors(cudaMalloc((void **)&lag1, inf->n_bytes_per_frame));
        //checkCudaErrors(cudaMalloc((void **)&lag2, inf->n_bytes_per_frame));
        //checkCudaErrors(cudaMalloc((void **)&tmp_haar, inf->nx2*inf->ny2*sizeof(Mdouble))); //TODO 3d
        //checkCudaErrors(cudaMalloc((void **)&tmp_lagr, inf->nx2*inf->ny2*sizeof(Mdouble))); //TODO 3d
        //checkCudaErrors(cudaMalloc((void **)&tmp_haar2, inf->nx2*inf->ny2*sizeof(Mdouble))); //TODO 3d

        //Init the lagrangian fields
        //step35::Kernels<Mdouble> kernel;
        //kernel.reset(lag1, inf->ext_num_pix);
        //kernel.reset(lag2, inf->ext_num_pix);
        //Generate Mpatch element used to signal threads to shut down
        cend=new Mpatch(0);

#ifdef TIME_KERNELS
        kernel_time = 0;
        //Write total dykstra kernel execution time at shutdown to file "kernel_time"
        out_time.open("kernel_time");
#endif

        //Create specialized constructStruct which does specialized setup based
        //on which Mpatch we're dealing with
        constructStruct<Mpatch, Mdouble, c > rs;
        rs.constructFun(this);
    }

    //@sect5{Destructor: CUDADriver}
    //Signal all threads to end, calls queue destructor
    ~CUDADriver () {

        //Create specialized constructStruct which does specialized destruction based
        //on which Mpatch we're dealing with
        constructStruct<Mpatch, Mdouble, c > rs;
        rs.destructFun(this);

        //Kill queue, it should be empty now, all worker threads are dead
        delete queue;

#ifdef TIME_KERNELS
        //Write total dykstra kernel execution time at shutdown to file "kernel_time"
        if (kernel_time)
            out_time << kernel_time << std::endl;
        out_time.close();
#endif

        //Cleanup
        delete[] fm1_h;
        //delete[] tmp_h;
        cufftDestroy(*plan_fft);
        cufftDestroy(*iplan_fft);
        delete plan_fft;
        delete iplan_fft;
        delete cend;
        delete inf;
    }

    //@sect5{Function: nostream_handler}
    //@brief Wrapper if we use the dyadyc_dykstra Kernel.
    void nostream_handler() {
#ifdef TIME_KERNELS
        cudaEvent_t start, stop;
        float time;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        cudaEventSynchronize(stop);
#endif

        step35::Kernels<Mdouble> kernels;

        if (inf->dim == 2) {
            //Do the approximate method shifted by 2^so. This choice of shifts enables us to omit some sets in later instances of the kernel
            for (int so=0; so<5; so++) {
                for (int z=0; z<inf->ext_depth; z++) {
                    kernels.dyadic_dykstra(inf->e_d,inf->ext_width, inf->ext_height, inf->ext_depth, 1 << so, 1 << so, z, so);
                }
            }
        }
        else { //TODO 3d
            for (int so=0; so<5; so++) {
                for (int z=0; z<inf->ext_depth; z++) {
                    kernels.dyadic_dykstra(inf->e_d.array().val(), inf->ext_width, inf->ext_height, inf->ext_depth, 1 << so, 1 << so, z, so);
                }
            }
        }

#ifdef TIME_KERNELS
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time, start, stop);
        kernel_time += time;
#endif
    }

    //@sect5{Function: nostream_handler2}
    //@brief Wrapper if we use the stoch_dykstra Kernel.
    void nostream_handler2() {
        int randomPower = 5;//rand() % 11;

        step35::Kernels<Mdouble> kernels;
            std::cout << "random number:" << randomPower << std::endl;

        if (inf->dim == 2) {
//            for(int i = 50000;i<1024+50000;i++)
//                std::cout << inf->e_d(i)<< " " ;
            //Do the approximate method shifted by 2^so. This choice of shifts enables us to omit some sets in later instances of the kernel
            for (int so=-1; so<5; so++) {
                for (int z=0; z<inf->ext_depth; z++) {
                    kernels.stoch_dykstra(inf->e_d,inf->ext_height, inf->ext_width, inf->ext_depth, (so<0)? 0: (1 << so), (so<0)? 0: (1 << so), z, so, randomPower);
                }
            }
        }
        else { //TODO 3d
            for (int so=0; so<5; so++) {
                for (int z=0; z<inf->ext_depth; z++) {
                    kernels.dyadic_dykstra(inf->e_d.array().val(), inf->ext_width, inf->ext_height, inf->ext_depth, 1 << so, 1 << so, z, so);
                }
            }
        }
    }

    //@sect5{Function: stream_handler}
    //@brief Handles one cuda stream if we use cluster<cset_small>, queue worker thread.
    //This function is executed by `stream_count` threads in parallel, each creates an cuda
    //stream and coordinates the work done on that stream. It removes items from the queue
    //and processes. End thread when terminate item is recieved from the queue.
    //@param rank dummy threadid used for debugging
    void stream_handler(int rank) {
        printf("Hello from worker thread %ld \n", (long int)rank); //for debug

        //CUDA device id has to be set for each thread
        checkCudaErrors(cudaSetDevice(inf->device_id));
        checkCudaErrors(cudaDeviceSynchronize());

        //Create CUDA stream
        cudaStream_t mystream;
        checkCudaErrors(cudaStreamCreate(&mystream));

        std::list<Mpatch*> cluster_vec;
        //Maximum number of elements we try to get from queue in one call
        const int max_num_cluster = 128;
        //Minimum number of elements we try to get from queue in one call, not guaranteed
        const int min_num_cluster = 12;

        //Allocate the maximum needed memory on host and device for
        //cluster information and q_offset

        //Array of device pointers which contains the information of the clusters
        int **cluster_info; //=new int*[4*max_num_cluster];
        int **cluster_info_d;
        checkCudaErrors(cudaHostAlloc((void **)&cluster_info, 4*max_num_cluster*sizeof(int*), cudaHostAllocWriteCombined));
        checkCudaErrors(cudaMalloc((void **)&cluster_info_d, 4*max_num_cluster*sizeof(int*)));

        //Where to find the $q$ values for eauch cluster
        int *q_offset; //=new int[max_num_cluster];
        int *q_offset_d;
        checkCudaErrors(cudaHostAlloc((void **)&q_offset, max_num_cluster*sizeof(int), cudaHostAllocWriteCombined));
        checkCudaErrors(cudaMalloc((void **)&q_offset_d, max_num_cluster*sizeof(int)));

        //All $q$ values
        Mdouble* d_q;
        checkCudaErrors(cudaMalloc((void **)&d_q, 1024*max_num_cluster*sizeof(Mdouble)));

        step35::Kernels<Mdouble> kernels;

#ifdef TIME_KERNELS
        //Collect GPU dykstra kernel total execution time
        cudaEvent_t start, stop;
        float time, total_time;
#endif

        //Get a list of min_num_cluster to max_num_cluster clusters
        cluster_vec = queue->get(min_num_cluster, max_num_cluster);
        //queue->get(min, max) has a timeout, check if we actually received something
        while ( !cluster_vec.size() ) {
            cluster_vec = queue->get(min_num_cluster, max_num_cluster);
        }

        //Big while loop until terminate signal is received
        while ( cluster_vec.front()->size != 0 ) {
#ifdef TIME_KERNELS
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
#endif
            //Number of elements we received
            int num_of_cluster=cluster_vec.size();

            //In the big $q$ array, where can I find the $q$ for the first pixel of my cluster?
            int offset=0;

            //Only used in the following for loop
            int k=0;
            //Fill array of device pointers
            for (auto it=cluster_vec.begin(), end=cluster_vec.end(); it!=end; ++it, k++) {
                //In the big $q$ array, where can I find the $q$ value for the first pixel of my cluster?
                q_offset[k]=offset;
                //Create an array of device pointers \n
                //First device pointer points to the array containing the edge lengths of the cset_small
                cluster_info[4*k]=(*it)->nvec_d;
                //Second the pointer to the horizontal start pixel positions
                cluster_info[4*k+1]=(*it)->ivec_d;
                //Third the pointer to the vertical start pixel positions
                cluster_info[4*k+2]=(*it)->jvec_d;
                //Fourth the pointer to the number of cset_small in one cluster
                cluster_info[4*k+3]=(*it)->pnum_d;

                offset+=(*it)->size;
            }

            //Copy $q$ for all frames to device
            offset=0;
            for (auto it=cluster_vec.begin(), end=cluster_vec.end(); it!=end; ++it) {
                checkCudaErrors(cudaMemcpyAsync(&(d_q[offset]), (*it)->qmat, (*it)->size*sizeof(Mdouble), cudaMemcpyHostToDevice, mystream));
                offset+=(*it)->size;
            }

            //Copy cluster_info to device
            checkCudaErrors(cudaMemcpyAsync(cluster_info_d, cluster_info, 4*num_of_cluster*sizeof(int*), cudaMemcpyHostToDevice, mystream));
            //Copy q_offset to device
            checkCudaErrors(cudaMemcpyAsync(q_offset_d, q_offset, num_of_cluster*sizeof(int), cudaMemcpyHostToDevice, mystream));

#ifdef TIME_KERNELS
            cudaEventRecord(start, mystream);
#endif
            //Calculate Dykstra's algorithm on all frames in the list of clusters we recieved
            kernels.dykstra(d_q,inf->e_d.array().val(), cluster_info_d, q_offset_d, num_of_cluster, inf->width, inf->height, 1024, &mystream);

#ifdef TIME_KERNELS
            cudaEventRecord(stop, mystream);
#endif

            //Copy $q$ back
            offset=0;
            for (typename std::list<Mpatch*>::iterator it=cluster_vec.begin(), end=cluster_vec.end(); it!=end; ++it) {
                checkCudaErrors(cudaMemcpyAsync( (*it)->qmat, &(d_q[offset]), (*it)->size*sizeof(Mdouble), cudaMemcpyDeviceToHost, mystream));
                offset+=(*it)->size;
            }

            //Blocks until stream has completed all operations
            cudaStreamSynchronize(mystream);
            getLastCudaError("__dysktra<<<>>> execution failed\n");

#ifdef TIME_KERNELS
            cudaEventElapsedTime(&time, start, stop);
            total_time += time;
#endif
            //Signal queue that task is done
            queue->task_done(cluster_vec);

            //Get a list of min_num_cluster to max_num_cluster clusters
            cluster_vec = queue->get(min_num_cluster, max_num_cluster);
            while ( !cluster_vec.size() ) {
                cluster_vec = queue->get(min_num_cluster, max_num_cluster);
            }
        } //End big while loop

        //Worker thread shutdown cleanup
        std::cout << "All done! Goodbye from thread " << rank << std::endl;

#ifdef TIME_KERNELS
        out_time << total_time << std::endl;
#endif

        //Insert terminate signal back to queue
        queue->add(cend);
        //Destroy stream
        cudaStreamDestroy(mystream);
        //Free allocated memory
        cudaFree(d_q);
        cudaFree(cluster_info_d);
        cudaFree(q_offset_d);
        cudaFreeHost(cluster_info);
        cudaFreeHost(q_offset);
    }

    //@sect5{Function: conv2}
    //@brief Cyclic convolution, based on FFT
    //@param in input array
    //@param out output array
      void conv2(Mdouble *in, Mdouble *out) {

#ifdef DOUBLE_PRECISION
        cufftExecD2Z(*plan_fft, in, fm1);
#else
        cufftExecR2C(*plan_fft, in, fm1);
#endif
        checkCudaErrors(cudaDeviceSynchronize());

        //Convolve, multiply in Fourier space
        step35::Kernels<Mdouble> kernel;
        kernel.element_norm_product(fm1, inf->fpsf_d, inf->ext_width, inf->ext_height, inf->ext_depth);

        //Transform back
        // FIXME: replace by SciPAL's FFT wrappers
#ifdef DOUBLE_PRECISION
        cufftExecZ2D(*iplan_fft, fm1, out);
#else
        cufftExecC2R(*iplan_fft, fm1, out);
#endif
        checkCudaErrors(cudaDeviceSynchronize());
        getLastCudaError("cufft error!\n");
    }

    //@sect5{Function: conv2}
    //@brief Cyclic convolution, based on FFT
    //@param in input array
    //@param out output array

    void conv2_ET(SciPAL::Vector<Mdouble,cublas> &in, SciPAL::Vector<Mdouble,cublas> &out) {
       //SciPAL FFT
        SciPAL::Vector<SciPAL::CudaComplex<Mdouble>,cublas> tmpfft_d(inf->ext_num_pix);
       SciPAL::CUDAFFT<Mdouble,2,SciPAL::TransformType<Mdouble>::FFTType_R2C,gpu_cuda>
               cuda_fft(inf->ext_height,inf->ext_width,tmpfft_d,in);

        SciPAL::Vector<SciPAL::CudaComplex<Mdouble>,cublas> tmpifft_d(inf->ext_num_pix);
       SciPAL::CUDAFFT<Mdouble,2,SciPAL::TransformType<Mdouble>::FFTType_C2C,gpu_cuda>
               cuda_ifft(inf->ext_height,inf->ext_width,tmpifft_d/*out*/,tmpfft_d);

        cuda_fft(tmpfft_d,in,FORWARD);
        //Convolve, multiply in Fourier space
        step35::Kernels<Mdouble> kernel;
        kernel.element_norm_product(tmpfft_d.array().val(), inf->fpsf_d, inf->ext_width, inf->ext_height, inf->ext_depth);

        cuda_ifft(tmpifft_d/*out*/,tmpfft_d,BACKWARD);

        kernel.real(out, tmpifft_d);
    }

    //@sect5{Function: projection_gauss}
    //@brief CPU implementation of dykstras algorithm with Gauss projection, performs projection for a single set
    //@param d_q pointer to vector q
    //@param im_h pointer to image
    //@param cs_h array of weights
    //@param width width of the image
    //@param height height of the image
    //@param n edge length of the set
    //@param i0 starting row of the set
    //@param j0 starting col of the set
    void projection_gauss(Mdouble* d_q, Mdouble* im_h, Mdouble* cs_h, const int width, const int height,
                          const int n, const int i0, const int j0) {
        //Do the difference and square
        int size=n*n;
        std::vector<Mdouble> f(size);
        Mdouble square_sum;
        int is,js;
        //Square the vector
        for (int i=0; i<size; i++) {
            is = i/n + i0;
            if ( is >= width )
                std::cout << "Error in gauss Projection, dimension too large " << std::endl;
            js = i%n + j0;
            f[i]=im_h[is*height+js]-d_q[i];
            d_q[i]=-f[i];
            square_sum=square_sum+f[i]*f[i];
        }
        square_sum*=cs_h[n];
        //Check the condition
        if ( square_sum > 1.0 ) {
            square_sum=1.0/sqrt(square_sum);
            for (int i=0; i<size; i++) {
                //Do the projection
                is = i/n + i0;
                js = i%n + j0;
                im_h[is*height+js]=f[i]*square_sum;
                d_q[i]+=f[i];
            }
        } else {
            for (int i=0; i<size; i++) {
                is = i/n + i0;
                js = i%n + j0;
                im_h[is*height+js]=f[i];
                d_q[i]+=f[i];
            }
        }
    }

    //@sect5{Function: iterate}
    //@brief Initiate Dykstra's algothrim once for all frames
    void iterate() {
        runStruct<Mpatch, Mdouble, c > rs;
        //The specialized runStruct does the actuall work
        rs.run_func(this);
    }

    //@sect5{Function: push_data}
    //@brief Push all data from host to device
    void push_data() {
        inf->im_d = inf->im_h;
        inf->x_d = inf->x_h;
        inf->z_d = inf->z_h;
        inf->e_d = inf->e_h;
        tmp_d = tmp_h;

//        checkCudaErrors(cudaMemcpy(inf->im_d.array().val(), &(inf->im_h[0]), inf->n_bytes_per_frame, cudaMemcpyHostToDevice));
//        checkCudaErrors(cudaMemcpy(inf->x_d.array().val(), &(inf->x_h[0]), inf->n_bytes_per_frame, cudaMemcpyHostToDevice));
//        checkCudaErrors(cudaMemcpy(inf->z_d.array().val(), &(inf->z_h[0]), inf->n_bytes_per_frame, cudaMemcpyHostToDevice));
//        checkCudaErrors(cudaMemcpy(inf->e_d.array().val(), &(inf->e_h[0]), inf->n_bytes_per_frame, cudaMemcpyHostToDevice));
//        checkCudaErrors(cudaMemcpy(tmp_d.array().val(), &(tmp_h[0]), inf->n_bytes_per_frame, cudaMemcpyHostToDevice));
    }

    //@sect5{Function: get_data}
    //@brief Pull all data from device to host
    void get_data() {
        inf->im_d.push_to(inf->im_h);
        inf->x_d.push_to(inf->x_h);
        inf->z_d.push_to(inf->z_h);
        inf->e_d.push_to(inf->e_h);
        tmp_d.push_to(tmp_h);
//        checkCudaErrors(cudaMemcpy(&(inf->im_h[0]), inf->im_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToHost));
//        checkCudaErrors(cudaMemcpy(&(inf->x_h[0]), inf->x_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToHost));
//        checkCudaErrors(cudaMemcpy(&(inf->z_h[0]), inf->z_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToHost));
//        checkCudaErrors(cudaMemcpy(&(inf->e_h[0]), inf->e_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToHost));
//        checkCudaErrors(cudaMemcpy(&(tmp_h[0]), tmp_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToHost));
    }

    // FIXME: use SciPAl vectors and ETs!!!
    //@sect5{Function: x_step}
    //@brief Performs approximative argmin with respect to $x$
    void x_step(const Mdouble rho1,const Mdouble rho2) {
        step35::Kernels<Mdouble> kernel;
        //$\text{tmp}_d=I$tmp2
        kernel.reset(tmp_d.array().val(), inf->ext_num_pix);
        kernel.sum(tmp_d.array().val(), inf->im_d.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        //$\text{tmp}_d=I-e$
        kernel.diff(tmp_d.array().val(), inf->e_d.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        conv2(inf->x_d.array().val(), tmp2_d.array().val());
        //$\text{tmp}_d=I-e-A*x$
        kernel.diff(tmp_d.array().val(), tmp2_d.array().val(), inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth);



        //$\text{tmp}_d=\left(I-e-A*x\right)*\rho_1$
        kernel.mult(tmp_d.array().val(), rho1, inf->ext_num_pix);
        //$\text{tmp}_d=\left((im-e-A*x\right)*\rho_1+\Upsilon_1$
        kernel.sum(tmp_d.array().val(), lag1.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);


        //$\text{tmp}_d=A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
        conv2(tmp_d.array().val(), tmp_d.array().val());

        //$\text{tmp2}_d=z$
        kernel.reset(tmp2_d.array().val(), inf->ext_num_pix);
        kernel.sum(tmp2_d.array().val(), inf->z_d.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        //$\text{tmp2}_d=z-x$ 

        kernel.diff(tmp2_d.array().val(), inf->x_d.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);
        //$\text{tmp2}_d=z-x$


        //$\text{tmp2}_d=\left((z-x\right)*\rho_2$
        kernel.mult(tmp2_d.array().val(), rho2, inf->ext_num_pix);


        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
        kernel.sum(tmp2_d.array().val(), tmp_d.array().val(), inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2$
        kernel.diff(tmp2_d.array().val(), lag2.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp2}_d=\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        kernel.mult(tmp2_d.array().val(), inf->gamma, inf->ext_num_pix);

        //$x=x+\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        kernel.sum(inf->x_d.array().val(), tmp2_d.array().val(), 0, inf->ext_width, inf->ext_height, inf->ext_depth);


    }
    // FIXME: use SciPAl vectors and ETs!!!
    //@sect5{Function: x_step}
    //@brief Performs approximative argmin with respect to $x$
    void x_step_ET(const Mdouble rho1,const Mdouble rho2) {
        step35::Kernels<Mdouble> kernel;

        //$\text{tmp}_d=I-e$
        tmp_d = inf->im_d - inf->e_d;

        conv2(inf->x_d.array().val(), tmp2_d.array().val());
        //$\text{tmp}_d=I-e-A*x$
        kernel.diff(tmp_d.array().val(), tmp2_d.array().val(), inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp}_d=\left((I-e-A*x\right)*\rho_1+\Upsilon_1$
        tmp_d = lag1 + (rho1*tmp_d); //but tmp_d = (rho1*tmp_d)+lag1: Does NOT work!!!

        //$\text{tmp}_d=A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
        conv2_ET(tmp_d, tmp_d);

        //$\text{tmp2}_d=z-x$
        tmp2_d = inf->z_d - inf->x_d;

        //$\text{tmp2}_d=\left((z-x\right)*\rho_2$
        tmp2_d *= rho2;

        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)$
        kernel.sum(tmp2_d.array().val(), tmp_d.array().val(), inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth);

        //$\text{tmp2}_d=\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2$
        tmp2_d -= lag2;

        //$\text{tmp2}_d=\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        tmp2_d *= inf->gamma;

        //$x=x+\left(\left(z-x\right)*\rho_2+A*\left(\left(I-e-A*x\right)*\rho_1+\Upsilon_1\right)-\Upsilon_2\right)*\gamma$
        inf->x_d += tmp2_d;

    }

    //@sect5{Function: update_lagrangian}
    //@brief Update the lagrangian estimates
    //$ \Upsilon_1^{r+1} = \Upsilon_1^r + \alpha_1 \left( I - A*x - \epsilon \right) $ \n
    //$\Upsilon_2^{r+1} = \Upsilon_2^r + \alpha_2 \left( x - z \right) $
    void update_lagrangian(Mdouble alpha1, Mdouble alpha2) {
        step35::Kernels<Mdouble> kernel;
        kernel.reset(tmp_d.array().val(), inf->ext_num_pix);
        kernel.sum(tmp_d.array().val(), inf->x_d.array().val(),0, inf->ext_width, inf->ext_height, inf->ext_depth);
        conv2(tmp_d.array().val(), tmp_d.array().val());
        //Update the lagrangian estimates
        kernel.update_lagrangian(lag1.array().val(), lag2.array().val(), inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth,
                                 alpha1, alpha2, inf->e_d.array().val(), inf->im_d.array().val(),tmp_d.array().val(), inf->x_d.array().val(), inf->z_d.array().val());
    }

    //@sect5{Function: dykstra_gauss}
    //@brief Wrapper around the whole projection procedure
    //@param rho parameter for the first constraint
    void dykstra_gauss(const Mdouble rho) {
        conv2(inf->x_d.array().val(),tmp_d.array().val());
        step35::Kernels<Mdouble> kernel;
        //the value to be projected is copied into e_d
        kernel.prepare_e(inf->e_d.array().val(), inf->im_d.array().val(), tmp_d.array().val(), lag1.array().val(), rho, inf->sigma, inf->ext_width, inf->ext_height, inf->ext_depth);
        //Perform the Dykstra Algotithm
        iterate();
    }

    //@sect5{Function: m_smoothing}
    //@brief Smoothing step
    //@param gamma Strength of the regularization, small choices render the algorithm unstable
    //TODO rho2
    void m_smoothing(Mdouble gamma, Mdouble rho2) {
        step35::Kernels<Mdouble> kernel;
        //Regularization by Haar Wavelet Coefficient sparsity
        if ( inf->regType == haar ) {
            if (inf->ext_depth > 1) {//TODO 3d
                std::cerr << "Haar regularisation with 3d images is not yet implemented!" << std::endl;
                std::abort();
            }
            //kernel.reset(tmp_haar,inf->nx2*inf->ny2);
            tmp_haar = SciPAL::Vector<Mdouble,cublas>(inf->nx2*inf->ny2);
            tmp_lagr = SciPAL::Vector<Mdouble,cublas>(inf->nx2*inf->ny2);

            //kernel.reset(tmp_lagr,inf->nx2*inf->ny2);

            //Copy $x$ into the bigger temp variable while conserving its shape
            for (int i=0; i<inf->nx2; i++) {
                if ( i < inf->ext_height) {
                    checkCudaErrors(cudaMemcpyAsync(&(tmp_haar.array().val()[i*inf->ny2]), &(inf->x_d.array().val()[i*inf->ext_height]),
                                                    inf->ext_width*sizeof(Mdouble), cudaMemcpyDeviceToDevice));
                    checkCudaErrors(cudaMemcpyAsync(&(tmp_lagr.array().val()[i*inf->ny2]), &(lag2.array().val()[i*inf->ext_height]),
                                                    inf->ext_width*sizeof(Mdouble), cudaMemcpyDeviceToDevice));
                }
            }
            checkCudaErrors(cudaDeviceSynchronize());

            //Forward 2D Haar Wavelet transform
            kernel.haar(tmp_haar.array().val(),tmp_haar2.array().val(),inf->ny2);
            kernel.haar(tmp_lagr.array().val(),tmp_haar2.array().val(),inf->ny2);
            kernel.soft_threshold(tmp_haar.array().val(),lag2.array().val(),tmp_haar.array().val(),rho2,gamma,inf->nx2*inf->ny2);
            //Backward 2D Haar Wavelet transform
            kernel.inverse_haar(tmp_haar.array().val(),tmp_haar2.array().val(),inf->ny2);

            //Copy back, pay attention not to mess up the shape
            for (int i=0; i<inf->nx2; i++) {
                if ( i < inf->ext_width )
                    checkCudaErrors(cudaMemcpyAsync(&(inf->z_d.array().val()[i*inf->ext_height]), &(tmp_haar.array().val()[i*inf->ny2]),
                                                    inf->ext_width*sizeof(Mdouble), cudaMemcpyDeviceToDevice));
            }
            checkCudaErrors(cudaDeviceSynchronize());
        }

        //Regularization by direct space sparsity
        if ( inf->regType == sparse ) {
            //checkCudaErrors(cudaMemcpy(inf->z_d, inf->x_d, inf->framesize, cudaMemcpyDeviceToDevice));
            kernel.soft_threshold(inf->z_d.array().val(),lag2.array().val(),inf->x_d.array().val(), rho2, gamma, inf->ext_num_pix);
            //kernel.tv_regularization(inf->x_d,inf->z_d,lag2,gamma,rho2,inf->ext_width,inf->ext_height,inf->ext_depth);
        }
        //Regularization by Fourier Space L_2 Norm
        if ( inf->regType == quadratic ) {
            checkCudaErrors(cudaMemcpy(inf->z_d.array().val(), inf->x_d.array().val(), inf->n_bytes_per_frame, cudaMemcpyDeviceToDevice));
            //the solution with smallest L_2 Norm is obtained, this corresponds to the pseudoinverse
            kernel.pseudo_inverse(inf->z_d.array().val(),lag2.array().val(),rho2,gamma,inf->ext_num_pix);
        }
    }
};
} // namespace step35 END

#endif // CUDADriver_STEP_35_H
