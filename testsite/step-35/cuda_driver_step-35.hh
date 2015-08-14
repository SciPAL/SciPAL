//@sect3{File: cuda_driver_step-35.hh}
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

//Helper classes for CUDADriver class
#ifndef CUDA_DRIVER_STEP_35_HH
#define CUDA_DRIVER_STEP_35_HH

//CUDA
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>

//std stuff
#include <list>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <algorithm>    // std::min

//CUDA cuFFT
#include <cufft.h>

//Our stuff
#include "cuda_driver_step-35.h"
#include "cuda_kernel_wrapper_step-35.cu.h"
#include "patch.h"
#include "cuda_helper.h"
#include "preprocessor_directives.h"
#include "ADMMParams.h"

#include <deal.II/lac/vector.h>

namespace step35 {
//@sect4{Class: info}
//
//This class contains the collected information about the image we're processing, like the width,
//height, maximum frame size, pointers to the image and $c_s$ on host and device side. The boolean array
//occupied has the size of the image, one element is true if this pixel is currently proccessed on the
//device and false otherwise.
template<typename Mdouble, typename Mcomplex, ParallelArch c>
class ImageInfo {
  public:
    //@param width width of the image
    //@param height height of the image
    //@param depth, in case of 3d, number of images
    const int width, height, depth;
    //@param ext_width width extended power of 2 to the next integer to next bigger multiple of 32
    //@param ext_height height extended to next bigger multiple of 32
    //@param ext_depth number of images extended to next bigger multiple of 32
    //@param nx2 ext_width extended to the next integer of two
    //@param ny2 ext_height extended to the next integer power of two
    //@param ext_num_pix equals ext_width*ext_height*ext_depth
    int ext_width, ext_height, nx2, ny2, ext_depth, ext_num_pix;
    //Residual Field on Device and Host
    //Mdouble *e_d; // FIXME,*e_h;
    SciPAL::Vector<Mdouble,cublas> e_d;
    //Noisy image on Device and Host
    // Mdouble *im_h; // FIXME
    // FIXME: storign references to vectors is bad design. Either use dealii::SmartPointers or rethink the design of this class. Preferably the latter.
    std::vector<Mdouble> &original_image;
    //Estimate and Device and Host
    //Mdouble *x_d; // FIXME,*x_h;
    SciPAL::Vector<Mdouble,cublas> x_d;
    //Smoothed Estimate on Device and Host
    //Mdouble *z_d; // FIXME,*z_h;
    SciPAL::Vector<Mdouble,cublas> z_d;
    //Vector of weights and point spread function on host
    // FIXME: replace by dealii::Vectors
    Mdouble *cs_h,*cs_d; // ,*psf_h;

    std::vector<Mdouble> x_h, e_h, z_h, im_h;


    dealii::Vector<Mdouble> psf_h;

    //Point spread function on device
    Mcomplex *fpsf_d;
    //Noisy image on device
    //Mdouble *im_d;
    SciPAL::Vector<Mdouble,cublas> im_d;
    //Smoothing parameter
    const Mdouble gamma;
    //psf width
    const int sigma;
    //Size for allocation of number and complex type arrays
    int n_bytes_per_frame, cframesize;
    //Regularisation type
    const regularisationType regType;
    //Dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    const int dim;
    //Current CUDA device id
    int device_id;

    //@sect5{Constructor: info}
    //@brief Class constructor
    //@param im pointer to the image on host
    //@param psf psf, will be fourier transformed
    //@param cs_vect array containing all $c_s$ values
    //@param nx width of the image
    //@param ny height of the image
    //@param cs_n maximum frame size
    //@param mgamma regularization parameter w.r.t. $x$
    //@param msigma estimated gauss variance
    //@param newreg regularisation Method
    //@param mdim dimension of the images (if TIFF-stack process each slice individually if dim = 2)
    ImageInfo(std::vector<Mdouble> & input_image, std::vector<Mdouble> &psf, Mdouble* cs_vec,
                      const int nx, const int ny, const int nz, const int cs_n,
                      Mdouble mgamma, int msigma, regularisationType newreg,
                      int mdim)
       :
         original_image(input_image), width(nx), height(ny), depth(nz),
         // FIXME:
         ext_width(32*((width+31)/32)),
         ext_height( 32*((height+31)/32) ), // height;
         ext_depth(depth), //32*((depth+31)/32) ),
         ext_num_pix(ext_width*ext_height*ext_depth),
         n_bytes_per_frame(ext_num_pix*sizeof(Mdouble)),
         cframesize(ext_depth*ext_width*(ext_height/2+1)*sizeof(Mcomplex)),
         psf_h(ext_num_pix),
         x_h(ext_num_pix),
         e_h(ext_num_pix),
         z_h(ext_num_pix),
         im_h(ext_num_pix),
         gamma(mgamma), sigma(msigma),regType(newreg),
         dim(mdim),
         e_d(ext_num_pix),
         x_d(ext_num_pix),
         z_d(ext_num_pix),
         im_d(ext_num_pix)
    {

        // FIXME: this is already in the gpuInfo object instantiated in the main function!
        //Set CUDA device id
        checkCudaErrors(cudaGetDevice(&this->device_id));

        //Calc padding sizes, has to be multiple of 32
        // FIXME: this is weird. CTors have the purpose to initialize as much as possible in the initializer line, see above.
        /*
        ext_width  = width;
        ext_height = height;
        ext_depth  = depth;
       while ( ext_height%32 != 0 )
            ext_height++;
        while ( ext_width%32 != 0 )
            ext_width++;
        if (ext_depth > 1) {
            while ( ext_depth%32 != 0 )
                ext_depth++;
        }*/


        //Number of total pixels
        // FIXME: bad practice
//        ext_num_pix = ext_width*ext_height*ext_depth;
  //      framesize = ext_num_pix*sizeof(Mdouble);
        // cframesize=ext_depth*ext_width*(ext_height/2+1)*sizeof(Mcomplex);
        // psf_h=new Mdouble[ext_num_pix];

        // FIXME: do padding when reading tiff -> add option to tiffreader
        //Padding
        for (int z=0; z<ext_depth; z++) {
            for (int i=0; i<ext_width; i++) {
                for (int j=0; j<ext_height; j++) {
                    if ( z < depth && i < width && j < height)
                        psf_h[z*ext_width*ext_height + i*ext_height + j] =
                                psf[z*width*height + i*height + j];
                    else
                        psf_h[z*ext_width*ext_height + i*ext_height + j] = 0;
                }
            }
        }

        //For Haar Regularisation, edge length has to be power of two
        nx2 = 1;
        ny2 = 1;
        while ( nx2 < ext_width || ny2 < ext_height  ) {
            nx2=nx2*2;//pow of two
            ny2=ny2*2;//pow of two
        }

        // FIXME: new is bad practice
        /*
        x_h=new Mdouble[ext_num_pix];
        //Now copy the device psf ft in here

        e_h  = new Mdouble[ext_num_pix];
        z_h  = new Mdouble[ext_num_pix];
        im_h = new Mdouble[ext_num_pix];
        */

        // FIXME: if you use std:vector or dealii::Vector the following is redundant!
        /*
        for (int z=0; z<ext_depth; z++) {
            for (int x=0; x<ext_width; x++) {
                for (int y=0; y<ext_height; y++) {
                    im_h[x*ext_height+y]=0;
                }
            }
        }
        */
        // FIXME: I do not want this written in components
        for (int z=0; z<depth; z++) {
            for (int x=0; x<width; x++) {
                for (int y=0; y<height; y++) {
                    if ( z < depth && x < width && y < height )
                        im_h[z*ext_width*ext_height + x*ext_height + y] =
                                input_image[z*width*height + x*height + y];
                }
            }
        }
        for (int i=0; i<ext_num_pix; i++) {
            x_h[i]=im_h[i];
            e_h[i]=0;
            z_h[i]=0;
        }
        // im_o=image;
        cs_h=cs_vec;
        // FIXME: is this a bug? loacl variable shadows class attribute.
        // Mdouble* cs_d;
        if ( c == gpu_cuda ) {
            //Allocate and copy needed arrays on gpu
            checkCudaErrors(cudaMalloc((void **)&cs_d, cs_n*sizeof(Mdouble)));
            //checkCudaErrors(cudaMalloc((void **)&x_d, n_bytes_per_frame));
            //checkCudaErrors(cudaMalloc((void **)&z_d, n_bytes_per_frame));
            //checkCudaErrors(cudaMalloc((void **)&e_d, n_bytes_per_frame));
            //checkCudaErrors(cudaMalloc((void **)&im_d, n_bytes_per_frame));
            Mdouble *fpsf_tmp_d;
            checkCudaErrors(cudaMalloc((void **)&fpsf_tmp_d, n_bytes_per_frame));
            checkCudaErrors(cudaMalloc((void **)&fpsf_d, cframesize));
            im_d = im_h;
            x_d = x_h;
            z_d = z_h;
            e_d = e_h;
            //checkCudaErrors(cudaMemcpyAsync(im_d.array().val(), &(im_h[0]), n_bytes_per_frame, cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpyAsync(fpsf_tmp_d, &(psf_h(0)), n_bytes_per_frame, cudaMemcpyHostToDevice));
            //checkCudaErrors(cudaMemcpyAsync(x_d.array().val(), &(x_h[0]), n_bytes_per_frame, cudaMemcpyHostToDevice));
            //checkCudaErrors(cudaMemcpyAsync(z_d.array().val(), &(z_h[0]), n_bytes_per_frame, cudaMemcpyHostToDevice));
            //checkCudaErrors(cudaMemcpy/*Async*/(e_d.array().val(), &(e_h[0]), n_bytes_per_frame, cudaMemcpyHostToDevice); //);
            checkCudaErrors(cudaMemcpyAsync(cs_d, cs_h, cs_n*sizeof(Mdouble), cudaMemcpyHostToDevice));
            //Copy $c_s$ to constant CUDA memory
            step35::Kernels<Mdouble> kernel;
            kernel.set_cs(cs_h, cs_n);

            checkCudaErrors(cudaDeviceSynchronize());


            //Fourier transformation of psf using cuFFT
            cufftHandle plan;
#ifdef DOUBLE_PRECISION
            cufftPlan3d(&plan, ext_depth, ext_width, ext_height, CUFFT_D2Z);
            cufftExecD2Z(plan, fpsf_tmp_d, fpsf_d);
#else
            cufftPlan3d(&plan, ext_depth, ext_width, ext_height,  CUFFT_R2C);
            cufftExecR2C(plan, fpsf_tmp_d, fpsf_d);
#endif

            //Cleanup
            cufftDestroy(plan);
            getLastCudaError("cufft error!\n");
            cudaFree(fpsf_tmp_d);
        }
    }

    //@sect5{Destructor: info}
    ~ImageInfo() {
        if ( c == gpu_cuda ) {
            cudaFree(cs_d);
            //cudaFree(e_d);
//            cudaFree(x_d);
//            cudaFree(z_d);
//            cudaFree(im_d);
            cudaFree(fpsf_d);
        }
    }

    //@sect5{Function get_result}
    //@brief Copy back image from device to host
    //@param im_host target where to write the image to
    void get_result(Mdouble* im_host) {
        if ( c == gpu_cuda ) {
            checkCudaErrors(cudaMemcpy(im_host, im_d, n_bytes_per_frame, cudaMemcpyDeviceToHost));
            for (int z=0; z<depth; z++) {
                for (int i=0; i<width; i++) {
                    for (int j=0; j<height; j++) {
                        original_image[z*width*height + i*height + j] =
                                im_host[z*ext_width*ext_height + i*ext_height + j];
                    }
                }
            }
        }
    }
};
// @sect4{Class: wqueue}
//
// This class implements a thread-safe FIFO queue with a maximum capacity.
template <typename Mpatch, typename Mdouble, ParallelArch c>
class wqueue {
    std::list<Mpatch*> m_queue;
    //Mutex to access m_queue
    std::mutex m_mutex;
    //Threads in add() will wait on this signal if queue if full
    std::condition_variable not_full;
    //Threads in get() will wait on this signal if queue if empty
    std::condition_variable not_empty;

    //List of items which are currently in m_queue and being processed by the worker threads
    std::list<Mpatch*> in_process;
    //Mutex to access in_process
    std::mutex p_mutex;
    //There might be a thread waiting to add items to the queue due to conflicts with the in_process queue
    std::condition_variable conflict;
  public:
    std::atomic<unsigned int> capacity;

    //@sect5{Constructor: wqueue}
    //@param mcapacity maximum capacity
    wqueue(unsigned int mcapacity) {
        capacity = mcapacity;
    }

    //@sect5{Function: add}
    //@brief Blocking append item to queue. This call is blocking, if the queue is full
    //       it will wait until one item is removed.
    //@param item appended to queue
    void add(Mpatch* item) {
        //Lock m_queue
        std::unique_lock<std::mutex> lock(m_mutex);

        //Blocking wait if queue is full
        not_full.wait(lock, [this]() {
            return m_queue.size() != capacity;
        });

        //Lock in_process
        std::lock_guard<std::mutex> guard(p_mutex);
        //Add to in_process list
        in_process.push_back(item);
        //Add to queue
        m_queue.push_back(item);

        //Signal waiting thread in get() that there's a new item
        not_empty.notify_one();
    }

    //@sect5{Function: blocking_add}
    //@brief Blocking wait until threre's no conflict with any other item in in_process
    //@param item item appended to queue
    void blocking_add(Mpatch* item) {
        //Lock in_process
        std::unique_lock<std::mutex> lock(p_mutex);
        //Block until there's no conflict
        conflict.wait(lock, [&]() {
            return this->__check_free(item);
        });

        lock.unlock();
        add(item);
    }

    //@sect5{Function: get}
    //@brief Blocking remove item from queue. This call is blocking, if the queue is
    //       empty it will wait until one item is added.
    //@return first item in queue.
    Mpatch* get() {
        std::unique_lock<std::mutex> lock(m_mutex);

        //Blocking wait if queue is empty
        not_empty.wait(lock, [this]() {
            return m_queue.size() != 0;
        });

        Mpatch* item = m_queue.front();
        m_queue.pop_front();

        //Signal waiting thread in add() that a place is available
        not_full.notify_one();
        return item;
    }

    //@sect5{Function: get(min, max)}
    //@brief Get and remove a list of items from queue. This call is blocking, if the queue is empty it will wait until one item is added.
    //@param min minimum number of elements it should preferably return. This is no strict condition, it may not be fullfilled.
    //           The call might be blocking for some time perdiod, after that it will return no more than max elements.
    //@param max maximum number of elements to be returned.
    //@return a list of size 0 to max items.
    std::list<Mpatch*> get(unsigned int min, unsigned int max) {
        std::unique_lock<std::mutex> lock(m_mutex);

        if (min) {
            //Blocking wait if queue has less entries than min for maximum 10 milliseconds
            not_empty.wait_for(lock, std::chrono::milliseconds(10), [&]() {
                return m_queue.size() > min-1;
            });
        }

        std::list<Mpatch*> list;
        typename std::list<Mpatch*>::iterator it = m_queue.begin();
        int size = std::min(max, (unsigned int) m_queue.size());

        //Check if queue was empty, if so return empty list
        if (!size) {
            return list;
        }

        //Move elements [0, size) from m_queue to list
        std::advance(it, size);
        list.splice(list.begin(), m_queue, m_queue.begin(), it);

        //Signal waiting thread in add() that a place is available
        not_full.notify_one();

        return list;
    }

    //@sect5{Function: size}
    //@return Current number of items in queue
    int size() {
        //Lock m_mutex
        std::lock_guard<std::mutex> guard(m_mutex);
        return m_queue.size();
    }

   //@sect5{Function: task_done}
   //@brief Removes item from in_process list
   //@param item item to remove from in_process
    void task_done(Mpatch* item) {
        //Lock in_process
        std::unique_lock<std::mutex> lock(p_mutex);
        in_process.remove(item);
        conflict.notify_one();
    }

    //Same for a list of items

    //@sect5{Function: task_done}
    //@brief Removes list of items from in_process list
    //@param list list of items to remove from in_process
    void task_done(std::list<Mpatch*> list) {
        //Lock in_process
        std::unique_lock<std::mutex> lock(p_mutex);
        for (typename std::list<Mpatch*>::iterator it=list.begin(), end=list.end(); it!=end; ++it) {
            in_process.remove(*it);
        }
        conflict.notify_one();
    }

  private:
    //@sect5{Function: __check_free}
    //@return return false if conflict \n
    //        return true if no conflict
    bool __check_free(Mpatch* item) {
        compareStruct<Mpatch, Mdouble> cs;
        int k=1;
        for (typename std::list<Mpatch*>::iterator it=in_process.begin(), end=in_process.end(); it!=end; ++it, k++) {
            if ( !cs.compare((*it), item) ) {
                //Conflict
                return false;
            }
        }

        //No conflict
        return true;
    }
};

}
#endif // CUDA_DRIVER_STEP_35_HH
