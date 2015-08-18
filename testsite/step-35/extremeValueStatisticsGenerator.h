//@sect3{File: extremeValueStatisticsGenerator.h}
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

#ifndef EXTREMVALSTATGEN_H_
#define EXTREMVALSTATGEN_H_

//std
#include<list>
#include<algorithm>
#include<sstream>

//Boost
#include <boost/mpl/modulus.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

//OpenMP
#ifdef USE_OMP
#include <omp.h>
#endif

//SciPAL
#include <base/ParallelArch.h>

//Our stuff
#include "cuda_driver_step-35.h"
#include "haar.hpp"
#include "patch.h"
#include "preprocessor_directives.h" // for DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION
const double HALF=0.5;
const double ONE=1.0;
const double TWO=2.0;
const double EIGHT=8.0;
#else
const float HALF=0.5;
const float ONE=1.0;
const float TWO=2.0;
const float EIGHT=8.0;
#endif

//@sect4{Struct: constructStruct}
//Template specialized struct for the extremeValueStatisticsGenerator class based on which Mpatch we're dealing with
template<typename T,typename Q, ParallelArch c> class extremeValueStatisticsGenerator;
template<typename Mpatch,typename T> struct constructStruct;
template<typename T>
//Specialisation for cset_small
struct constructStruct<cset_small,T> {
    //@sect5{Function: construct}
    //@brief This function generates a linear list of clusters with max. edge length n per patch, maxsize max_size per cluster
    //@param n maximum edge length per patch
    //@param maximum cluster size in px
    //@param nx width of the image
    //@param ny height of the image
    cluster<cset_small,T>* construct(int n, int max_size, int nx, int ny, int nz) {
        nz = nz; //suppress compiler warning
        cluster<cset_small,T> *root,*akt;
        akt=NULL;
        root=NULL;
        compareStruct<cset_small,T> cs;
        compareStruct<cluster<cset_small,T>,T> csc;
        //Loop over allowed sizes
        for (int size = 1 ; size<= n ; size++) {
            //Loop over allowed offsets
#ifdef USE_OMP
            #pragma omp parallel for
#endif
            for (int offset=0; offset<size; offset++) {
                //Thread private variables, the root cluster is appended to the
                //global last element when the thread finishes
                cluster<cset_small,T> *myakt,*myroot,*cprev;
                cprev=NULL;
                myakt=new cluster<cset_small,T>(max_size);
                myroot=myakt;
                cset_small *tmp,*prev;
                prev=NULL;
                tmp=NULL;
                std::list<int> indexlist;
                for (int i=offset; i<=nx-size; i+=size) {
                    if ( myakt->size > 0 ) {
                        //Check if there is a collision. There should be none ny construction
                        if ( cprev != NULL && !csc.compare(cprev,myakt) )
                            std::cout << "Warning, clusters after new since newline overlap\n";
                        cprev=myakt;
                        myakt->finalize();
                        myakt->next=new cluster<cset_small,T>(max_size);
                        myakt=myakt->next;
                    }
                    for (int j=offset; j<=ny-size; j+=size) {
                        indexlist.clear();
                        for (int i0=i; i0<i+size; i0++) {
                            for (int j0=j; j0<j+size; j0++) {
                                indexlist.push_back(i0*ny+j0);
                            }
                        }
                        prev=tmp;
                        tmp=new cset_small(indexlist,i,j,size);
                        //Check if there is a collision. There should be none ny construction
                        if ( prev != NULL  && cs.compare(prev,tmp)) {
                            std::cout << " Warning: subsequent patches overlap\n";
                        }
                        if ( myakt->add_patch(tmp) ) {
                            //When we get here it means we could NOT add the patch to the cluster
                            myakt->finalize();
                            if ( cprev != NULL && !csc.compare(cprev,myakt) )
                                std::cout << "Warning, clusters after new since full overlap\n";
                            cprev=myakt;
                            myakt->next=new cluster<cset_small,T>(max_size);
                            myakt=myakt->next;
                            if ( myakt->add_patch(tmp) ) {
                                //When we get here something went really wrong
                                std::cout << "Warning: patch dropped during cluster creation, too big\n";
                            }
                        }
                        tmp=tmp->next;
                    }
                }
                //Calculate the cluster info based on its patches. This makes cluster size estimation faster later on
                myakt->finalize();
                //Append the root cluster of this thread to the global linked list. Avoid collisions by doing this critical
#ifdef USE_OMP
                #pragma omp critical
#endif
                {
                    if ( akt != NULL ) {
                        akt->next=myroot;
                        akt=myakt;
                    } else {
                        root=myroot;
                        akt=myakt;
                    }
                }
            }
        }
        return root;
    }

};

//Specialisation for cset_dyadic
template<typename T>
struct constructStruct<cset_dyadic,T> {
    //@sect5{Function: construct}
    //@brief This function generates a linear list of clusters with max. edge length n per patch, maxsize max_size per cluster
    //@param n dummy
    //@param max_size dummy
    //@param nx image width
    //@param ny image height
    cluster<cset_dyadic,T>* construct(int n, int max_size, int nx, int ny,int nz) {
        n = n; //supress compiler warning
        nz = nz; //supress compiler warning

        std::list<int> indexlist;
        cset_dyadic *ctmp,*croot;
        cluster<cset_dyadic,T> *root;
        //This is just a dummy object used to keep a uniform structure interface
        root=new cluster<cset_dyadic,T>(max_size);
        bool first=true;
        //assume that the first layer is enough to make a good statistic
        const int z=0;
        //Loop over all squares with edge length 2^n up to n=6
        for (int shift =1; shift<=32; shift*=2) {
            for (int i=0; i<nx; i+=shift) {
                for (int j=0; j<ny; j+=shift) {
                    for (int i0=i; i0<i+shift; i0++) {
                        for (int j0=j; j0<j+shift; j0++) {
                            indexlist.push_back(i0*ny+j0+z*nx*ny);
                        }
                    }
                    if ( first ) {
                        croot=new cset_dyadic(indexlist,i,j,z,shift);
                        root->add_patch(croot);
                        ctmp=croot;
                        first=false;
                    } else {
                        ctmp=new cset_dyadic(indexlist,i,j,z,shift);
                    }
                    ctmp=ctmp->next;
                    indexlist.clear();
                }
            }
        }
        return root;
    }
};

//@sect4{Class: extremeValueStatisticsGenerator}
//@brief Generates the patches and puts them in clusters if needed,
//       calculates the weights $c_s$
template<typename T,typename Q, ParallelArch c>
//T patch type \n
//Q number type \n
//c gpu_cuda or cpu \n
class extremeValueStatisticsGenerator {
  public:
    Q *imageXY;
    Q nsize;

    //For T=cset_dyadic \n
    //First element in linear list of patches
    T *croot;
    //Current element in linear list of patches
    T *cpoint;

    //For T=cluster<cset_small, Q> \n
    //First element in linear list of cluster
    cluster<T,Q> *cluster_root;
    //Array of image dimensions
    int *n2;
    //Width of psf
    int sigma;
    //Width of image
    int nx;
    //Padded width of image
    int nx2,ny2,nz2;
    //Maximum number of iterations
    int max_it;
    //Height of image
    int ny;
    int nz;
    //Reduced height after R2C FFT
    int nyd;
    //Used in statistics generation
    int step;
    int split;

    //@sect5{Constructor}
    //Generates the multiresolution sets. In case of cset_dyadic sets are created for the quantile
    //estimation in order to calculate the c_s. In case of cset_small the sets are joined
    // in clusters.
    //@param nnx image width
    //@param nny image height
    //@param imnew image
    //@param new_sigma width of psf
    //@param new_step step to be used
    //@param max_size maximum patch edge length
    extremeValueStatisticsGenerator(int nnx, int nny,int nnz, std::vector<Q> & imnew, int new_sigma, int new_step, int max_size) {
        sigma=new_sigma;
        split=1;
        nx=nnx;
        ny=nny;
        nz=nnz;
        //Pad to a multiple of 32, actually useful for the approximative method
        while ( nx%32 != 0 )
            nx++;
        while ( ny%32 != 0 )
            ny++;
        while ( nz%32 != 0 )
            nz++;
        imageXY=&imnew[0];

        //Set up the convex sets
        constructStruct<T,Q> cs;
        cluster_root=cs.construct(new_step,max_size,nx,ny,nz);
        croot=cluster_root->patches.front();
    }

    //@sect5{Function: tsg}
    //@brief Used to evaluate the fourth root transform extreme value statistics
    //@param e field to evaluate
    //@param sigma gaussian noise standard deviation
    //@param indexlist list of positions to evaluate
    //@param size size of e
    Q tsg(Q *e,const Q sigma,int* indexlist,const Q size) {
        Q t=0;
        for (int i=0; i<size; i++) {
            t=t+boost::math::pow<2>(e[indexlist[i]]);
        }
        t=t/boost::math::pow<2>(sigma);
        return (sqrt(sqrt(t))-sqrt(sqrt((size)-HALF)))*sqrt(EIGHT*sqrt(size));
    }

    //@sect5{Function: get_quantile_gauss}
    //@brief This function generates the statistics necessary for an equalized generation of the weights $c_s$
    //@return alpha quantile of the statistics
    //@param si standard deviation of gaussian noise
    //@param alpha quantlie to calculate
    //@param mr multiresolution depth
    //@param nx width
    //@param ny height
    //@param qalpha_ret pointer to store quantile in
    Q* get_quantile_gauss(const Q si, const Q alpha, const int mr, const int nx, const int ny, /*const*/ int nz, Q* qalpha_ret) {
        nz = nz; // Supress compiler warning

        Q *cs_h=new Q[mr];
        Q mrreal;
        if ( alpha < 0 || alpha > ONE || si < 0) {
            std::cerr << "Invalid input to distribution simulation" << std::endl;
            return NULL;
        }
        //First check if the needed value is already written down in the bibfile in the working directory
        std::ifstream bibfile("/home/neal.hermer/Desktop/cuda_prakt_2015/SciPAL/testsite/build-step-35-Desktop-Release/bibfile_gauss.txt");
        if ( bibfile ) {
            int tmr,tnx,tny;
            Q talpha,tres,tsi;
            while ( !bibfile.eof() ) {
                bibfile >> tsi >> tmr >> tnx >> tny >> talpha >> tres;
                if ( !bibfile.eof() ) {
                    if ( si == tsi && tmr == mr && tnx == nx && tny == ny && talpha == alpha ) {
                        Q sigmas,mus;
                        qalpha_ret[0]=tres;
                        for (int i=0; i<mr; i++) {
                            mrreal=pow(2,i);
                            //calculate the weights
                            sigmas=ONE/(EIGHT*sqrt((Q)((mrreal+1)*(mrreal+1))));
                            mus=sqrt(sqrt((Q)((mrreal+1)*(mrreal+1))-HALF));
                            cs_h[i]=ONE/boost::math::pow<4>(tres*sigmas+mus);
                        }
                        return cs_h;
                    }
                }
            }
        }
        bibfile.close();

        //Since the necessary value is not in the bibfile start the calculations
        const int max_samples=1000;
        boost::mt19937 rng;
        boost::normal_distribution<> nd(0.0, si);
        boost::variate_generator<boost::mt19937&,
              boost::normal_distribution<> > var_nor(rng, nd);
        int w=nx-2*sigma;
        int h=ny-2*sigma;
        int samples=std::min(max_samples,w*h);
        Q *e=new Q[nx*ny];
        Q *k=new Q[w*h];
        Q tmp;
        Q qalpha=0;
        Q max=0;
        int steps=100;//a histogram with this detail is produced

        std::ofstream out("histogram_gauss.txt");
        //By default widht*height simulations of the statistics, in 3D case take only the first layer
        int scount=0;
        for (int i=0; i<w; i++) {
            for (int j=0; j<h; j++) {
                for (int ii=0; ii<nx; ii++) {
                    for (int jj=0; jj<ny; jj++) {
                        //simulate iid normal data
                        e[ii*ny+jj]=var_nor();
                    }
                }
                k[scount]=0;
                cpoint=croot;
#ifdef USE_OMP
                #pragma omp parallel
#endif
                {
#ifdef USE_OMP
                    int ID=omp_get_thread_num();
                    int n_ids = omp_get_num_threads();
#else
      int ID = 0;
      int n_ids = 1;
#endif

                    T *mycp;
                    mycp=croot;
                    // Why does this loop only run until the ID of the current thread?
                    for (int id=0; id<ID; id++) {
                        mycp=mycp->next;
                        if ( mycp == NULL)
                            break;
                    }
                    while ( mycp != NULL ) {
                        tmp=tsg(e,sigma,mycp->indexes,mycp->size);
                        //Synch with other threads
#ifdef USE_OMP
                        #pragma omp critical
#endif
                        {
                            if ( tmp > k[i*h+j] )
                                k[scount]=tmp;
                        }
                        mycp=mycp->next;
                        //Calculate offset in working queue for this thread
                        for (int id=1; id < n_ids; id++) {
                            if ( mycp != NULL)
                                mycp=mycp->next;
                            else
                                break;
                        }
                    }
                }
                if ( k[scount] > max )
                    max=k[scount];
                std::cout << "Calculating samples " << scount << "/" << samples << std::endl;
                scount++;
                if ( scount >= samples )
                    break;
            }
            if ( scount >= samples )
                break;
        }

        //Statistics ready, produce an output and evaluate the quantile
        Q *stat=new Q[steps];
        Q dk=max/((Q)steps);
        int counted=0;
        int limit=(int)round(alpha*((Q)(samples)));
        for (int i=0; i<steps; i++) {
            std::cout << "Calculating Histogram " << i << "/" << steps << std::endl;
            stat[i]=0;
            for (int ii=0; ii<samples; ii++) {
                if ( k[ii] >= dk*((Q)i) && k[ii] < dk*((Q)(i+1)) ) {
                    stat[i]++;
                    counted++;
                    if ( counted >= limit && qalpha == 0 )
                        qalpha=k[ii];
                }

            }
            out << dk*(HALF+(Q)i) << " " << stat[i] << std::endl;
        }

        //Since the result is not in the bib file write it there
        std::ofstream bib_out("bibfile_gauss.txt",std::ofstream::app);
        bib_out << si << " " <<  mr << " " << nx << " " << ny << " " << alpha << " " << qalpha << std::endl;
        Q sigmas,mus;
        qalpha_ret[0]=qalpha;
        for (int i=0; i<mr; i++) {
            mrreal=pow(2,i);
            sigmas=ONE/(EIGHT*sqrt((Q)((mrreal+1)*(mrreal+1))));
            mus=sqrt(sqrt((Q)((mrreal+1)*(mrreal+1))-HALF));
            cs_h[i]=ONE/boost::math::pow<4>(qalpha*sigmas+mus);
        }

        return cs_h;
    }
};
#endif /* EXTREMVALSTATGEN_H_ */
