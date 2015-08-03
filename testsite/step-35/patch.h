//@sect3{File: patch.h}
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

#ifndef PATCH_H
#define PATCH_H

//std
#include<list>
#include<vector>

//CUDA
#include <cuda_runtime.h>

//Our stuff
#include "cuda_helper.h"

//@sect4{Class: cset_small}
//@brief Holds one patch of the image for exact calculation
//@param i0 starting line
//@param j0 starting row
//@param n edge length of patch
class cset_small {
  public:
    //Array of positions in the set
    int* indexes;
    //The csets are ordered in a linked list
    cset_small *next;
    unsigned int i0,j0,size,m,n;
    //@sect5{Constructor}
    //@param indexlist list of positions to form the patch on
    //@param new_i0 new initial line
    //@param new_j0 new initial row
    //@param new_n new edge length
    cset_small(std::list<int> indexlist, const unsigned int new_i0, const unsigned int new_j0, const unsigned int new_n) {
        //Init next as NULL pointer, used to check for end of list
        next=NULL;
        n=new_n;
        m=1;
        //In order to sum over an array the next power of 2 bigger than its number of pixels is needed
        while ( m < n )
            m=m*2;
        size=indexlist.size();
        //Indexes are handed in as list, and written to an array
        indexes=new int[size];
        i0=new_i0;
        j0=new_j0;
        int i=0;
        for (std::list<int>::iterator it=indexlist.begin(); it!=indexlist.end(); ++it) {
            indexes[i]=*it;
            i++;
        }
    }
    //@sect5{Deconstructor}
    ~cset_small() {
        delete[] indexes;
    }
};

//@sect4{Class: cset_dyadic}
//Contains major chunks of size 32x32 pixel for approximative dysktra \n
//This class is fed to the statistics generation tool
class cset_dyadic {
  public:
    //Array of image indizes
    int* indexes;
    //Pointer to next element in list
    cset_dyadic *next;
    //starting point in image, edge length, size, next power of 2 bigger than size
    unsigned int i0,j0,n,size,m,z;
    //@sect5{Constructor}
    //@param indexlist list of indexes in row-major format
    //@param new_i0 origin row
    //@param new_j0 new origin col
    //@param new_n new edge length
    cset_dyadic(std::list<int> indexlist, const int new_i0, const int new_j0, const int z0, const int new_n) {
        next=NULL;
        n=new_n;
        m=1;
        while ( m < n )
            m=m*2;
        size=indexlist.size();
        //Indexlist is given as list and written to array
        indexes=new int[size];
        i0=new_i0;
        j0=new_j0;
        z=z0;
        int i=0;
        for (std::list<int>::iterator it=indexlist.begin(); it!=indexlist.end(); ++it) {
            indexes[i]=*it;
            i++;
        }
    }
    //@sect5{Constructor}
    //This constructor provides consistency in the patches interface during construction
    cset_dyadic(int zero) {
        next=NULL;
        n=0;
        m=0;
        size=0;
        j0=0;
        i0=zero;
        z=0;
        indexes=NULL;
    }
    //@sect5{Destructor}
    ~cset_dyadic() {
        delete[] indexes;
    }

    //@sect5{Function: reset_q}
    //Set $q$ to zero, provides consistency in the class interface with cset_small
    int reset_q() {
        return 0;
    }
};

//@sect4{Struct: compareStruct}
//Specialized struct to compare if two objects of type Mpatch are
//conflicting and shouldn't be calculated at the same time. \n
//A struct is necessary to have a unique interface in the driver class for both cases
template<typename Mpatch,typename T> class cluster;
template<typename Mpatch,typename T> struct compareStruct;
template<typename T>
//Specialisation for cset_small
struct compareStruct<cset_small,T> {
    //@sect5{Function: compare}
    //@return true if no overlap (no conflict),
    //@param c1 first cset_small to be compared
    //@param c2 second cset_small to be compared
    //@param pt cluster where c1 resides and c2 is meant to be added
    //        false if overlap (conflict)
    bool compare(cset_small *c1, cset_small *c2, cluster<cset_small,T> *pt) {
        if ( ( c1->i0+c1->n-1 < c2->i0 || c1->j0+c1->n-1 < c2->j0 ||
                c1->i0 > c2->i0+c2->n-1 || c1->j0 > c2->j0+c2->n-1 )
                && pt->size+(c2->n*c2->n) <= pt->maxsize ) {
            return true;
        } else {
            return false;
        }
    }

    //@sect5{Function: compare}
    //@return true if no overlap (no conflict),
    //@param c1 first cset_small to be compared
    //@param c2 second cset_small to be compared
    //        false if overlap (conflict)
    bool compare(cset_small *c1, cset_small *c2) {
        if  ( c1->i0+c1->n-1 < c2->i0 || c1->j0+c1->n-1 < c2->j0 ||
                c1->i0 > c2->i0+c2->n-1 || c1->j0 > c2->j0+c2->n-1  )
            return true;
        else
            return false;
    }

    //@sect5{Function: getBoundarys}
    //@param ivec array to store initial row points in
    //@param jvec array to store initial col points in
    //@param nvec array to store edge length in
    //@param it pointer to current cset_small object
    //set start pixel i0,j0 and frame size n
    void getBoundarys(int *ivec, int *jvec, int *nvec, int i, cset_small *it) {
        ivec[i]=it->i0;
        jvec[i]=it->j0;
        nvec[i]=it->n;
    }

    //@sect5{Function: boundaryUpdate}
    //@param pt pointer to current cluster
    //@param c1 pointer to current cset_small
    //calculate a big rectangle which includes all cset_small
    void boundaryUpdate(cluster<cset_small,T> *pt,cset_small *c1) {
        pt->i0=std::min(pt->i0,c1->i0);
        pt->j0=std::min(pt->j0,c1->j0);
        pt->i1=std::max(pt->i1,c1->i0+c1->n-1);
        pt->j1=std::max(pt->j1,c1->j0+c1->n-1);
        pt->size+=(c1->n*c1->n);
    }
};

//Specialisation for cluster<cset_small>
template<typename T>
struct compareStruct<cluster<cset_small,T>,T> {
    //@sect5{Function: compare}
    //@param c1 first cluster to be compared
    //@param c2 second cluster to be compared
    //@return true if no overlap (no conflict),
    //        false if overlap (conflict)
    bool compare(cluster<cset_small,T> *c1,cluster<cset_small,T> *c2) {
        if  ( c1->i1 < c2->i0 || c1->j1 < c2->j0 || c1->i0 > c2->i1 || c1->j0 > c2->j1 )
            return true;
        else
            return false;
    }
};

//Specialisation for cset_dyadic
template<typename T>
struct compareStruct<cset_dyadic,T> {

    //@sect5{Function: compare}
    //@param c1 first cset_dyadic to be compared
    //@param c2 second cset_dyadic to be compared
    //@param pt cluster c1 resides in
    //@return true if no overlap (no conflict),
    //        false if overlap (conflict)
    bool compare(cset_dyadic *c1,cset_dyadic *c2,cluster<cset_dyadic,T> *pt) {
        if ( ( c1->i0+c1->n-1 < c2->i0 || c1->j0+c1->n-1 < c2->j0 || c1->i0 > c2->i0+c2->n-1 || c1->j0 > c2->j0+c2->n-1 ) && pt->size+(c2->n*c2->n) <= pt->maxsize )
            return true;
        else
            return false;
    }

    //@sect5{Function: boundaryUpdate}
    //@param pt cluster of cset_dyadic
    //@param c1 cset_dyadic inside the cluster
    //just for completeness of the specialization, does the same as for cset_small
    void boundaryUpdate(cluster<cset_dyadic,T> *pt,cset_dyadic *c1) {
        if ( c1->i0 < pt->i0 )
            pt->i0=c1->i0;
        if ( c1->j0 < pt->j0 )
            pt->j0=c1->j0;
        if ( c1->i0+c1->n > pt->i1 )
            pt->i1=c1->i0+c1->n;
        if ( c1->j0+c1->n > pt->j1 )
            pt->j1=c1->j0+c1->n;
        pt->size=pt->size+(c1->n*c1->n);
    }
};

//@sect4{Class: cluster}
//@brief A vector of Mpatches with additional functionality.
//We use this for Mpatch=cset_small to get a list of cset_small with
//a limited total maximum number of pixels. maxsize should can not be larger
//than 1024, the CUDA implementation will apply Dykstra's projection algorithm
//to all cset_small in one cluster object.
//As we use one thread per pixel, one cluster object will be handled by one threadblock.
//In the CUDADriver we will generate a list of cluster objects so that we can start a
//kernel with serveral threadblocks each independently dealing with a cluster object.
template<typename Mpatch,typename T>
class cluster {
  public:
    //The smallest rectangle in the plane including all points of the cluster is drawn
    //from (i0,j0) to (i1,j1) \n
    //size is the current size \n
    //maxsize is the maximum size
    unsigned int i0,j0,i1,j1,size,maxsize;
    //Vector of patches in the cluster
    std::vector<Mpatch*> patches;
    //Structure to enable comparison of templatised patches
    compareStruct<Mpatch ,T> cs;
    //The clusters are ordered as linked list, hence a pointer to next is needed
    cluster<Mpatch,T> *next;
    //Each cluster holds the variables $q_s$ for all of its patches
    T *qmat;
    //Pointer on device pointing to number of patches integer
    int *pnum_d;
    //Arrays holding position and size info for all patches
    int *ivec,*jvec,*nvec;
    //Arrays on device holding the patch infos,
    //copying the arrays at the beginning instead of giving it as argument
    //prevents small and inefficient memcopys to the device
    int *ivec_d,*jvec_d,*nvec_d;

    //@sect5{Constructor}
    //@brief Set variables to something, real work is done by add_patch
    //@param new_maxsize total maximum number of pixels allowed in the cluster
    cluster(int new_maxsize) {
        //init a negative area A=(i1-i0)*(j1-j0) < 0 in this setting \n
        //this can be used to identify improperly clusters created but not filled with anything
        i0=10000;
        j0=10000;
        i1=0;
        j1=0;
        size=0;
        maxsize=new_maxsize;
        next=NULL;
        qmat=NULL;
        pnum_d=0;
    }

    //@sect5{Destructor}
    //Clean up memory on host and device
    ~cluster() {
        cudaFreeHost(qmat);
        delete[] ivec;
        delete[] jvec;
        delete[] nvec;

        cudaFree(ivec_d);
        cudaFree(jvec_d);
        cudaFree(nvec_d);
        cudaFree(pnum_d);
    }

    //@sect5{Function: add_patch}
    //@brief Adds patch to cluster, returns true if not possible due to conflict
    //@return true = failure, false = success
    bool add_patch(Mpatch *new_patch) {
        //Add patch address to list
        for (typename std::vector<Mpatch*>::iterator it=patches.begin(); it!=patches.end(); ++it) {
            if ( !cs.compare((*it),new_patch,this) ) {
                //return true if adding fails
                return true;
            }
        }
        patches.push_back(new_patch);
        //Update the borders and the size
        cs.boundaryUpdate(this,new_patch);
        return false;
    }

    //@sect5{Function: finalize}
    //@brief Needed after cluster is "complete", update class variales and copy class object
    //       information to device so that we don't always have to do that
    void finalize() {
        //Allocate the needed memory and give each patch a little chunk
        //page-locked memory instead of qmat=new T[size]; as we need to copy it often
        checkCudaErrors(cudaHostAlloc((void **)&qmat, size*sizeof(T), cudaHostAllocDefault));

        ivec=new int[patches.size()];
        jvec=new int[patches.size()];
        nvec=new int[patches.size()];

        compareStruct<Mpatch,T> cs;
        int i=0;
        for (typename std::vector<Mpatch*>::iterator it=patches.begin(); it!=patches.end(); i++, ++it) {
            cs.getBoundarys(ivec,jvec,nvec,i,*it);
        }

        int pnum_h[1];
        pnum_h[0]=patches.size();

        //We store store a copy of the cluster information on the device
        //so that we don't need to copy it every time
        cudaMalloc((void **)&ivec_d, patches.size()*sizeof(int));
        cudaMalloc((void **)&jvec_d, patches.size()*sizeof(int));
        cudaMalloc((void **)&nvec_d, patches.size()*sizeof(int));
        cudaMalloc((void **)&pnum_d, sizeof(int));
        cudaMemcpyAsync(pnum_d, pnum_h, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpyAsync(ivec_d, ivec, patches.size()*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpyAsync(jvec_d, jvec, patches.size()*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpyAsync(nvec_d, nvec, patches.size()*sizeof(int), cudaMemcpyHostToDevice);
    }

    //@sect5{Function: reset}
    //@brief Set all $q$ values to zero
    void reset() {
        for (unsigned int i=0; i<size; i++)
            qmat[i]=0;
    }
};

#endif // PATCH_H
