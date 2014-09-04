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


#ifndef RNDWRAPPER_H
#define RNDWRAPPER_H
#include <curand.h>

//Boost includes for RNG for cpu
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/nondet_random.hpp>
#include <base/VTraits.h>
namespace SciPAL{
// @sect3{struct: RANDTraits}
//This struct holds the execution commands for different precisions.
template<typename T> struct RANDTraits;

template<>
struct RANDTraits<float>
{
    static void exec(curandGenerator_t &gen, float *dst, int size)
    {
        curandGenerateUniform(gen, dst, size);
    }
};

template<>
struct RANDTraits<double>
{
    static void exec(curandGenerator_t &gen, double *dst, int size)
    {
        curandGenerateUniformDouble(gen, dst, size);
    }
};

// @sect3{class: CUDARNDbase}
// This base class holds some fundamental informations needed for the FFT
// e.g. dimension, datatypes.
template <typename T, ParallelArch arch>
class CUDARNDbase {

public:
    void operator()(typename SciPAL::VTraits<T, cpu>::Vector &d_rndnumbers);
    curandGenerator_t gen;
    boost::mt19937 gen_cpu;
private:

};

// @sect3{class: CUDARND}
// Wrapper class for random number generator
template <typename T, ParallelArch arch>
struct CUDARND
 : CUDARNDbase<T,arch> {};

//Specializations for cuda / cpu
template <typename T>
struct CUDARND<T, cpu>
 : CUDARNDbase<T, cpu>
{
    CUDARND(int seed)
    {
        this->gen_cpu.seed(seed);
    }
    void operator()(typename SciPAL::VTraits<T, cpu>::Vector &d_rndnumbers)
    {
        boost::uniform_real<> range(0.0,1.0);
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rnd(this->gen_cpu, range);
        for(int i=0; i<d_rndnumbers.size(); i++)
        {
            //works with this Vectors only on cpu!
            d_rndnumbers.array().val()[i] = rnd();
        }
    }
};


//Specialization for CUDA.
template <typename T>
struct CUDARND<T, gpu_cuda>
        : CUDARNDbase<T, gpu_cuda>
{
    CUDARND(int seed)
    {
        curandCreateGenerator(&this->gen,CURAND_RNG_PSEUDO_DEFAULT);
        curandSetPseudoRandomGeneratorSeed(this->gen,seed);
    }
    ~CUDARND()
    {
        curandDestroyGenerator(this->gen);
    }

    void operator()(typename VTraits<T,gpu_cuda>::Vector &d_rndnumbers)
    {

        if(d_rndnumbers.size()!=0)
        {
            RANDTraits<T>::exec(this->gen, d_rndnumbers.array().val(), d_rndnumbers.size() );
        }
        else
        {
            std::cout<<"Destination vector for random numbers not initialized."<<std::endl;
        }
    }
};//end struct


}
#endif // RNDWRAPPER_H
