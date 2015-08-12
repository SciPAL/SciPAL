#ifndef STACK_H
#define STACK_H

#include <lac/expression_template.h>
#include <lac/expression_templates_host.h>
#include <complex>
#include <float.h>

//SciPAL implementation for expression templates

//#include <lac/SciPAL_kernels_wrapper.cu.h>

#include <lac/Array.h>
#include <lac/Shape.h>

#include <lac/BlasVectorOperations.h>

#include <base/PrecisionTraits.h>
#include <base/Zero_One_Traits.h>

#include <base/ForewardDeclarations.h>

namespace SciPAL {

// @sect3{Class: Stack}
//!
//! A class for aggregating SciPAL entities like Vector Matrix and Views
template<typename T>
class Stack
        :
        public  dealii::Subscriptor
{
public:
    typedef T Type;

    typedef typename Type::Number Number;

    typedef typename Type::value_type value_type;

    typedef typename Type::blas_wrapper_type las_wrapper_type;

    typedef const Type& ConstHandle;

private:
    size_t n_elements;
    //this variable keeps track, if an element of the stack is initialized
    std::vector<bool> initialized;
    std::vector<Type*> elements;

public:
    Stack():n_elements(0), initialized(0)
    {}

    Stack(size_t _n_elements)
        :n_elements(_n_elements), initialized(n_elements, false), elements(n_elements)
    {
        for(size_t ii = 0 ; ii < n_elements; ii++)
            elements[ii] = new Type();
    }

    Stack(size_t _n_elements, size_t init_elements)
        :n_elements(_n_elements), initialized(n_elements), elements(n_elements)
    {
        for(size_t ii = 0 ; ii < n_elements; ii++)
            elements[ii] = new Type(init_elements);
    }

    ~Stack() {}

    const size_t size()
    {
        return n_elements;
    }

    void resize(size_t _n_elements)
    {

        n_elements = _n_elements;

        elements.resize(n_elements);
        initialized.resize(n_elements);
        for(size_t ii = 0 ; ii < n_elements; ii++)
            elements[ii] = new Type();
    }

    //copy from other stack
    Stack& operator =(Stack &other)
    {
        if( n_elements != other.size())
            this->resize(other.size());

        for(size_t ii=0; ii<n_elements; ii++)
            (*elements[ii]) = other[ii];

        return *this;
    }

    //copy from std::vector
    template<typename T2>
    Stack& operator =(std::vector<T2> &other)
    {
        if( n_elements != other.size())
            this->resize(other.size());

        for(size_t ii=0; ii<n_elements; ii++)
            (*elements[ii]) = other[ii];

        return *this;
    }


    Type& operator [](size_t ind)
    {
        return (*elements[ind]);
    }

     const Type&  operator [](size_t ind) const
    {
        return (*elements[ind]);
    }

};
}

#endif // STACK_H

