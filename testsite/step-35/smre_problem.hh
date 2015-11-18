#ifndef SMRE_PROBLEM_HH
#define SMRE_PROBLEM_HH

#include <deal.II/base/point.h>

#include <step-35/ADMMParams.h>

namespace step35 {





// @sect3{Class: ImageDoFHandler}
//
// Images or stacks of them (aka movies) are bounded domains in 2D or 3D. As there numerical representation is
// 1-dimensional, we need a structure which encapsulates the mapping from the multidimensional coordinates to the
// linear memory address. In analogy to finite element methods this is called a degrees-of-freedom handler.
// template<int dim>
class ImageDoFHandler {

public:
    static const int dim = 3;
    typedef dealii::Point<dim> Point;

    // The default CTor initializes everything to zero and requires to call reinit() at some later time.
    ImageDoFHandler() {}

    // The CTor takes the true dimensions of the image stack
    // and additionally a flag, whether these should be extended for accomodating some padding
    // in order to meet particular rquirements concerning memory alignment.
    // If the padding multiplier is set to a positive value, then this value is used for
    // computing the leading dimension for the fastest running index, that is x.
    // On input the z component of the bounding box has to be at least 1 as a 2D image represents one plane in 3D.
    ImageDoFHandler(const Point& bounding_box, int padding_mult=-1)
        :
          bbox(bounding_box),
          padding_multiplier(padding_mult)
    {
        this->reinit(bounding_box, padding_multiplier);
    }

    // The copy constructor for copying the DoF information around easily.
    ImageDoFHandler(const ImageDoFHandler & other)
        :
          bbox(other.bbox),
          padding_multiplier(other.padding_multiplier)
    {
        this->reinit(other.bbox, other.padding_multiplier);
    }

    // The reinit function does the same as the constructor.
    // On input the z component of the bounding box has to be at least 1 as a 2D image represents one plane in 3D.
    void reinit(const Point& bounding_box, int padding_mult=-1)
    {
        this->bbox = bounding_box;
        this->padding_multiplier = padding_mult;

        _leading_dim_x = int(bbox[0]);
        _leading_dim_y = int(bbox[1]);

        if (padding_multiplier > 0)
        {
            _leading_dim_x = padding_multiplier * (int(bbox[0] + padding_multiplier - 1)/padding_multiplier);
            _leading_dim_y = padding_multiplier * (int(bbox[1] + padding_multiplier - 1)/padding_multiplier);
        }

        // Last but not least, check whether all values are positive.
        AssertThrow(bbox[0] > 0, dealii::ExcMessage("x-dim is zero or negative."));
        AssertThrow(bbox[1] > 0, dealii::ExcMessage("y-dim is zero or negative."));
        AssertThrow(_leading_dim_y > 0, dealii::ExcMessage("y-dim is zero or negative."));
        AssertThrow(bbox[2] > 0, dealii::ExcMessage("z-dim is zero or negative."));

        _n_dofs_per_plane = _leading_dim_x * _leading_dim_y;
        _n_dofs = _n_dofs_per_plane * n_dofs_z();
        _n_active_dofs = n_dofs_x() * n_dofs_y() * n_dofs_z();
    }


    // We still need a few helper functions for convencience.
    int n_dofs() const { return _n_dofs; }

    int n_active_dofs() const { return _n_active_dofs; }

    // The true width of an image in the 2D plane.
    int n_dofs_x() const { return bbox[0]; }

    // The width after padding.
    int pwidth() const { return n_dofs_x_padded(); }

    // The true height of an image in the 2D plane and its padded counterpart.
    int n_dofs_y() const { return bbox[1]; }
    int pheight() const { return n_dofs_y_padded(); }

    // The height of a stack of images. Also designated as depth.
    int n_dofs_z() const { return bbox[2]; }
    int pdepth() const { return n_dofs_z(); }

    // Repeats of pwidth() and pheight() for situations when talking about the number of DoFs is more intuitive.
    int n_dofs_x_padded() const { return _leading_dim_x; }

    int n_dofs_y_padded() const { return _leading_dim_y; }



    // For the lexicographic ordering we assume @p y to the fastest running index, followed by @p x and then @p z.
     // Within the 2D plane ordering is thus column-major.
    int global_index(const int x, const int y, const int z=0) const
    {
        // Exceptions are expensive, hence we test the correctness of the input only in DEBUG mode.
        Assert(x >= 0 && x < _leading_dim_x, dealii::ExcIndexRange(x, 0, _leading_dim_x));
        Assert(y >= 0 && y < _leading_dim_y, dealii::ExcIndexRange(y, 0, _leading_dim_y));
        Assert(z >= 0 && z < bbox[2], dealii::ExcIndexRange(z, 0, bbox[2]));

        return z * _n_dofs_per_plane
                +
                x * _leading_dim_y
                +
                y;
    }

protected:
    int _leading_dim_x, _leading_dim_y, _n_dofs, _n_dofs_per_plane, _n_active_dofs;
    Point bbox;
    int padding_multiplier;
};




//@sect3{Class: SMREProblem}
//
// The overall problem of reconstructing an image using the SMRE approach is put into a class
// which administrates memory management, I/O and running the solver.
// From a structural point of view it is very similar to the examples of stationary PDEs in deal.II.
// To some extent this is no spurprise as mathematically we try to solve an inverse problem
// which basically follows the same pattern of solving a PDE. The only difference is that applying the
// operator representing the inverse problem usually means to compute a convolution or evaluating an
// integral equation rather than just computing derivatives.
// Much of the things in this class can be generalized to fit other use cases as well adn could be factored
// out into a base class InverseProblem. However, for the time being we do not do this and reserve this for a later step.
template <typename NumberType>
class SMREProblem {

    SMREProblem();

    void run(const ADMMParams & prm);

protected:
    void setup_system();

    void assemble();

    void solve();

    void output();

    NumberType dt, reg_strength;

};

}


// @sect3{Function: setup_system}
//
// We allocate the memory for the various solution vectors (i.e. images) and set the initial condition.
template <typename NumberType>
void step35::SMREProblem<NumberType>::setup_system()
{

}


template <typename NumberType>
void step35::SMREProblem<NumberType>::assemble() {}


template <typename NumberType>
void step35::SMREProblem<NumberType>::solve() {}


template <typename NumberType>
void step35::SMREProblem<NumberType>::output() {}


// @sect3{Function: run}
//
// This function drives the whole solurtion process.
template <typename NumberType>
void step35::SMREProblem<NumberType>::run(const ADMMParams &prm)
{


    this->setup_system();

    this->solve();

    this->output();


}






#endif // SMRE_PROBLEM_HH

