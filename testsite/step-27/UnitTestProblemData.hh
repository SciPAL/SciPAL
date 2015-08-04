#ifndef UNITTESTPROBLEMDATA_HH
#define UNITTESTPROBLEMDATA_HH

#include <deal.II/base/function.h>

#include <QDir>

//TODO: Runtime parameter
#define DK_PROTEIN 4.



namespace step27 {

// Because of the extensive list of test cases just for the FEM-BEM coupling
// we put the classes into a separate namespace in a separate file.
// In the long term this allows to separate the particular physical problem from the
// development aspects and their needs.
namespace UnitTests {



template <int dim>
class ProteinDipole : public dealii::Function<dim>
{

protected:
    struct Data {

        Data(const std::string &prm_path, const QDir &log_dir);

        void declare(dealii::ParameterHandler &prm);

        void get(dealii::ParameterHandler &prm);


        const double dipole_length;
        const double dipole_strength;

        // Centers of the charges.
        const dealii::Point<dim>  p0, p1;
    };

    Data data;

public:
    ProteinDipole (const std::string & prm_path, const QDir &log_dir,
                   const unsigned int n_components = 1)
        :
          dealii::Function<dim>(n_components),
          data(prm_path, log_dir)
    {}

    virtual double value (const dealii::Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;

    virtual void gradient_list(const std::vector<dealii::Point<dim> > &points,
                               std::vector<dealii::Tensor<1, dim> > &gradients,
                               const unsigned int component = 0) const;

    virtual dealii::Tensor<1,dim> gradient (const dealii::Point<dim> & p,
                                    const unsigned int component = 0) const;
};



// @sect3{The <code>HarmonicTest</code> class template}


template <int dim>
class HarmonicSum  : public dealii::Function<dim>
{

protected:
    static double strength() { return .1; }


    // Dummy structure. Needed for compatibility reasons.
    // Can be filled with soem fancy parameters in a future version.
    struct Data {};

public:
    HarmonicSum (const std::string & /*prm_path*/, const QDir &/*log_dir*/, const unsigned int n_components = 1)
        :
          dealii::Function<dim>(n_components) {}

    virtual double value (const dealii::Point<dim>   &p,
                          const unsigned int  component = 0) const
    {
        double sum = component;

        sum = 2*p(0);

        for (int d = 1; d < dim; ++d)
            sum += p(d);

        return strength() * sum;
    }

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const
    {
        const unsigned int n_points = points.size();

        Assert (values.size() == n_points,
                dealii::ExcDimensionMismatch (values.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            values[i] = HarmonicSum<dim>::value (points[i]);
    }


    virtual dealii::Tensor<1,dim> gradient (const dealii::Point<dim> &/* p*/,
                                    const unsigned int component = 0) const
    {
        dealii::Tensor<1,dim> grad;

        grad[0] = component; // dummy asignment

        switch(dim)
        {
        case 2: {
            grad[0] = 2.;
            grad[1] = 1.;
            return grad;
        }
        case 3:{
            grad[0] = 2.;
            grad[1] = 1.;
            grad[2] = 1.;

            return strength() * grad;
        }
        default:
            Assert(false, dealii::ExcInternalError());
            return dealii::Point<dim>();
        }
    }


    virtual void gradient_list (const std::vector<dealii::Point<dim> > &points,
                                std::vector<dealii::Tensor<1,dim> > & gradients,
                                const unsigned int component = 0)  const
    {
        const unsigned int n_points = points.size();

        Assert (gradients.size() == n_points,
                dealii::ExcDimensionMismatch (gradients.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            gradients[i] = HarmonicSum<dim>::gradient (points[i]);
    }
};




template <int dim>
class HarmonicProd : public dealii::Function<dim>
{
protected:
     static double strength() { return .01; }


     // Dummy structure. Needed for compatibility reasons.
     // Can be filled with soem fancy parameters in a future version.
     struct Data {};

public:

    HarmonicProd (const std::string & /*prm_path*/, const QDir &/*log_dir*/, const unsigned int n_components = 1)
        :
          dealii::Function<dim>(n_components) {}

    virtual double value (const dealii::Point<dim>   &p,
                          const unsigned int  component = 0) const
    {
        double prod = p(0);
        for (int d = 1; d < dim; ++d)
        {
            prod *=  p(d);
        }
        return strength() * prod;
    }

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const
    {
        const unsigned int n_points = points.size();

        Assert (values.size() == n_points,
                dealii::ExcDimensionMismatch (values.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            values[i] = HarmonicProd<dim>::value (points[i]);
    }


    virtual dealii::Tensor<1,dim> gradient (const dealii::Point<dim> & p,
                                    const unsigned int component = 0) const
    {
        dealii::Tensor<1,dim> grad;

        grad[0] = component;

        switch(dim)
        {
        case 2: {
            grad[0] = p(1);
            grad[1] = p(0);
            return  strength() * grad;
        }
        case 3:{
            grad[0] =       p(1) * p(2);
            grad[1] = p(0) *        p(2);
            grad[2] = p(0) * p(1);

            return strength() * grad;
        }
        default:
            Assert(false, dealii::ExcInternalError());
            return dealii::Point<dim>();
        }
    }


    virtual void gradient_list (const std::vector<dealii::Point<dim> > &points,
                                std::vector<dealii::Tensor<1,dim> > & gradients,
                                const unsigned int component = 0)  const
    {
        const unsigned int n_points = points.size();

        Assert (gradients.size() == n_points,
                dealii::ExcDimensionMismatch (gradients.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            gradients[i] = HarmonicProd<dim>::gradient (points[i]);
    }
};



template <int dim, typename F, typename G /*must be derived from dealii::Function<dim> */ >
class HarmonicExpr : protected F, protected G //, HarmonicProd<dim>, ProteinDipole<dim>
{
    typedef HarmonicExpr<dim, F, G> Expr;

public:
    HarmonicExpr (const std::string & prm_path, const QDir &log_dir, const unsigned int n_components = 1)
        :
          F(prm_path, log_dir, n_components),
          G(prm_path, log_dir, n_components)
        /* ProteinDipole<dim>(n_components) */{}

    virtual double value (const dealii::Point<dim>   &p,
                          const unsigned int  component = 0) const
    {
        return this->F::value(p, component)
                +
                this->G::value(p, component);
        /*+
                  this->ProteinDipole<dim>::value(p, component)*/;
    }

    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const
    {
        const unsigned int n_points = points.size();

        Assert (values.size() == n_points,
                dealii::ExcDimensionMismatch (values.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            values[i] = Expr::value (points[i]);
    }


    virtual dealii::Tensor<1,dim> gradient (const dealii::Point<dim> & p,
                                    const unsigned int component = 0) const
    {
        dealii::Tensor<1,dim> grad;

        switch(dim)
        {
        case 2:
        case 3:{

            grad =  this->F::gradient(p, component);

            grad += this->G::gradient(p, component);

            //            grad += this->ProteinDipole<dim>::gradient(p, component);

            return grad;

        }
        default:
            Assert(false, dealii::ExcInternalError());
            return dealii::Point<dim>();
        }
    }


    virtual void gradient_list (const std::vector<dealii::Point<dim> > &points,
                                std::vector<dealii::Tensor<1,dim> > & gradients,
                                const unsigned int component = 0)  const
    {
        const unsigned int n_points = points.size();

        Assert (gradients.size() == n_points,
                dealii::ExcDimensionMismatch (gradients.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            gradients[i] = Expr::gradient (points[i]);
    }
};



template<int dim>
struct Sum {

    typedef dealii::ZeroFunction<dim> NewtonPotential;

    typedef HarmonicSum<dim> Type;
};


template<int dim>
struct Prod {

    typedef dealii::ZeroFunction<dim> NewtonPotential;

    typedef HarmonicProd<dim> Type;
};



template<int dim>
struct Dipole {

    typedef ProteinDipole<dim> NewtonPotential;

    typedef ProteinDipole<dim> Type;
};



template<int dim>
struct SumProd {

    typedef dealii::ZeroFunction<dim> NewtonPotential;

    typedef
    HarmonicExpr<dim,
    HarmonicSum<dim>,
    HarmonicProd<dim> > Type;

};

template<int dim>
struct SumProdDipole {

    typedef ProteinDipole<dim> NewtonPotential;

    typedef
    HarmonicExpr<dim,
        HarmonicSum<dim>,
        HarmonicExpr<dim,
            HarmonicProd<dim>,
            ProteinDipole<dim> >
    > Type;

};



template<int dim, typename Expr >
class HarmonicTest
        :
        public  dealii::Function<dim> {

    typename Expr::Type expr;

public:

   // typename Expr::NewtonPotential newton_potential;


    HarmonicTest (const std::string & prm_path,
                  const QDir &log_dir,
                  const unsigned int n_components = 1)
        :
          dealii::Function<dim>(n_components),
          expr(prm_path, log_dir, n_components) // ,
       //   newton_potential(n_components)
    {}


    virtual double value (const dealii::Point<dim>   &p,
                          const unsigned int  component = 0) const
    {
        return expr.value(p, component);
    }




    virtual void value_list (const std::vector<dealii::Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const
    {
        const unsigned int n_points = points.size();

        Assert (values.size() == n_points,
                dealii::ExcDimensionMismatch (values.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        // dummy statement
        values[0] = component;

        for (unsigned int i=0; i<n_points; ++i)
            values[i] = expr.value (points[i]);
    }



    virtual dealii::Tensor<1,dim> gradient (const dealii::Point<dim> & p,
                                    const unsigned int component = 0) const
    {
        return expr.gradient(p, component);
    }



    virtual void gradient_list (const std::vector<dealii::Point<dim> > &points,
                                std::vector<dealii::Tensor<1,dim> > & gradients,
                                const unsigned int component = 0)  const
    {
        const unsigned int n_points = points.size();

        Assert (gradients.size() == n_points,
                dealii::ExcDimensionMismatch (gradients.size(), n_points));

        Assert (component == 0,
                dealii::ExcIndexRange (component, 0, 1));

        for (unsigned int i=0; i<n_points; ++i)
            gradients[i] = expr.gradient (points[i]);
    }

};









template<int dim>
void ProteinDipole<dim>::Data::declare( dealii::ParameterHandler & prm)
{
    prm.enter_subsection("Protein dipole data");

    prm.declare_entry("Dipole strength", "1.",
                      dealii::Patterns::Double(),
                      "The number of elementary charges.");

    prm.declare_entry("Dipole length", "1.",
                      dealii::Patterns::Double(),
                      "The length "
                      "is measured in the same units as "
                      "the dimensions of the FEM mesh. "
                      "Thus the numbers given here depend on how the physical model "
                      "has been made dimensionless.");

    prm.declare_entry("Position of positive charge", "Not yet implemented",
                      dealii::Patterns::Anything(),
                      "In a future version the positions of the charges "
                      "can be set directly");

    prm.declare_entry("Position of negative charge", "Not yet implemented",
                      dealii::Patterns::Anything(),
                      "In a future version the positions of the charges "
                      "can be set directly");

    prm.leave_subsection();
}



template<int dim>
void ProteinDipole<dim>::Data::get(dealii::ParameterHandler &prm)
{
    prm.enter_subsection("Protein dipole data");

    const_cast<double&> (dipole_strength) = prm.get_double("Dipole strength");

    const_cast<double&> (dipole_length) = prm.get_double("Dipole length");

    prm.leave_subsection();
}


template<int dim>
ProteinDipole<dim>::Data::Data(const std::string & prm_path, const QDir &log_dir)
   :
      dipole_length(1.), // dpl),
      dipole_strength(0) // dps / dk_p)
{

    // TOFDO: read params from file
    dealii::ParameterHandler prm_handler;

    std::string prm_file (prm_path + "/protein_dipole_data.prm");

    this->declare(prm_handler);

    prm_handler.read_input(prm_file);

    this->get(prm_handler);

    std::ofstream log_dipole_data_prm( (log_dir.absolutePath() + QDir::separator() + QFileInfo(prm_file.c_str()).fileName() + ".log").toStdString().c_str() );
    prm_handler.print_parameters (log_dipole_data_prm,
                                  dealii::ParameterHandler::Text);

    // Charges are placed along the x axis symmetric w.r.t. the origin.
    // Because the centers are declared as constant but have to be set
    // in the body of the constructor we have to temporarily remove constness.
    const_cast<dealii::Point<dim>&>(p0)(0) = (+ this->dipole_length * .5 /* this is for converting length to radius */
             /* make sure, that it is inside the sphere. */ );
    const_cast<dealii::Point<dim>&>(p1)(0) = (- this->dipole_length * .5);

    std::cout << "Dipole charge positions : " << p0 << ", " << p1 << "\n"
              << "dipole strength : " << this->dipole_strength << "\n"
              << "dipole length : " << this->dipole_length << std::endl;
}


template <int dim>
double ProteinDipole<dim>::value (const dealii::Point<dim> &p,
                                  const unsigned int) const
{
    double value = 0.;

    double r_0 = p.distance(data.p0);
    double r_1 = p.distance(data.p1);

    // Only for charges of finite distance the dipole moment is nonzero.
    if (std::fabs(this->data.dipole_length) > 1e-14)
        value = this->data.dipole_strength *(1./r_0 - 1./r_1 ); // / (4 * numbers::PI);

    return value;
}



template <int dim>
void ProteinDipole<dim>::value_list (const std::vector<dealii::Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int              component) const
{
    const unsigned int n_points = points.size();

    Assert (values.size() == n_points,
            dealii::ExcDimensionMismatch (values.size(), n_points));

    Assert (component == 0,
            dealii::ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
        values[i] = ProteinDipole<dim>::value (points[i]);
}


template <int dim>
dealii::Tensor<1,dim> ProteinDipole<dim>::gradient (const dealii::Point<dim> & p,
                                            const unsigned int /*component*/) const
{
    dealii::Tensor<1,dim> grad;

    switch(dim)
    {
    case 3:{
        grad = 0.;

        // Only for charges of finite distance the dipole moment is nonzero.
        if (std::fabs(this->data.dipole_length) > 1e-14)
        {
            double r_0 = p.distance(data.p0);
            double r_1 = p.distance(data.p1);

            grad -= (p - data.p0)/std::pow(r_0, 3.);
            grad += (p - data.p1)/std::pow(r_1, 3.);

            grad *=  this->data.dipole_strength; // / (4 * numbers::PI);
        }

        return grad;
    }
    default:
        Assert(false, dealii::ExcInternalError());
        return dealii::Point<dim>();
    }
}



template <int dim>
void ProteinDipole<dim>::gradient_list (const std::vector<dealii::Point<dim> > &points,
                                        std::vector<dealii::Tensor<1,dim> > & gradients,
                                        const unsigned int component)  const
{
    const unsigned int n_points = points.size();

    Assert (gradients.size() == n_points,
            dealii::ExcDimensionMismatch (gradients.size(), n_points));

    Assert (component == 0,
            dealii::ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
        gradients[i] = ProteinDipole<dim>::gradient (points[i]);
}



} // namespace UnitTests END


} // namespace step27 END





#endif // UNITTESTPROBLEMDATA_HH
