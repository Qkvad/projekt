#ifndef DUNE_PARAMETERC_HH
#define DUNE_PARAMETERC_HH

#include <cmath>

template<typename GV, typename RF>
class ParameterC
{
  const GV gv;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParameterC( const GV gv_ ) : gv(gv_)
  {
  }

  std::string name() const {return "A";};

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 1 : 0;
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if(xglobal[0]>0 && xglobal[0]<1 && xglobal[1]>0 && xglobal[1]<1) I*= 5;
    if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]<0 && xglobal[1]>-1) I*= 5;
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
//    typename Traits::DomainType xglobal = e.geometry().global(x);
//    typename Traits::RangeFieldType norm = xglobal.two_norm2();
//    // izracunati egzaktno koristeci Laplacea egzaktnog rjesenja

//    double r = std::sqrt(xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1]);
//    double theta =  atan(xglobal[1]/xglobal[0]);

//    double delta = 0.5354409456;
//    double a, b;
//    double K = 1.0;
//    if(xglobal[0]>0 && xglobal[0]<1 && xglobal[1]>0 && xglobal[1]<1) {
//        a = 0.4472135955;
//        b = 1.0;
//        K = 5.0;
//    }
//    else if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]<0 && xglobal[1]>-1) {
//        a = -0.7453559925;
//        b = 2.333333333;
//        K = 5.0;
//    }
//    else if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]>0 && xglobal[1]<1) {
//        a = -0.9441175905;
//        b = 0.55555555555;
//    }
//    else  {
//        a = -2.401702653;
//        b = -0.4814814814;
//    }
//    double temp = a*sin(delta*theta) + b*cos(delta*theta);
//    double egz = pow(r,delta) * temp;
//    double Laplace = delta * (delta - 1) * pow(r,delta-2) *temp - pow(delta, 2) * pow(r, delta) * temp;
//    return - K * Laplace + c(e,x) * egz;
      return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::RangeFieldType norm = xglobal.two_norm2();


    // egzaktno rjesenje
    double r = std::sqrt(xglobal[0]*xglobal[0] + xglobal[1]*xglobal[1]);
    double theta =  atan(xglobal[1]/xglobal[0]);

    double delta = 0.5354409456;
    double a, b;
    if(xglobal[0]>0 && xglobal[0]<1 && xglobal[1]>0 && xglobal[1]<1) {
        a = 0.4472135955;
        b = 1.0;
    }
    else if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]<0 && xglobal[1]>-1) {
        a = -0.7453559925;
        b = 2.333333333;
    }
    else if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]>0 && xglobal[1]<1) {
        a = -0.9441175905;
        b = 0.55555555555;
    }
    else  {
        a = -2.401702653;
        b = -0.4814814814;
    }
    double temp = a*sin(delta*theta) + b*cos(delta*theta);
    return pow(r,delta) * temp;
  }

  //! Neumann boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

};


#endif // DUNE_PARAMETERC_HH
