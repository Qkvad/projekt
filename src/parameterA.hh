#ifndef DUNE_PARAMETERA_HH
#define DUNE_PARAMETERA_HH

template<typename GV, typename RF>
class ParameterA
{
  const GV gv;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef RF RangeFieldType;
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  ParameterA( const GV gv_ ) : gv(gv_)
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
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::RangeFieldType norm = xglobal.two_norm2();
    // izracunati egzaktno koristeci Laplacea egzaktnog rjesenja

    // Primjer  1
//    double egz = exp(-xglobal[0]-xglobal[1]*xglobal[1]);
//    double Laplace = egz + (4*xglobal[1]*xglobal[1] - 2)*egz;
//    return -Laplace + c(e,x)*egz;

    // Primjer 2
    double pow2x = xglobal[0]*xglobal[0];
    double pow2y = xglobal[1]*xglobal[1];
    double pow3x = pow2x*xglobal[0];
    double pow3y = pow2y*xglobal[1];
    double pow4x = pow3x*xglobal[0];
    double pow4y = pow3y*xglobal[1];
    double egz = xglobal[0]*(xglobal[0]-1)*xglobal[1]*(xglobal[1]-1)*exp(-pow2x-pow2y);
    double Laplace = -2*(2*pow4x - 2*pow3x -5*pow2x + 3*xglobal[0] + 1)*xglobal[1]*(xglobal[1]-1)*exp(-pow2x-pow2y)
            - 2*(2*pow4y - 2*pow3y -5*pow2y + 3*xglobal[1] + 1)*xglobal[0]*(xglobal[0]-1)*exp(-pow2x-pow2y);
    return -Laplace + c(e,x)*egz;

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


    // egzaktna rjesenja
    // Primjer 1
    //return exp(-xglobal[0]-xglobal[1]*xglobal[1]);

    // Primjer 2
    double pow2x = xglobal[0]*xglobal[0];
    double pow2y = xglobal[1]*xglobal[1];
    return xglobal[0]*(xglobal[0]-1)*xglobal[1]*(xglobal[1]-1)*exp(-pow2x-pow2y);
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


#endif // DUNE_PARAMETERA_HH
