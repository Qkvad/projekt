#ifndef DUNE_PARSOLVE_PROBLEMB_HH
#define DUNE_PARSOLVE_PROBLEMB_HH

// function for defining the scalar diffusion coefficient
template<typename GV, typename RF>
class k_A
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> >,
      k_A<GV,RF> >
{
public:
  typedef RF RFType;
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,k_A<GV,RF> > BaseT;

  k_A (const GV& gv_)
        : gv(gv_)
  {
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
        y = 1.0;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

private:
  const GV& gv;
};

// function for defining the diffusion tensor
template<typename GV, typename RF>
class K_A
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_A<GV,RF> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,
      GV::dimension*GV::dimension,Dune::FieldMatrix<RF,GV::dimension,GV::dimension> >,
      K_A<GV,RF> > BaseT;

  K_A (const GV& gv_) : gv(gv_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j)
          y[i][i] = 1.0;
        else
          y[i][j] = 0.0;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }
};

// function for defining the Helmholtz term
template<typename GV, typename RF>
class A0_A
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  A0_A<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,A0_A<GV,RF> > BaseT;

  A0_A (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                                                          typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

// function for defining the source term
template<typename GV, typename RF>
class F_A
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  F_A<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,F_A<GV,RF> > BaseT;

  F_A (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                                                          typename Traits::RangeType& y) const
  {
      double pow2x = x[0]*x[0];
      double pow2y = x[1]*x[1];
      double pow3x = pow2x*x[0];
      double pow3y = pow2y*x[1];
      double pow4x = pow3x*x[0];
      double pow4y = pow3y*x[1];
      double egz = x[0]*(x[0]-1)*x[1]*(x[1]-1)*exp(-pow2x-pow2y);
      double Laplace = -2*(2*pow4x - 2*pow3x -5*pow2x + 3*x[0] + 1)*x[1]*(x[1]-1)*exp(-pow2x-pow2y)
              - 2*(2*pow4y - 2*pow3y -5*pow2y + 3*x[1] + 1)*x[0]*(x[0]-1)*exp(-pow2x-pow2y);
      y = -Laplace + c(e,x)*egz;
  }
};



// constraints parameter class for selecting boundary condition type
class BCTypeParam_B
  : public Dune::PDELab::FluxConstraintsParameters,
        public Dune::PDELab::DirichletConstraintsParameters
        /*@\label{bcp:base}@*/
{
public:

  template<typename I>
  bool isNeumann(
                                   const I & intersection,   /*@\label{bcp:name}@*/
                                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                                   ) const
  {
    //Dune::FieldVector<typename I::ctype, I::dimension>
    //  xg = intersection.geometry().global( coord );
        return false;  // Dirichlet b.c. on ALL boundaries!
  }

  template<typename I>
  bool isDirichlet(
                                   const I & intersection,   /*@\label{bcp:name}@*/
                                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                                   ) const
  {
        return !isNeumann( intersection, coord );
  }

};


// boundary grid function selecting boundary conditions
template<typename GV>
class B_A
  : public Dune::PDELab::BoundaryGridFunctionBase<Dune::PDELab::
                                                  BoundaryGridFunctionTraits<
                                                    GV,Dune::PDELab::DiffusionBoundaryCondition::Type,1,
                                                    Dune::FieldVector<
                                                      Dune::PDELab::DiffusionBoundaryCondition::Type,1> >,
                                                  B_A<GV> >
{
  const GV& gv;

public:
  typedef Dune::PDELab::DiffusionBoundaryCondition BC;
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV,
    Dune::PDELab::DiffusionBoundaryCondition::Type,1,
    Dune::FieldVector<Dune::PDELab::DiffusionBoundaryCondition::Type,1> > Traits;
  typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,B_A<GV> > BaseT;

  B_A (const GV& gv_) : gv(gv_) {}

  template<typename I>
  inline void evaluate (const Dune::PDELab::IntersectionGeometry<I>& ig,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = BC::Dirichlet;
  }

  //! get a reference to the GridView
  inline const GV& getGridView ()
  {
    return gv;
  }
};

// function for Dirichlet boundary conditions and initialization
template<typename GV, typename RF>
class G_A
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  G_A<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,G_A<GV,RF> > BaseT;

  G_A (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                                                          typename Traits::RangeType& y) const
  {
        y = 0;
        for (int i=0; i<GV::dimension; i++)
          if (x[i]<1E-12 || x[i]>1-1E-12)
                {
                  double pow2x = x[0]*x[0];
                  double pow2y = x[1]*x[1];
                  y = x[0]*(x[0]-1)*x[1]*(x[1]-1)*exp(-pow2x-pow2y);

                }
        //y = exp(-x[0]-(x[1]*x[1]));
        return;
  }
};

// function for defining the flux boundary condition
template<typename GV, typename RF>
class J_A
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  J_A<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,J_A<GV,RF> > BaseT;

  J_A (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                                                          typename Traits::RangeType& y) const
  {
        y = 0;
        return;
  }
};

// flux as velocity field for the mixed method
template<typename GV, typename RF>
class V_A
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension>,
                                                                                                          V_A<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,V_A<GV,RF> > BaseT;

  V_A (const GV& gv) : BaseT(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                                                          typename Traits::RangeType& y) const
  {
    y = 0.0;
  }
};

#endif
