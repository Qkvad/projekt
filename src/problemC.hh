#ifndef DUNE_PARSOLVE_PROBLEMC_HH
#define DUNE_PARSOLVE_PROBLEMC_HH

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
    typename Traits::DomainType xglobal = e.geometry().global(x);
    y = 1.0;
    if(xglobal[0]>0 && xglobal[0]<1 && xglobal[1]>0 && xglobal[1]<1) y = 5;
    if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]<0 && xglobal[1]>-1) y = 5; 
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
    typename Traits::DomainType xglobal = e.geometry().global(x);
    for (int i=0; i<GV::dimension; i++)
      for (int j=0; j<GV::dimension; j++)
        if (i==j){
            y[i][i] = 1.0;
            if(xglobal[0]>0 && xglobal[0]<1 && xglobal[1]>0 && xglobal[1]<1) y[i][i] = 5;
            if(xglobal[0]<0 && xglobal[0]>-1 && xglobal[1]<0 && xglobal[1]>-1) y[i][i] = 5;
        }
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
//    double r = std::sqrt(x[0]*x[0] + x[1]*x[1]);
//    double theta =  atan(x[1]/x[0]);

//    double delta = 0.5354409456;
//    double a, b;
//    double K = 1.0;
//    if(x[0]>0 && x[0]<1 && x[1]>0 && x[1]<1) {
//        a = 0.4472135955;
//        b = 1.0;
//        K = 5.0;
//    }
//    else if(x[0]<0 && x[0]>-1 && x[1]<0 && x[1]>-1) {
//        a = -0.7453559925;
//        b = 2.333333333;
//        K = 5.0;
//    }
//    else if(x[0]<0 && x[0]>-1 && x[1]>0 && x[1]<1) {
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
//    y = - K * Laplace + c(e,x) * egz;
      y=0.0;
  }
};



// constraints parameter class for selecting boundary condition type
class BCTypeParam_C
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
                    double r = std::sqrt(x[0]*x[0] + x[1]*x[1]);
                    double theta =  atan(x[1]/x[0]);

                    double delta = 0.5354409456;
                    double a, b;
                    if(x[0]>0 && x[0]<1 && x[1]>0 && x[1]<1) {
                        a = 0.4472135955;
                        b = 1.0;
                    }
                    else if(x[0]<0 && x[0]>-1 && x[1]<0 && x[1]>-1) {
                        a = -0.7453559925;
                        b = 2.333333333;
                    }
                    else if(x[0]<0 && x[0]>-1 && x[1]>0 && x[1]<1) {
                        a = -0.9441175905;
                        b = 0.55555555555;
                    }
                    else  {
                        a = -2.401702653;
                        b = -0.4814814814;
                    }
                    double temp = a*sin(delta*theta) + b*cos(delta*theta);
                    y = pow(r,delta) * temp;

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
