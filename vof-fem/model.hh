#ifndef MODEL_HH
#define MODEL_HH

template< class Problem >
struct Model
{
  using ProblemType = Problem;

  using DomainType =
    typename ProblemType::DomainType;
  using RangeType =
    typename ProblemType::RangeType;
  using ctype =
    typename ProblemType::ctype;

  Model ()
  : problem_( ProblemType() )
  {}

  ctype numericalFlux ( const DomainType &normal,
                        double t,
                        const DomainType &x,
                        const RangeType &uInside, const RangeType &uOutside,
                        RangeType &flux ) const
  {
    // Local Lax-Friedrichs Flux
    RangeType fInside, fOutside;
    ctype ws = std::max( analyticalFlux( normal, t, x, uInside,  fInside  ),
                         analyticalFlux( normal, t, x, uOutside, fOutside ) );
    flux = fInside + fOutside;
    flux.axpy( ws, uInside - uOutside );
    flux *= 0.5;

    /*
    std::cout 
              << std::setw( 6 ) << (double) normal
              << std::setw( 12 ) << (double) x
              << std::setw( 12 ) << (double) uInside
              << std::setw( 12 ) << (double) fInside
              << std::setw( 12 ) << (double) uOutside
              << std::setw( 12 ) << (double) fOutside 
              << std::setw( 6 ) << (double) ws 
              << std::setw( 12 ) << (double) flux << std::endl;
    */
    return ws;
  }

  ctype boundaryFlux ( const DomainType &normal,
                       double t,
                       const DomainType &x,
                       const RangeType &u,
                       RangeType &flux ) const
  {
    RangeType uBnd;
    problem().evaluate( x, t, uBnd );
    return numericalFlux( normal, t, x, u, uBnd, flux );
  }

  ctype analyticalFlux ( const DomainType &normal,
                         double t,
                         const DomainType &x,
                         const RangeType &u,
                         RangeType &flux ) const
  {
    return problem().analyticalFlux( x, t, normal, u, flux );
  }

  const ProblemType& problem () const { return problem_; }

private:
  const ProblemType problem_;
};

#endif // MODEL_HH
