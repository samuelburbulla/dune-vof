#ifndef SLOPE_HH
#define SLOPE_HH

#include <cmath>

#include <dune/common/fvector.hh>

template< class ctype, int dim > struct Slope;

template< class ctype >
struct Slope< ctype, 2 >
{
  using DomainType = Dune::FieldVector< ctype, 2 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  Slope ( double angle = 0.5 * M_PI )
    : normal_{ std::sin( angle ), -std::cos( angle ) }
  {}

  void evaluate ( const DomainType& x, RangeType& u ) const
  {
    const auto val = normal_ * ( x - DomainType( { 0.5, 0.5 } ) );
    u = val < 0.0 ? 1.0 : 0.0;
  }

  void evaluate ( const DomainType& x, double t, RangeType& u ) const
  {
    evaluate( x, u );
  }

  void velocityField ( const DomainType& x, double t, DomainType& v ) const
  {
    v = 0.0;
  }

  double curvature ( const DomainType& x ) const
  {
    return 0.0;
  }

  const DomainType& normal() const { return normal_; }

private:
  DomainType normal_;
};

#endif // #ifndef SLOPE_HH
