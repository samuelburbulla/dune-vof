#ifndef LINEARWALL_HH
#define LINEARWALL_HH

#include <cmath>

#include <dune/common/fvector.hh>

template < class ctype, int dim >
struct LinearWall
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    u = ( x[0] <= pos + t * speed_  ) ? 1.0 : 0.0;
  }

  RangeType jump ( double t ) const
  {
    return pos + t * speed_;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot = 0;
    rot[ 0 ] = speed_;
  }

private:
  double speed_ = 0.02;
  double pos = 0.3;

};

#endif //#ifndef LINEARWALL_HH
