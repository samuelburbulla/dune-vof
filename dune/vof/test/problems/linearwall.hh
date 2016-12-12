#ifndef LINEARWALL_HH
#define LINEARWALL_HH

#include <cmath>

template < class ctype, int dim >
struct LinearWall
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    double pos = 0.4;
    double width = 0.3;

    u = ( std::abs( x[0] - t * speed_ - pos ) <= width ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot[ 0 ] = speed_;
    rot[ 1 ] = 0.0;
  }

private:
  double speed_ = 0.02;

};

#endif //#ifndef LINEARWALL_HH
