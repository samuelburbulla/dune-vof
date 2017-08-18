#ifndef HILLS_HH
#define HILLS_HH

#include <cmath>

#include <dune/common/fvector.hh>

template < class ctype, int dim >
struct Hills
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, RangeType &u ) const
  {
    return evaluate ( x, 0.0, u );
  }

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    double z = 0.5;
    for ( int i = 0; i < dim-1; ++i )
      z -= std::abs( x[ i ] - 0.5 );

    z += 0.25;
    u = ( x[dim-1] < z ) ? RangeType( 1.0 ) : RangeType( 0.0 );
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const {}

};

#endif //#ifndef HILLS_HH
