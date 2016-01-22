#ifndef PROBLEM_SLOTTEDCYLINDER_HH
#define PROBLEM_SLOTTEDCYLINDER_HH

//- dune-common includes
#include <dune/common/fvector.hh>


template < class ctype, int dim >
struct SlottedCylinder
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, RangeType &u ) const
  {
    evaluate( x, 0.0, u );
  }

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    DomainType center{ 0.5, 0.5};

    DomainType xPrime = center;
    DomainType tmp = x - center;
    xPrime.axpy( std::cos( ( 2.0 * M_PI / 10.0 ) * t ), tmp );
    tmp = { -tmp[ 1 ], tmp[ 0 ] };
    xPrime.axpy( std::sin( ( 2.0 * M_PI / 10.0 ) * t ), tmp );

    u = RangeType( 0.0 );

    double slotWidth = 0.075;
    double dist = ( xPrime - center ).two_norm();
    if ( dist < 0.4 )
      u = RangeType( 1.0 );

    dist = ( xPrime - ( center + DomainType{ 0.0, slotWidth } ) ).two_norm();
    if ( dist < slotWidth )
      u = RangeType( 0.0 );

    if ( xPrime[ 1 ] > center[ 1 ] + slotWidth && std::abs( xPrime[ 0 ] - center[ 0 ] ) < slotWidth )
      u = RangeType( 0.0 );
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot = x;
    rot -= DomainType{ 0.5, 0.5 };

    std::swap( rot[ 0 ], rot[ 1 ] );
    rot[ 1 ] = -rot[ 1 ];
    rot *= 2 * M_PI / 10;
  }

  static double maxVelocity() { return 2 * M_PI / 10 * sqrt( 2.0 ); };

};

#endif // PROBLEM_SLOTTEDCYLINDER_HH
