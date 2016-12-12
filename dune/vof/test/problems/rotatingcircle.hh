#ifndef ROTATINGCIRCLE_HH
#define ROTATINGCIRCLE_HH

#include <cmath>

template < class ctype, int dim >
struct RotatingCircle {};


template < class ctype >
struct RotatingCircle < ctype, 2 >
{
  using DomainType = Dune::FieldVector< ctype, 2 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  DomainType center( const double t ) const
  {
    DomainType center = rotationCenter();
    DomainType offset { 0, 0 };
    DomainType normal { 0, 0 };

    center.axpy( std::cos( ( 2 * M_PI / 10 ) * t ), offset );
    center.axpy( std::sin( ( 2 * M_PI / 10 ) * t ), normal );

    return center;
  }

  ctype radius( const double t ) const
  {
    return 0.4;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &r ) const
  {
    DomainType center = rotationCenter();

    r[ 0 ] = x[ 1 ] - center[ 1 ];
    r[ 1 ] = center[ 0 ] - x[ 0 ];

    r *= 2 * M_PI / 10;
  }

  ctype curvature( const DomainType &x ) const
  {
    return 1.0 / radius( 0.0 );
  }

private:
  DomainType rotationCenter () const
  {
    return DomainType { 0.5, 0.5 };
  }

};

template < class ctype >
struct RotatingCircle < ctype, 3 >
{
  using DomainType = Dune::FieldVector< ctype, 3 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = 3 };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, RangeType &u ) const
  {
    return evaluate ( x, 0.0, u );
  }

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    DomainType center{ 0.5, 0.5, 0.5 };
    center.axpy( std::cos( (2 * M_PI / 10)*t ), DomainType{ 0.25, 0.0, 0.0 } );
    center.axpy( std::sin( (2 * M_PI / 10)*t ), DomainType{ 0.0,  0.25, 0.0 } );
    center.axpy( - std::sin( (2 * M_PI / 10)*t ), DomainType{ 0.0, 0.0, 0.25 } );

    double dist = ( x - center ).two_norm();

    u = ( dist < 0.15 ) ? RangeType( 1.0 ) : RangeType( 0.0 );
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    DomainType c = x;
    c -= DomainType{ 0.5, 0.5, 0.5 };

    rot[ 0 ] = - c[ 1 ] + c[ 2 ];
    rot[ 1 ] = c[ 0 ];
    rot[ 2 ] = - c[ 0 ];
    rot *= 2 * M_PI / 10 / std::sqrt(2);
  }

};

#endif //#ifndef ROTATINGCIRCLE_HH
