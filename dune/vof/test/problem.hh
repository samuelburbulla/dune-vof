#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>

template < class ctype, int dim >
struct RotatingCircle {};


template < class ctype >
struct RotatingCircle < ctype, 2 >
{
  using DomainType = Dune::FieldVector< ctype, 2 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  DomainType center( const double t )
  {
    DomainType center { 0.5, 0.5 };
    DomainType offset { 0, 0.25 };
    DomainType normal { 0.25, 0 };

    center.axpy( std::cos( ( 2 * M_PI / 10 ) * t ), offset );
    center.axpy( std::sin( ( 2 * M_PI / 10 ) * t ), normal );

    return center;
  }

  ctype radius( const double t )
  {
    return 0.15;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &r )
  {
    r[ 0 ] = x[ 1 ] - 0.5;
    r[ 1 ] = 0.5 - x[ 0 ];

    r *= 2 * M_PI / 10;
  }
};
#endif //#ifndef PROBLEM_HH
