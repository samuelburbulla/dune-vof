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

  template< class DomainVector >
  void evaluate ( const DomainVector &x, const double t, double &u )
  {
    DomainVector center { 0.5, 0.5 };
    DomainVector offset { 0, 0.25 };
    DomainVector normal { -0.25, 0 };

    center.axpy( std::cos( ( 2 * M_PI / 10 ) * t ), offset );
    center.axpy( - std::sin( ( 2 * M_PI / 10 ) * t ), normal );

    double dist = ( x - center ).two_norm();

    u = (dist < 0.15) ? 1.0 : 0.0;
  }

  template< class DomainVector >
  void velocityField ( const DomainVector &x, const double t, DomainVector &r )
  {
    r[ 0 ] = x[ 1 ] - 0.5;
    r[ 1 ] = 0.5 - x[ 0 ];

    r *= 2 * M_PI / 10;
  }
};
#endif //#ifndef PROBLEM_HH
