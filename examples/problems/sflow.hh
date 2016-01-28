#ifndef PROBLEM_SFLOW_HH
#define PROBLEM_SFLOW_HH

//- dune-common includes
#include <dune/common/fvector.hh>


template < class ctype, int dim >
struct SFlow
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
    DomainType y{ 0.5, 0.5 };
    double radius = 0.25;

    y -= x;
    u = ( y.two_norm2() <= radius * radius ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &v ) const
  {
    DomainType y{ 4 * x[ 0 ] - 2,
                  4 * x[ 1 ] - 2 };

    v = {  0.25 * ( y[ 0 ] + y[ 1 ] * y[ 1 ] * y[ 1 ] ),
          -0.25 * ( y[ 1 ] + y[ 0 ] * y[ 0 ] * y[ 0 ] ) };
  }

  static double maxVelocity() { return sqrt(12.5); };

};

#endif // PROBLEM_SFLOW_HH
