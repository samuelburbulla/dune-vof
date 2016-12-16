#ifndef PROBLEM_SFLOW_HH
#define PROBLEM_SFLOW_HH

//- dune-common includes
#include <dune/common/fvector.hh>


template < class ctype, int dim >
struct SFlow {};


template < class ctype >
struct SFlow< ctype, 2 >
{
  using DomainType = Dune::FieldVector< ctype, 2 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = 2 };
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

template < class ctype >
struct SFlow< ctype, 3 >
{
  using DomainType = Dune::FieldVector< ctype, 3 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = 3 };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, RangeType &u ) const
  {
    evaluate( x, 0.0, u );
  }

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    DomainType y{ 0.5, 0.5, 0.5 };
    double radius = 0.25;

    y -= x;
    u = ( y.two_norm2() <= radius * radius ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &v ) const
  {
    DomainType y{ 4 * x[ 0 ] - 2,
                  4 * x[ 1 ] - 2,
                  0.0 };

    v = {  0.25 * ( y[ 0 ] + y[ 1 ] * y[ 1 ] * y[ 1 ] ),
          -0.25 * ( y[ 1 ] + y[ 0 ] * y[ 0 ] * y[ 0 ] ),
          0.0 };
  }

  static double maxVelocity() { return sqrt(12.5); };

};

#endif // PROBLEM_SFLOW_HH
