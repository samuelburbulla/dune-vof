#ifndef PROBLEM_LINEARWALL_HH
#define PROBLEM_LINEARWALL_HH

//- dune-common includes
#include <dune/common/fvector.hh>


template < class ctype, int dim >
struct LinearWall
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;
  using ctype = RangeField;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, RangeType &u ) const
  {
    evaluate( x, 0.0, u );
  }

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    double pos = 0.1;
    double width = 0.05;

    u = ( std::abs( x[0] - t * maxVelocity() - pos ) <= width ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot[ 0 ] = maxVelocity();
    rot[ 1 ] = 0.0;
  }

  static double maxVelocity() { return 0.08; };

};

#endif // PROBLEM_LINEARWALL_HH
