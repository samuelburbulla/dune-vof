#ifndef PROBLEM_LINEARWALL_HH
#define PROBLEM_LINEARWALL_HH

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>

template< class FunctionSpace >
struct LinearWall
{
  static_assert( Dune::AlwaysFalse< FunctionSpace >::value, "No Specialization found" );
};

template< class DomainField, class RangeField >
struct LinearWall< Dune::Fem::FunctionSpace< DomainField, RangeField, 2, 1 > >
{
  using FunctionSpaceType = Dune::Fem::FunctionSpace< DomainField, RangeField, 2, 1 >;

  using DomainType = typename FunctionSpaceType::DomainType;
  using RangeType = typename FunctionSpaceType::RangeType;
  using ctype = RangeField;

  enum { dimDomain = 1 };
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

template< class DomainField, class RangeField >
struct LinearWall< Dune::Fem::FunctionSpace< DomainField, RangeField, 3, 1 > >
{
  using FunctionSpaceType = Dune::Fem::FunctionSpace< DomainField, RangeField, 3, 1 >;

  using DomainType = typename FunctionSpaceType::DomainType;
  using RangeType = typename FunctionSpaceType::RangeType;
  using ctype = RangeField;

  enum { dimDomain = 3 };
  enum { dimRange = 1 };

  static_assert( Dune::AlwaysFalse< DomainField >::value , "No Specialization found." );
};

#endif // PROBLEM_LINEARWALL_HH
