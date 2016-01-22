#ifndef PROBLEM_ROTATINGCIRCLE_HH
#define PROBLEM_ROTATINGCIRCLE_HH

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>

template< class FunctionSpace >
struct RotatingCircle
{
  static_assert( Dune::AlwaysFalse< FunctionSpace >::value, "No Specialization found" );
};

template< class DomainField, class RangeField >
struct RotatingCircle< Dune::Fem::FunctionSpace< DomainField, RangeField, 2, 1 > >
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
    DomainType center{ 0.5, 0.5 };
    center.axpy( std::cos( (2 * M_PI / 10)*t ), DomainType{ 0.0,  0.25 } );
    center.axpy( std::sin( (2 * M_PI / 10)*t ), DomainType{ 0.25, 0.0  } );

    double dist = ( x - center ).two_norm();

    u = ( dist < 0.15 ) ? RangeType( 1.0 ) : RangeType( 0.0 );
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

template< class DomainField, class RangeField >
struct RotatingCircle< Dune::Fem::FunctionSpace< DomainField, RangeField, 3, 1 > >
{
  using FunctionSpaceType = Dune::Fem::FunctionSpace< DomainField, RangeField, 3, 1 >;

  using DomainType = typename FunctionSpaceType::DomainType;
  using RangeType = typename FunctionSpaceType::RangeType;
  using ctype = RangeField;

  enum { dimDomain = 3 };
  enum { dimRange = 1 };

  static_assert( Dune::AlwaysFalse< DomainField >::value , "No Specialization found." );
};

#endif // PROBLEM_ROTATINGCIRCLE_HH
