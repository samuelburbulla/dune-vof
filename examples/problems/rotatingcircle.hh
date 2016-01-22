#ifndef PROBLEM_ROTATINGCIRCLE_HH
#define PROBLEM_ROTATINGCIRCLE_HH

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>


template < class ctype, int dim >
struct RotatingCircle
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

#endif // #ifndef PROBLEM_ROTATINGCIRCLE_HH
