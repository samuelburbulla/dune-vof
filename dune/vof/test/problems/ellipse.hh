#ifndef ELLIPSE_HH
#define ELLIPSE_HH

#include <cmath>
#include <array>

#include <dune/common/fvector.hh>


template< class ctype, int dim >
struct Ellipse
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  static constexpr int dimDomain = dim;

  Ellipse ( const std::array< DomainType, dim >& axis, const std::array< double, dim >& radii )
    : axis_( axis ), radii_( radii )
  {}

  void evaluate ( const DomainType& x, RangeType& u ) const
  {
    double val = 0;
    for( int i = 0; i < dimDomain; ++i )
    {
      const auto tmp = axis_[ i ] * ( x - DomainType( 0.5 ) );
      val += tmp * tmp / ( radii_[ i ] * radii_[ i ] );
    }

    u = val < 1.0 ? 1.0 : 0.0;
  }

  void evaluate ( const DomainType& x, double t, RangeType& u ) const
  {
    evaluate( x, u );
  }

  void velocityField ( const DomainType& x, double t, DomainType& v ) const
  {
    v = 0.0;
  }

  double curvature ( const DomainType& x ) const
  {
    double curv = 0.0;

    const auto tmp = axis_[ 0 ] * ( x - DomainType( 0.5 ) );
    curv += tmp * tmp * radii_[ 1 ] * radii_[ 1 ] * radii_[ 1 ] * radii_[ 1 ];

    const auto tmp2 = axis_[ 1 ] * ( x - DomainType( 0.5 ) );
    curv += tmp2 * tmp2 * radii_[ 0 ] * radii_[ 0 ] * radii_[ 0 ] * radii_[ 0 ];

    curv = curv * curv * curv;
    curv = std::sqrt( curv );

    for( int i = 0; i < dimDomain; ++i )
      curv /= radii_[ i ] * radii_[ i ] * radii_[ i ] * radii_[ i ];

    return 1.0 / curv;
  }

  std::array< DomainType, dim > axis_;
  std::array< double, dim > radii_;
};

#endif // #ifndef ELLIPSE_HH
