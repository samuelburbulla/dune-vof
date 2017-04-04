#ifndef ELLIPSE_HH
#define ELLIPSE_HH

#include <cmath>
#include <array>
#include <numeric>

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
      const auto tmp = ( axis_[ i ] * ( x - DomainType( 0.5 ) ) ) / radii_[ i ];
      val += tmp * tmp;
    }

    u = val < 1.0 ? 1.0 : 0.0;
  }

  DomainType center () const { return DomainType( 0.5 ); }

  decltype(auto) referenceMap () const
  {
    return [ this ] ( const DomainType& x )
    {
      DomainType y;
      for ( int i = 0; i < dimDomain; ++i )
        y[ i ] = ( this->axis_[ i ] * ( x - this->center() ) ) / this->radii_[ i ];
      return y;
    };
  }

  double volumeElement () const
  {
    return std::accumulate( radii_.begin(), radii_.end(), 1.0, std::multiplies<>() );
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
    double radius = 0.0;
    radius += std::pow( axis_[ 0 ] * ( x - DomainType( 0.5 ) ), 2 ) * std::pow( radii_[ 1 ], 4 );
    radius += std::pow( axis_[ 1 ] * ( x - DomainType( 0.5 ) ), 2 ) * std::pow( radii_[ 0 ], 4 );
    radius = std::pow( radius, 3 );

    return std::pow( radii_[ 0 ] * radii_[ 1 ], 4 ) / std::sqrt( radius );
  }

  std::array< DomainType, dim > axis_;
  std::array< double, dim > radii_;
};

#endif // #ifndef ELLIPSE_HH
