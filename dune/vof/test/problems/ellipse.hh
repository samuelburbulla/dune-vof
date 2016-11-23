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
      val += tmp * tmp / radii_[ i ];
    }

    u = val < 1.0 ? 1.0 : 0.0;
  }

  void evaluate ( const DomainType& x, double t, RangeType& u ) const
  {
    evaluate( x, u );
  }

  void velocityField ( const DomainType&, double t, DomainType& v ) const
  {
    v = 0.0;
  }

  std::array< DomainType, dim > axis_;
  std::array< double, dim > radii_;
};

#endif // #ifndef ELLIPSE_HH
