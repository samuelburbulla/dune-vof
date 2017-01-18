# ifndef DUNE_VOF_AVERAGE_HH
# define DUNE_VOF_AVERAGE_HH

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/vof/interpolation.hh>

#include "problems/ellipse.hh"
#include "problems/slope.hh"
#include "../geometry/intersect.hh"
#include "../geometry/halfspace.hh"

namespace Dune
{
  namespace VoF
  {

    // Average
    // -------

    template< class DF, class F >
    void average ( DF &u, const F &f, const double time = 0.0 )
    {
      typedef typename DF::GridView::ctype ctype;
      using RangeType = Dune::FieldVector< ctype, 1 >;

      std::cout << " -- average using quadrature" << std::endl;

      for( const auto& entity : elements ( u.gridView() ) )
      {
        const auto geo = entity.geometry();

        // get quadrature rule of order p
        int p = 19;
        const auto &rule = Dune::QuadratureRules< ctype, DF::GridView::dimension >::rule( geo.type(), p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, "Requested quadrature order not available." );

        // compute approximate integral
        ctype result = 0;
        for ( const auto qp : rule )
        {
          RangeType u;
          f.evaluate( geo.global( qp.position() ), time, u );
          result += u * qp.weight() * geo.integrationElement( qp.position() );
        }

        u[ entity ] = result / geo.volume();
      }
    }

    template< class DF >
    void average ( DF &u, const Ellipse< double, 2 >& e )
    {
      std::cout << " -- average using intersection" << std::endl;
      circleInterpolation( e.referenceMap(), e.volumeElement(), u );
    }

    template< class DF >
    void average ( DF &uh, const Slope< double, 2 >& s )
    {
      using Coordinate = FieldVector< double, 2 >;
      std::cout << " -- average using intersection" << std::endl;
      uh.clear();

      HalfSpace< Coordinate > halfspace ( s.normal(), Coordinate( { 0.5, 0.5 } ) );

      for ( const auto& entity : elements( uh.gridView(), Partitions::interior ) )
      {
        const auto& geo = entity.geometry();
        Dune::VoF::Polygon< Coordinate > polygon = Dune::VoF::makePolytope( geo );

        auto it = intersect( polygon, halfspace );
        auto part = static_cast< typename decltype( it )::Result > ( it );

        uh[ entity ] = part.volume() / geo.volume();
      }
    }


  } // namespace VoF

} // namespace Dune

# endif // #ifndef DUNE_VOF_AVERAGE_HH
