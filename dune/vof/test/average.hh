# ifndef DUNE_VOF_AVERAGE_HH
# define DUNE_VOF_AVERAGE_HH

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

#include "interpolation.hh"
#include "recursiveinterpolation.hh"
#include "problems/ellipse.hh"
#include "problems/rotatingcircle.hh"
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
    void averageRecursive ( DF &u, const F &f, const double time = 0.0 )
    {
      if( u.gridView().comm().rank() == 0 )
        std::cout << " -- average using recursive algorithm" << std::endl;

      RecursiveInterpolation< typename std::remove_reference< decltype( u.gridView() ) >::type > interpolation ( u.gridView(), 5 );
      auto function = [ &f, time ]( const auto &x ) { FieldVector< double, 1 > u; f.evaluate( x, u ); return u; };
      interpolation( function, u );
    }


    template< class DF, class F >
    void average ( DF &u, const F &f, const double time = 0.0, const double x = 0.0 )
    {
      typedef typename DF::GridView::ctype ctype;
      using RangeType = Dune::FieldVector< ctype, 1 >;

      if( u.gridView().comm().rank() == 0 )
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
    void average ( DF &u, const Ellipse< double, 2 >& e, const double time = 0.0 )
    {
      if( u.gridView().comm().rank() == 0 )
        std::cout << " -- average using intersection" << std::endl;
      circleInterpolation( e.referenceMap(), e.volumeElement(), u );
    }

    template< class DF >
    void average ( DF &u, const RotatingCircle< double, 2 >& c, const double time = 0.0 )
    {
      if( u.gridView().comm().rank() == 0 )
        std::cout << " -- average using intersection" << std::endl;
      circleInterpolation( c.center( time ), c.radius( time ), u );
    }

    template< class DF >
    void average ( DF &uh, const Slope< double, 2 >& s, const double time = 0.0 )
    {
      using Coordinate = FieldVector< double, 2 >;
      if( uh.gridView().comm().rank() == 0 )
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
