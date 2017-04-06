#ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
#define DUNE_VOF_GEOMETRY_ALGORITHM_HH

#include <cassert>

#include <functional>
#include <limits>
#include <vector>

#include <dune/vof/brents.hh>
#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/3d/polygonwithdirections.hh>
#include <dune/vof/geometry/3d/rotation.hh>


namespace Dune {

  namespace VoF {


    /**
     * \ingroup Geometry
     * \brief calculate fraction of the volume of the intersection of a half space with a polygon
     *
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    double getVolumeFraction ( const Polygon< Coord > &polygon, const HalfSpace< Coord > &halfSpace )
    {
      Polygon< Coord > intersection = intersect( std::cref( polygon ), std::cref( halfSpace ) );

      return intersection.volume() / polygon.volume();
    }

    template< class Coord >
    double getVolumeFraction ( const Polyhedron< Coord > &polyhedron, const HalfSpace< Coord > &halfSpace )
    {
      Polyhedron< Coord > intersection = intersect( std::cref( polyhedron ), std::cref( halfSpace ) );

      return intersection.volume() / polyhedron.volume();
    }


    // locateHalfSpace
    // ---------------

    /**
     * \ingroup Geometry
     * \brief locate half space with a given normal that intersects a polygon so it has a given volume fraction
     *
     * \param polygon   polygon
     * \param normal    normal vector
     * \param fraction  volume fraction
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    auto locateHalfSpace ( const Polygon< Coord >& polygon, const Coord& normal, double fraction ) -> HalfSpace< Coord >
    {
      using ctype = typename HalfSpace< Coord >::ctype;
      using limits = std::numeric_limits< ctype >;

      ctype pMin = 0, pMax = 0;
      ctype volume,
            volMin = limits::lowest(), // lower bound
            volMax = limits::max(); // upper bound


      // upper/lower bound for
      for( int i = 0; i < polygon.size(); ++i )
      {
        ctype dist = -( normal * polygon.vertex( i ) );
        HalfSpace< Coord > hs( normal, dist );

        volume = getVolumeFraction( polygon, hs );

        if( ( volume <= volMax ) && ( volume >= fraction ) )
        {
          pMax = dist;
          volMax = volume;
        }
        if( ( volume >= volMin ) && ( volume <= fraction ) )
        {
          pMin = dist;
          volMin = volume;
        }
      }
      assert( volMin <= volMax );

      if ( volMax == limits::max() )
        return HalfSpace< Coord >( normal, pMin );
      else if ( volMin == limits::lowest() )
        return HalfSpace< Coord >( normal, pMax );
      else if ( volMin == volMax )
        return HalfSpace< Coord >( normal, pMin );
      else
        return HalfSpace< Coord >( normal, brentsMethod( [ &polygon, &normal, &fraction ] ( ctype p ) -> ctype {
                                                            HalfSpace< Coord > hs( normal, p );
                                                            return ( getVolumeFraction( polygon, hs ) - fraction );
                                                          }, pMin, pMax, 1e-12 ) );
    }






    template< class Coord >
    auto locateHalfSpace ( const Polyhedron< Coord >& cell, const Coord& innerNormal, double fraction ) -> HalfSpace< Coord >
    {
      using std::abs;
      assert( abs( innerNormal.two_norm() - static_cast< double > ( 1.0 ) ) < std::numeric_limits< double >::epsilon() );

      Coord outerNormal ( innerNormal );
      outerNormal *= -1.0;

      fraction = ( fraction > 1.0 ) ? 1.0 : fraction;
      fraction = ( fraction < 0.0 ) ? 0.0 : fraction;


      bool inverseMode = false;
      /*
      if ( fraction > 0.6 )
      {
        fraction = 1.0 - fraction;
        outerNormal *= -1.0;
        inverseMode = true;
      }
      */

      double givenVolume = fraction * cell.volume();

      Polyhedron< Coord > polyhedron ( rotateToReferenceFrame ( outerNormal, cell ) );

      // Compute and sort plane constants
      // --------------------------------
      const std::size_t N = polyhedron.nodes().size();
      std::vector< double > d ( N );
      std::vector< std::size_t > sortedNodesIds ( N );

      for ( std::size_t i = 0; i < N; ++i )
      {
        d[ i ] = polyhedron.node( i )[ 2 ];
        sortedNodesIds[ i ] = i;
      }

      std::sort( sortedNodesIds.begin(), sortedNodesIds.end(), [ &d ]( const std::size_t i, const std::size_t j ) { return d[ i ] < d[ j ]; } );
      std::sort( d.begin(), d.end() );

      std::vector< double > dUnique ( d );
      auto it = std::unique( dUnique.begin(), dUnique.end(), []( const double x, const double y ) { return std::abs( x - y ) < 1e-10; } );
      dUnique.resize( std::distance( dUnique.begin(), it ) );


      PolygonWithDirections< Polyhedron< Coord > > polygon ( polyhedron );
      polygon.initialize( sortedNodesIds, d );


      // Bracketing
      // ----------
      double Vd_k = 0, Vd_k1;
      for ( std::size_t k = 0; k < dUnique.size() - 1; ++k )
      {

        // Compute cooefficients
        const std::size_t Nn = polygon.nodes().size();
        double A = 0.0;
        for ( std::size_t n = 0; n < Nn; ++n )
        {
          assert ( !polygon.directions( n ).empty() );

          const std::size_t epsn = polygon.directions( n ).size() - 1;

          Coord v1 = polygon.directionVector( n, epsn );
          Coord v2 = polygon.directionVector( (n+1) % Nn, 0 );
          A -= ( v1[0] * v2[1] - v1[1] * v2[0] ) / ( v1[2] * v2[2] );

          for ( std::size_t l = 0; l < epsn; ++l )
          {
            Coord w1 = polygon.directionVector( n, l );
            Coord w2 = polygon.directionVector( n, l+1 );
            A -= ( w1[0] * w2[1] - w1[1] * w2[0] ) / ( w1[2] * w2[2] );
          }
        }
        assert ( A == A );

        double B = 0.0;
        for ( std::size_t i = 0; i < polygon.edges().size(); ++i )
        {
          Coord normal = polygon.correspondingFace( i ).outerNormal();
          double norm = project( normal ).two_norm();

          if ( norm > 0 )
            B -= polygon.edge( i ).volume() * normal[2] / norm;
        }
        assert ( B == B );

        double C = polygon.volume();
        assert ( C == C );

        double h = dUnique[ k+1 ] - dUnique[ k ];
        auto V = [ A, B, C ] ( const double h ) { return C * h + B * h * h / 2.0 + A * h * h * h / 6.0; };
        Vd_k1 = Vd_k + V( h );

        if ( std::abs( Vd_k1 - givenVolume ) < 1e-14 )
        {
          double distance = inverseMode ? -dUnique[ k+1 ] : dUnique[ k+1 ];
          return HalfSpace< Coord > ( innerNormal, distance );
        }

        else if ( Vd_k1 > givenVolume )
        {
          const auto& Vh_V = [ &V, givenVolume, Vd_k ] ( const double h ) { return V( h ) - ( givenVolume - Vd_k ); };
          double distance = dUnique[ k ] + Dune::VoF::brentsMethod ( Vh_V, 0.0, h, 1e-14 );
          distance = inverseMode ? -distance : distance;
          return HalfSpace< Coord > ( innerNormal, distance );
        }

        else
          polygon.evolveToNextPolygon( k, h, dUnique );

        Vd_k = Vd_k1;
      }

      return HalfSpace< Coord > ( innerNormal, inverseMode ? -dUnique[ dUnique.size() - 1 ] : dUnique[ dUnique.size() - 1 ] );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
