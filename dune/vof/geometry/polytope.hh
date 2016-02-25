#ifndef DUNE_VOF_GEOMETRY_POLYTOPE_HH
#define DUNE_VOF_GEOMETRY_POLYTOPE_HH

#include <cassert>
#include <vector>

#include <dune/common/deprecated.hh>


namespace Dune {

  namespace VoF {

     // OLD CODE / INTERFACE

     // konvexes 2D-Polygon
    template< class DomainVector >
    struct DUNE_DEPRECATED_MSG( "Use Polytope when it is implemented.") Polygon2D
    {
      Polygon2D< DomainVector >() {}

      const DomainVector operator[] ( const int i ) const { return points[ i % points.size() ]; }

      void addVertex ( const DomainVector &vertex, const bool correctOrder = false, const double TOL = 1e-12 )
      {
        std::size_t n = points.size();

        // no specific order needed or vertex explicitly given in correct order
        if( n == 0 || n == 1 || correctOrder )
        {
          points.push_back( vertex );
        }
        else
        {
          // insert new vertex in counterclockwise order
          for( std::size_t i = 0; i < n; ++i )
          {
            auto normal = points[ (i+1)%n ] - points[ i ];
            std::swap( normal[ 0 ], normal[ 1 ] );
            normal[ 0 ] *= -1.0;

            auto center = points[ i ];
            center += points[ (i+1)%n ];
            center *= 0.5;

            center -= vertex;

            if( normal * center > 0 )
            {
              points.insert( points.begin() + i + 1, vertex );
              break;
            }
          }
        }
      }

      const std::size_t corners () const { return points.size(); }

      double volume () const
      {
        double sum = 0;
        int n = points.size();
        for( int i = 0; i < n; i++ )
          sum += ( points[ i ][ 1 ] + points[ (i+1)%n ][ 1 ] ) * ( points[ i ][ 0 ] - points[ (i+1)%n ][ 0 ] );
        return sum / 2.0;
      }

      void clear () { points.clear(); }

    private:
      std::vector< DomainVector > points;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYTOPE_HH
