#ifndef DUNE_VOF_GEOMETRY_POLYTOP_HH
#define DUNE_VOF_GEOMETRY_POLYTOP_HH

#include <cassert>
#include <vector>

#include <dune/common/deprecated.hh>

// TODO:
//  - concretise interface of AlgebraicPolytope / Polytope
//  - add datat representation and implementation
// (- remove Polygon2d )

namespace Dune {

  namespace VoF {

    // forward declarations
    // --------------------

    template< int > class AlgebraicPolytope;
    template< class, int > class Polytope;


    // AlgebraicPolytop
    // ----------------

    template< int dim >
    class AlgebraicPolytope
    {
      using This = AlgebraicPolytope< dim >;
      static_assert( dim >= 0, "Dimension must be >= 0." );

    public:
      static constexpr int dimension = dim;

      template< int codim >
      struct Codim
      {
        using Polytope = AlgebraicPolytope< dim-codim >;
      };

      // constructors ?
      // ...

      int size ( int codim = 0 ) const
      {
        DUNE_THROW( NotImplemented, "size( ... ) not yet imeplemented." );
        return 0;
      }

      template< int codim >
      const typename Codim< codim >::Polytope& subPolytope ( int i ) const
      {
        assert( i < size( codim ) );
      }

    private:
      // data members ... ?
      // ...
    };


    // Polytop
    // -------

    template< class Coord, int dim >
    class Polytope
    {
      using This = Polytope< Coord, dim >;
      static_assert( dim >= 0, "Dimension must be >= 0." );

    public:
      template< int codim >
      struct Codim {
        using SubTopology = AlgebraicPolytope< dim - codim >;
        using SubPolytope  = Polytope< Coord, dim-codim >;
      };
      using Topology = typename Codim< 0 >::SubTopology;

      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;
      using Container = std::vector< Coordinate >;

      static constexpr int dimension = dim;
      static constexpr int dimensionworld = Coordinate::dimension;

      Polytope ( const Topology& topology, const Container& vertices )
      : topology_( topology ), vertices_( vertices.begin(), vertices.end() )
      {
        assert( size( 0 ) == vertices().size() );
      };

      int size ( int codim = 0 ) const { return topology().size( codim ); }
      const Topology& topology () const { return topology_; }

      // subTopology, subPolytope ... ?
      // ...

      const Coordinate& vertex ( int i ) const
      {
        assert( i < size( 0 ) );
        return vertices_[ i ];
      }

      ctype volume () const
      {
        DUNE_THROW( NotImplemented, "volume() not yet immplented" );
        return 0.0;
      }

      Coordinate centroid () const
      {
        DUNE_THROW( NotImplemented, " cnetroid() Not yet immplented" );
        return 0.0;
      }

      // old interface ...
      // -----------------

      int corners () const DUNE_DEPRECATED_MSG( "Use size( 0 ) instead." ) { return topology().size( 0 ); }

      const Coordinate& operator[] ( int i ) const DUNE_DEPRECATED_MSG( "Use vertex( ... ) and topology information instead." )
      {
        return vertices_[ i % vertices_.size() ];
      }

      void addVertex ( const Coordinate&, bool = false, double = 1e-12 ) DUNE_DEPRECATED_MSG( "Functionality removed from interface." ) {}
      void clear () DUNE_DEPRECATED_MSG( "Functionality removed from interface." ) {}

    private:
      const Container& vertices ()  const { return vertices_; }

      const Topology& topology_;
      Container vertices_;
    };



     // OLD CODE / INTERFACE

     // konvexes 2D-Polygon
    template< class DomainVector >
    struct DUNE_DEPRECATED_MSG( "Use Polytope< ... > when it is implemented." ) Polygon2D
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

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYTOP_HH
