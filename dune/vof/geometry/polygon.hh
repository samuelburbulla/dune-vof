#ifndef DUNE_VOF_GEOMETRY_POLYGON_HH
#define DUNE_VOF_GEOMETRY_POLYGON_HH

#include <cassert>

#include <array>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>


namespace Dune {

  namespace VoF {

    template< class > class Polygon;
    template< class > class Line;


    template< class Coord >
    class Polygon
    {
      using This = Polygon< Coord >;
      static_assert(Coord::dimension == 2, "Dimension must be == 2." );

    public:
      using Coordinate = Coord;
      using Edge = Line< Coordinate >;
      using ctype = typename Coordinate::value_type;
      using Container = std::vector< Coordinate >;

      static constexpr int dimension = 2;
      static constexpr int dimensionworld = Coordinate::dimension;
      static_assert( dimension <= dimensionworld, "dimension larger than dimensionworld." );

      Polygon () = default;

      Polygon ( Container vertices )
      : vertices_( std::move( vertices ) )
      {
        assert( vertices_.size() > 0 );
      }

      const Coordinate& vertex ( int i ) const
      {
        assert( i < size( 2 ) );
        return vertices()[ i ];
      }

      Edge edge ( int i ) const
      {
        assert( i < size( 1 ) );
        return Edge{ vertex( i ), vertex( (i+1)%size( 1 ) ) };
      }

      ctype volume () const
      {
        double sum = 0;
        int n = size( 1 );
        for( int i = 0; i < n; i++ )
          sum += ( vertex( i )[ 1 ] + vertex( (i+1)%n )[ 1 ] ) * ( vertex( i )[ 0 ] - vertex( (i+1)%n )[ 0 ] );
        return sum / 2.0;
      }

      Coordinate centroid () const
      {
        DUNE_THROW( NotImplemented, "Polygon.centroid() not yet implemented." );
      }

      // Implementation defined

      int size ( int codim = 2 ) const { return (codim == 0) ? 1 : vertices().size(); }

    private:
      const Container& vertices ()  const { return vertices_; }

      Container vertices_;
    };


    template< class Coord >
    class Line
    {
      using This = Line< Coord >;
      static_assert(Coord::dimension == 2, "Dimension must be == 2." );

    public:
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;
      using Container = std::vector< Coordinate >;

      static constexpr int dimension = 1;
      static constexpr int dimensionworld = Coordinate::dimension;
      static_assert( dimension <= dimensionworld, "dimension larger than dimensionworld." );

      Line ()
      {}

      Line ( Container vertices )
      : vertices_( std::move( vertices ) )
      {
        assert( vertices_.size() == 2 );
      };

      Line ( const Coordinate& v0, const Coordinate& v1 )
      : vertices_{ v0, v1 }
      {};

      const Coordinate& vertex ( int i ) const
      {
        assert( i < size( 1 ) );
        return vertices()[ i ];
      }

      ctype volume () const
      {
        return ( vertex( 0 ) - vertex( 1 ) ).two_norm();
      }

      Coordinate centroid () const
      {
        auto c = vertex( 0 ) + vertex( 1 );
        c *= 0.5;
        return c;
      }

      // Implementation defined

      int size ( int codim = 1 ) const { return ( codim == 0 ) ? 1 : vertices().size() ; }

    private:
      const Container& vertices ()  const { return vertices_; }

      Container vertices_;
    };


    template< class Coord >
    static inline auto make_polygon( std::vector< Coord > vertices ) -> Polygon< Coord >
    {
      return Polygon< Coord >( std::move( vertices) );
    }

    template< class Geometry >
    static inline auto make_polygon( const Geometry& geometry ) -> Polygon< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if( type.isTriangle() )
        return make_polygon( Container{ geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 2 ) } );
      else if ( type.isQuadrilateral() )
        return make_polygon( Container{ geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 3 ), geometry.corner( 2 ) } );
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYGON_HH
