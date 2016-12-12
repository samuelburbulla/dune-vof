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


    /**
     * \ingroup geo2d
     * \brief (oriented) polygon
     *
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    class Polygon
    {
      using This = Polygon< Coord >;

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

      /**
       * \brief i-th vertex
       */
      const Coordinate& vertex ( int i ) const
      {
        assert( i < size( 2 ) );
        return vertices()[ i ];
      }

      /**
       * \brief i-th edge
       */
      Edge edge ( int i ) const
      {
        assert( i < size( 1 ) );
        return Edge{ vertex( i ), vertex( (i+1)%size( 1 ) ) };
      }

      /**
       * \brief signed volume
       */
      ctype volume () const
      {
        double sum = 0;
        int n = size( 1 );
        for( int i = 0; i < n; i++ )
          sum += ( vertex( i )[ 1 ] + vertex( (i+1)%n )[ 1 ] ) * ( vertex( i )[ 0 ] - vertex( (i+1)%n )[ 0 ] );
        return sum / 2.0;
      }

      /**
       * \brief centroid
       */
      Coordinate centroid () const
      {
        //DUNE_THROW( NotImplemented, "Polygon.centroid() not yet implemented." );
        Coordinate c ( 0 );
        for( int i = 0; i < size(); i++ )
          c += vertex( i );
        c /= size();
        return c;
      }

      // Implementation defined

      int size ( int codim = 2 ) const { return (codim == 0) ? 1 : vertices().size(); }

    private:
      const Container& vertices ()  const { return vertices_; }

      Container vertices_;
    };

    /**
     * \ingroup geo2d
     * \brief (oriented) line
     *
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    class Line
    {
      using This = Line< Coord >;

    public:
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;
      using Container = std::array< Coordinate, 2 >;

      static constexpr int dimension = 1;
      static constexpr int dimensionworld = Coordinate::dimension;
      static_assert( dimension <= dimensionworld, "dimension larger than dimensionworld." );

      Line ()
        : vertices_{ Coordinate( 0 ), Coordinate( 0 ) }
      {}

      Line ( Container vertices )
      : vertices_( std::move( vertices ) )
      {}

      Line ( const Coordinate& v0, const Coordinate& v1 )
      : vertices_{ v0, v1 }
      {}

      bool operator== ( const This &other) const
      {
        return vertices_ == other.vertices_;
      }

      const Coordinate& vertex ( std::size_t i ) const
      {
        assert( i < size( 1 ) );
        return vertices()[ i ];
      }

      /**
       * \brief signed volume
       */
      ctype volume () const
      {
        return ( vertex( 0 ) - vertex( 1 ) ).two_norm();
      }

      /**
       * \brief centroid
       */
      Coordinate centroid () const
      {
        auto c = vertex( 0 ) + vertex( 1 );
        c *= 0.5;
        return c;
      }

      // Implementation defined

      std::size_t size ( int codim = 1 ) const { return ( codim == 0 ) ? 1 : 2; }

    private:
      const Container& vertices ()  const { return vertices_; }

      Container vertices_;
    };


    /**
     * \ingroup geo2d
     * \brief generate polygon
     *
     * \param vertices vertices assumed to be in order
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    static inline auto makePolygon( std::vector< Coord > vertices ) -> Polygon< Coord >
    {
      return Polygon< Coord >( std::move( vertices) );
    }

    /**
     * \ingroup geo2d
     * \brief generate polygon
     *
     * \param geometry dune geometry
     */
    template< class Geometry >
    static inline auto makePolygon( const Geometry& geometry ) -> Polygon< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if( type.isTriangle() )
        return makePolygon( Container{ geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 2 ) } );
      else if ( type.isQuadrilateral() )
        return makePolygon( Container{ geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 3 ), geometry.corner( 2 ) } );
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }

    template< class Geometry, class Map >
    static inline auto makePolygon( const Geometry& geometry, Map&& map )
      -> Polygon< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if( type.isTriangle() )
        return makePolygon( Container{ map( geometry.corner( 0 ) ), map( geometry.corner( 1 ) ), map( geometry.corner( 2 ) ) } );
      else if ( type.isQuadrilateral() )
        return makePolygon( Container{ map( geometry.corner( 0 ) ), map( geometry.corner( 1 ) ), map( geometry.corner( 3 ) ), map( geometry.corner( 2 ) ) } );
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYGON_HH
