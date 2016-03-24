#ifndef DUNE_VOF_GEOMETRY_POLYTOP_HH
#define DUNE_VOF_GEOMETRY_POLYTOP_HH

#include <cassert>
#include <vector>

#include <dune/common/deprecated.hh>

/*
 * TODO:
 * - properly implement AlgebraicPolytope / Polytope
 * - ...
 */

namespace Dune {

  namespace VoF {

    // forward declarations
    // --------------------

    template< int > class AlgebraicPolytope;
    template< class, int > class Polytope;


    // AlgebraicPolytope
    // -----------------

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

      int size ( int codim = dim ) const
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


    // Polytope
    // --------

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
      static_assert( dimension <= dimensionworld, "dimension larger than dimensionworld." );

      Polytope ( const Topology& topology, const Container& vertices )
      : topology_( topology ), vertices_( vertices.begin(), vertices.end() )
      {
        assert( size( dim ) == vertices().size() );
      };

      int size ( int codim = dim ) const { return topology().size( codim ); }
      const Topology& topology () const { return topology_; }

      // subTopology, subPolytope ... ?
      // ...

      const Coordinate& vertex ( int i ) const
      {
        assert( i < size( dim ) );
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

    private:
      const Container& vertices ()  const { return vertices_; }

      const Topology& topology_;
      Container vertices_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYTOP_HH
