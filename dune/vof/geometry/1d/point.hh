#ifndef DUNE_VOF_GEOMETRY_1D_POINT_HH
#define DUNE_VOF_GEOMETRY_1D_POINT_HH

#include <cassert>
#include <utility>

namespace Dune {

  namespace VoF {

    template< class > class Point;

    /**
     * \ingroup geo1d
     * \brief point (with correct members)
     *
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    class Point
    {
      using This = Point< Coord >;

    public:
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;

      static constexpr int dimension = 1;
      static constexpr int dimensionworld = Coordinate::dimension;

      Point () = default;

      Point ( const Coordinate &vertex ) : vertex_( std::move( vertex ) ) {}

      bool operator== ( const This &other) const
      {
        return vertex_ == other.vertex( 0 );
      }

      /**
       * \brief i-th vertex
       */
      const Coordinate& vertex ( int i = 0 ) const
      {
        assert( i == 0 );
        return vertex_;
      }

      Coordinate& vertex ( int i = 0 )
      {
        assert( i == 0 );
        return vertex_;
      }

      /**
       * \brief signed volume
       */
      ctype volume () const
      {
        return 1.0;
      }

      /**
       * \brief centroid
       */
      Coordinate centroid () const
      {
        return vertex_;
      }

      // Implementation defined
      int size () const { return 1; }

    private:
      Coordinate vertex_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_1D_POINT_HH
