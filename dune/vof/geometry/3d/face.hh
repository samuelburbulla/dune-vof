#ifndef DUNE_VOF_GEOMETRY_FACE_HH
#define DUNE_VOF_GEOMETRY_FACE_HH

/* c++ includes */
#include <vector>


namespace Dune {

  namespace VoF {

    template < class Coord >
    struct Face
    {
      Face () = default;

      Face ( const std::vector< Coord >& nodes ) : nodes_ ( nodes ) {}

      std::size_t size () const { return nodes_.size(); }

      const Coord& vertex ( const std::size_t index ) const { return nodes_[ index ]; }

      Coord centroid () const
      {
        Coord center ( 0.0 );
        for ( std::size_t i = 0; i < nodes_.size(); ++i )
          center += nodes_[ i ];
        center *= 1.0 / nodes_.size();
        return center;
      }

    private:
      const std::vector< Coord > nodes_;
    };

  } // namespace VoF

} // namespace Dune

#endif // DUNE_VOF_GEOMETRY_FACE_HH
