#ifndef DUNE_VOF_GEOMETRY_FACE_HH
#define DUNE_VOF_GEOMETRY_FACE_HH

/* c++ includes */
#include <cmath>
#include <vector>


namespace Dune {

  namespace VoF {

    template < class Coord >
    struct Face
    {
      Face () {};

      Face ( const std::vector< Coord >& nodes ) : nodes_ ( nodes ) {}

      operator bool () const { return !nodes_.empty(); }

      bool operator== ( const Face< Coord > &other) const
      {
        if ( other.size() != this->size() )
          return false;

        std::size_t i0;
        for ( std::size_t i = 0; i < other.size(); ++i )
          if ( other.vertex( i ) == this->vertex( 0 ) )
          {
            i0 = i;
            break;
          }

        for ( std::size_t i = 0; i < other.size(); ++i )
          if ( other.vertex( (i + i0)%other.size() ) != this->vertex( i ) )
            return false;

        return true;
      }

      std::size_t size () const { return nodes_.size(); }

      const Coord& vertex ( const std::size_t index ) const { return nodes_[ index ]; }

      Coord centroid () const
      {
        Coord center ( 0.0 );
        for ( std::size_t i = 0; i < nodes_.size(); ++i )
          center += nodes_[ i ];
        center *= 1.0 / ( nodes_.size() );
        return center;
      }

      double volume () const
      {
        double sum = 0;
        Coord normal = generalizedCrossProduct( vertex( 2 ) - vertex( 1 ), vertex( 0 ) - vertex( 1 ) );
        normal /= normal.two_norm();

        int n = size();
        for( int i = 0; i < n; i++ )
        {
          Coord edge = vertex( (i+1)%n ) - vertex( i );
          Coord center = vertex( (i+1)%n ) + vertex( i );
          center *= 0.5;

          Coord outerNormal = generalizedCrossProduct( edge, normal );
          outerNormal /= outerNormal.two_norm();

          sum += edge.two_norm() * ( center * outerNormal );
        }
        assert( !std::isnan( sum ) );
        return std::abs( sum ) / 2.0;
      }

    private:
      const std::vector< Coord > nodes_;
    };

  } // namespace VoF

} // namespace Dune

#endif // DUNE_VOF_GEOMETRY_FACE_HH
