#ifndef DUNE_VOF_GEOMETRY_3D_FACE_HH
#define DUNE_VOF_GEOMETRY_3D_FACE_HH

/* c++ includes */
#include <cmath>
#include <vector>


namespace Dune {

  namespace VoF {

    template < class Coord >
    struct Face
    {
      Face () {};

      Face ( const std::vector< Coord >& nodes ) : nodes_ ( nodes )
      {
        assert( nodes.size() > 2 );
        faceUnitNormal_ = generalizedCrossProduct( vertex( 2 ) - vertex( 1 ), vertex( 0 ) - vertex( 1 ) );
        faceUnitNormal_ /= faceUnitNormal_.two_norm();
      }

      operator bool () const { return !nodes_.empty(); }

      bool operator== ( const Face< Coord > &other) const
      {
        if ( other.size() != this->size() )
          return false;

        std::size_t i0 = 0;
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

      const Coord& vertex ( const std::size_t index ) const { return nodes_[ index % size() ]; }


      double triangleVolume( const std::vector< Coord > &n ) const
      {
        assert ( n.size() == 3 );
        double sum = 0;
        for( std::size_t i = 0; i < 3; ++i )
        {
          Coord center = n[ (i+1)%3 ] ;
          center += n[ i ];
          center *= 0.5;

          Coord edge = n[ (i+1)%3 ];
          edge -= n[ i ];
          sum += center * generalizedCrossProduct( edge, faceUnitNormal_ );
        }
        return std::abs( sum ) / 2.0;
      }

      Coord centroid () const
      {
        Coord c ( 0.0 );
        for ( std::size_t i = 0; i < size(); ++i )
          c += nodes_[ i ];
        c /= static_cast< double >( size() );

        Coord centroid ( 0.0 );

        for( std::size_t i = 0; i < size(); ++i )
        {
          double volTriang = triangleVolume( { vertex( i ), vertex( i+1 ), c } );
          Coord cTriang = c;
          cTriang += vertex( i );
          cTriang += vertex( i+1 );
          cTriang /= 3.0;
          centroid.axpy( volTriang, cTriang );
        }

        centroid /= volume();
        return centroid;
      }

      double volume () const
      {
        double sum = 0;
        for( std::size_t i = 0; i < size(); ++i )
          sum += edgeCenter( i ) * edgeOuterNormal( i );
        return std::abs( sum ) / 2.0;
      }

    private:
      const Coord edgeOuterNormal( const std::size_t i ) const
      {
        Coord edge = vertex( i+1 );
        edge -= vertex( i );
        return generalizedCrossProduct( edge, faceUnitNormal_ );
      }

      const Coord edgeCenter( const std::size_t i ) const
      {
        Coord center = vertex( i+1 );
        center += vertex( i );
        center *= 0.5;
        return center;
      }

      const std::vector< Coord > nodes_;
      Coord faceUnitNormal_;
    };

  } // namespace VoF

} // namespace Dune

#endif // DUNE_VOF_GEOMETRY_3D_FACE_HH
