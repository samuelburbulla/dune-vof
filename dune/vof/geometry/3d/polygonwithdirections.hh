#ifndef DUNE_VOF_GEOMETRY_3D_POLYGONWITHDIRECTIONS_HH
#define DUNE_VOF_GEOMETRY_3D_POLYGONWITHDIRECTIONS_HH

/* c++ includes */
#include <algorithm>
#include <vector>
#include <limits>

/* dune includes */
#include "dune/common/fvector.hh"

/* local includes */
#include "polyhedron.hh"

namespace Dune {

  namespace VoF {

    template < class ctype >
    Dune::FieldVector< ctype, 2 > project ( Dune::FieldVector< ctype, 3 > v )
    {
      return Dune::FieldVector< ctype, 2 > { v[0], v[1] };
    }


    template < class Polyhedron >
    class PolygonWithDirections
    {
    public:
      using Coordinate = typename Dune::FieldVector< double, 2 >;
      using E = typename Dune::VoF::__impl::Edge< PolygonWithDirections > ;
      using Coordinate3D = typename Polyhedron::Coordinate;
      using Direction = typename Polyhedron::E;
      using Face = typename Polyhedron::F;

      PolygonWithDirections ( const Polyhedron& parent ) : parent_ ( parent ) {}

      const std::vector< Coordinate >& nodes() const { return nodes_; }

      const Coordinate& node ( std::size_t index ) const { return nodes_[ index ]; }

      const std::vector< std::size_t >& nodeIds() const { return nodeIds_; }

      const std::vector< E >& edges () const { return edges_; }

      const E& edge ( const std::size_t index ) const { return edges_[ index ]; }

      const std::vector< Face >& correspondingFaces () const { return correspondingFaces_; }

      const Face& correspondingFace ( const std::size_t index ) const { return correspondingFaces_[ index ]; }

      const std::vector< std::vector< Direction > >& directions() const { return directions_; }

      const std::vector< Direction >& directions( const std::size_t n ) const { return directions_[ n ]; }


      double volume() const
      {
        double vol = 0;
        for ( const auto& edge : edges() )
          vol += edge.center() * edge.outerNormal() * edge.volume();
        return std::abs( vol / 2.0 );
      }

      Coordinate3D directionVector ( const std::size_t n, const std::size_t k ) const
      {
        assert ( n < directions_.size() );
        assert ( k < directions_[ n ].size() );
        Coordinate3D directionVector = directions_[ n ][ k ].node( 1 );
        directionVector -= directions_[ n ][ k ].node( 0 );
        const double norm = project( directionVector ).two_norm();
        if ( norm > 0 ) directionVector /= norm;     // TODO: should not happen, prechoose directionVectors
        return directionVector;
      }

      void initialize ( const std::vector< std::size_t >& ids, const std::vector< double >& d )
      {
        nodes_.clear();
        edges_.clear();
        nodeIds_.clear();
        directions_.clear();
        correspondingFaces_.clear();

        std::size_t N = 1;
        for ( ; d[ N ] == d[ 0 ]; N++ );

        switch( N )
        {
        // Case 1: One single point forms the polygon.
        case 1:
        {
          nodes_.push_back( project( parent_.node( ids[ 0 ] ) ) );
          nodeIds_.push_back( ids[ 0 ] );

          auto directions = parent_.attachedEdgesToNode( ids[ 0 ] );
          std::sort( directions.begin(), directions.end(), *this );
          directions_.push_back( directions );
          break;
        }
        // Case 2: An edge forms the polygon.
        case 2:
        {
          nodes_.push_back( project( parent_.node( ids[ 0 ] ) ) );
          nodes_.push_back( project( parent_.node( ids[ 1 ] ) ) );
          nodeIds_.push_back( ids[ 0 ] );
          nodeIds_.push_back( ids[ 1 ] );

          std::vector< Direction > directions1 ( parent_.attachedEdgesToNode( ids[ 0 ] ) );
          std::vector< Direction > directions2 ( parent_.attachedEdgesToNode( ids[ 1 ] ) );

          for ( std::size_t i = 0; i < directions1.size(); ++i )
            for ( std::size_t j = 0; j < directions2.size(); ++j )
              if ( directions1[ i ].nodeId(0) == directions2[ j ].nodeId(1)
                && directions1[ i ].nodeId(1) == directions2[ j ].nodeId(0) )
              {
                directions1.erase( directions1.begin() + i );
                directions2.erase( directions2.begin() + j );
              }

          std::sort( directions1.begin(), directions1.end(), *this );
          std::sort( directions2.begin(), directions2.end(), *this );

          directions_.push_back( directions1 );
          directions_.push_back( directions2 );
          break;
        }
        // Case 3: A face forms the polygon.
        default:
        {
          for ( const Face& face : parent_.faces() )
          {
            bool allFound = true;
            for ( const std::size_t& nodeId : face.nodeIds() )
              allFound = allFound && ( std::find( ids.begin(), ids.begin() + N, nodeId ) != ids.begin() + N );

            if ( allFound )
            {
              for ( const std::size_t nodeId : face.nodeIds() )
              {
                nodeIds_.push_back( nodeId );
                nodes_.push_back( project( parent_.node( nodeId ) ) );
              }
              break;
            }
          }

          for ( std::size_t k = 0; k < N; ++k )
          {
            std::vector< Direction > directions = parent_.attachedEdgesToNode( nodeIds_[ k ] );

            for ( std::size_t j = 0; j < N; ++j )
              for ( std::size_t i = 0; i < directions.size(); ++i )
                if ( directions[ i ].nodeId( 1 ) == nodeIds_[ j ] )
                {
                  directions.erase( directions.begin() + i );
                  --i;
                }

            std::sort( directions.begin(), directions.end(), *this );
            directions_.push_back( directions );
          }
          break;
        }
        } // end switch

        // Build edges and find corresponding faces.
        if ( nodeIds_.size() != 1 )
          for ( std::size_t i = 0; i < nodeIds_.size(); ++i )
          {
            edges_.push_back( E ( this, {{ i, (i+1) % nodeIds_.size() }} ) );
            correspondingFaces_.push_back( parent_.attachedFaceToEdge ( Direction ( &parent_, {{ nodeIds_[ (i+1) % nodeIds_.size() ], nodeIds_[ i ] }} ) ) );
          }

      }


      void evolveToNextPolygon ( const std::size_t k, const double h, const std::vector< double >& dUnique )
      {
        std::vector< Coordinate > newNodes;
        std::vector< std::size_t > newNodeIds;
        std::vector< std::vector< Direction > > newDirections;
        correspondingFaces_.clear();
        edges_.clear();

        for ( std::size_t n = 0; n < nodes().size(); ++n )
          for ( std::size_t i = 0; i < directions( n ).size(); ++i )
          {
            const Direction& direction = directions( n )[ i ];

            if ( std::abs( direction.node( 1 )[ 2 ] - dUnique[ k+1 ] ) < 1e-10 )
            {
              const std::size_t nodeId = direction.nodeId( 1 );

              // If multiple edges converge to the same node, only store one of them
              // ===================================================================
              if ( newNodeIds.size() > 0 )
              {
                if ( nodeId == newNodeIds[ newNodeIds.size()-1 ] )
                {
                  correspondingFaces_.erase( correspondingFaces_.end() - 1 );
                  correspondingFaces_.push_back( parent_.attachedFaceToEdge( direction ) );
                  continue;
                }
              }

              // Add new node from parent and attached information
              // =================================================
              newNodes.push_back( project( parent_.node( nodeId ) ) );
              newNodeIds.push_back( nodeId );

              std::vector< Direction > dirs;
              for ( const auto& edge : parent_.attachedEdgesToNode( nodeId ) )
                if ( edge.node( 1 )[ 2 ] > dUnique[ k ] )
                  dirs.push_back( edge );

              std::sort( dirs.begin(), dirs.end(), *this );
              newDirections.push_back( dirs );

              correspondingFaces_.push_back( parent_.attachedFaceToEdge( direction ) );
            }
            else if ( direction.node( 1 )[ 2 ] > dUnique[ k+1 ] )
            {
              // Add new built node and attached information
              // ===========================================
              const Coordinate3D dv = directionVector( n, i );

              Coordinate dv2D = project( dv );
              dv2D *= h / dv[2];

              Coordinate newNode = node( n );
              newNode += dv2D;

              newNodes.push_back( newNode );
              newNodeIds.push_back( -1 );
              newDirections.push_back( { direction } );

              correspondingFaces_.push_back( parent_.attachedFaceToEdge( direction ) );
            }
          }

        if ( newNodeIds[ 0 ] == newNodeIds[ newNodeIds.size()-1 ] && newNodeIds[ 0 ] != std::size_t( -1 ) )
        {
          newNodes.erase( newNodes.begin() + newNodes.size() - 1 );
          newNodeIds.erase( newNodeIds.begin() + newNodeIds.size() - 1 );
          newDirections.erase( newDirections.begin() + newDirections.size() - 1 );
          correspondingFaces_.erase( correspondingFaces_.begin() + correspondingFaces_.size() - 1 );
        }

        for ( auto &dirs : newDirections )
          for ( std::size_t i = 0; i < dirs.size(); ++i )
            for ( const std::size_t nodeId : newNodeIds )
              if ( dirs[ i ].nodeId( 1 ) == nodeId )
              {
                dirs.erase( dirs.begin() + i );
                --i;
              }

        // Build edges.
        if ( newNodes.size() != 1 )
          for ( std::size_t i = 0; i < newNodes.size(); ++i )
            edges_.push_back( E ( this, {{ i, (i+1) % newNodes.size() }} ) );


        // Fix issue with corresponding face, if two nodes of parent lie next to each other and form an existing edge. Then choose next face.
        if ( newNodeIds.size() > 2 )
          for ( std::size_t i = 0; i < edges_.size(); ++i )
          {
            const E& edge = edges_[ i ];
            std::size_t id0 = newNodeIds[ edge.nodeId(0) ];
            std::size_t id1 = newNodeIds[ edge.nodeId(1) ];

            if ( id0 != std::size_t(-1) && id1 != std::size_t(-1) )
            {
              const std::vector< Direction >& attachedEdges = parent_.attachedEdgesToNode( id0 );
              const auto pos = std::find( attachedEdges.begin(), attachedEdges.end(), Direction( &parent_, {{ id0, id1 }} ) );
              if ( pos != attachedEdges.end() )
              {
                correspondingFaces_.erase( correspondingFaces_.begin() + i );
                correspondingFaces_.insert( correspondingFaces_.begin() + i, parent_.attachedFaceToEdge( Direction( &parent_, {{ id1, id0 }} ) ) );
              }
            }
          }

        assert( !directions_[ 0 ].empty() );

        nodes_ = newNodes;
        nodeIds_ = newNodeIds;
        directions_ = newDirections;


      };

      // Comparator for sorting directions.
      bool operator() ( Direction& lhs, Direction& rhs )
      {
         for ( const auto& edge : parent_.attachedFaceToEdge( lhs ).edges() )
          if ( edge.nodeId(0) == rhs.nodeId(1) && edge.nodeId(1) == rhs.nodeId(0) )
            return true;

        for ( const auto& edge : parent_.attachedFaceToEdge( rhs ).edges() )
          if ( edge.nodeId(0) == lhs.nodeId(1) && edge.nodeId(1) == lhs.nodeId(0) )
            return false;

        assert( false );
        return false;
      }


    private:
      const Polyhedron parent_;
      std::vector< Coordinate > nodes_;
      std::vector< E > edges_;
      std::vector< std::size_t > nodeIds_;
      std::vector< std::vector< Direction > > directions_;
      std::vector< Face > correspondingFaces_;
    };

  } // namespace VoF

} // namespace Dune

#endif // DUNE_VOF_GEOMETRY_3D_POLYGONWITHDIRECTIONS_HH
