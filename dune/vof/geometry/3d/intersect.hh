#ifndef DUNE_VOF_GEOMETRY_3D_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_3D_INTERSECT_HH

#include <cassert>
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <map>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/3d/face.hh>

namespace Dune {

  namespace VoF {

    /**
     * \brief namespace containing specific implementations
     */
    namespace __impl {


      template< class Coord >
      auto intersect ( const Dune::VoF::Face< Coord >& face, const HalfSpace< Coord >& halfSpace ) -> Dune::VoF::Face< Coord >
      {
        if ( !halfSpace )
          return Dune::VoF::Face< Coord >();

        std::vector< Coord > container;
        container.reserve( face.size() );

        for( std::size_t i = 0; i < face.size(); ++i )
        {
          Line< Coord > edge { face.vertex( i ), face.vertex( i%face.size() ) };

          auto l0 = halfSpace.levelSet( edge.vertex( 0 ) );
          auto l1 = halfSpace.levelSet( edge.vertex( 1 ) );

          if( l0 > 0.0 )
            container.push_back( edge.vertex( 0 ) );

          if ( ( l0 > 0.0 ) ^ ( l1 > 0.0 ) )
          {
            Coord point;
            point.axpy( -l1 / ( l0 - l1 ), edge.vertex( 0 ) );
            point.axpy(  l0 / ( l0 - l1 ), edge.vertex( 1 ) );

            container.push_back( point );
          }
        }

        if ( container.empty() )
          return Dune::VoF::Face< Coord >();
        else
          return Dune::VoF::Face< Coord >( std::move( container ) );
      }

      /**
       * \ingroup geo3d
       * \brief implementation for an intersection between a polyhedron and a half space
       *
       * \tparam  Coord  type of the global coordinate
       * \return  the intersection polygon
       */
      template< class Coord >
      auto intersect ( const Polyhedron< Coord >& polyhedron, const HalfSpace< Coord >& halfSpace ) -> Polyhedron< Coord >
      {
        if ( !halfSpace )
          return Polyhedron< Coord >();

        using Coordinate = Coord;
        using Face = typename std::vector< size_t >;
        using Edge = typename std::array< size_t, 2 >;

        const double eps = std::numeric_limits< double >::epsilon();

        std::vector< Coordinate > nodes;
        std::vector< Edge > edges;
        std::vector< Face > faces;

        std::vector< bool > isInner ( polyhedron.nodes().size(), false );
        std::vector< bool > isBoundary ( polyhedron.nodes().size(), false );
        std::vector< Coordinate > newNodes;
        Face intersectionFace;

        using ReverseEdge = typename std::pair< std::size_t, std::size_t >;
        using MapValueType = typename std::pair< ReverseEdge, std::size_t >;
        std::map< ReverseEdge, std::size_t > newNodesMap;

        std::size_t n = 0, onBoundary = 0;
        std::vector< std::size_t > p;

        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( halfSpace.levelSet( polyhedron.node( i ) ) > - eps )
          {
            isInner[ i ] = true;
            p.push_back( n );
            n++;
            if ( std::abs( halfSpace.levelSet( polyhedron.node( i ) ) ) < eps )
            {
              isBoundary[ i ] = true;
              onBoundary++;
            }
          }
          else
            p.push_back( -1 );

        // Handle trivial cases
        if ( n == 0 )
          return Polyhedron< Coord >();
        if ( n == 1 && onBoundary == 1 )
          return Polyhedron< Coord >( {}, std::vector< std::array< std::size_t, 2 > > {}, { polyhedron.node( p[0] ) } );
        if ( n == 2 && onBoundary == 2 )
          return Polyhedron< Coord >( {}, std::vector< std::array< std::size_t, 2 > > { {{ 0, 1 }}, {{ 1, 0 }} }, { polyhedron.node( p[0] ), polyhedron.node( p[1] ) } );

        // Non-trivial case
        for ( const auto& face : polyhedron.faces() )
        {
          Face newFace;
          std::size_t lastIsId = -1;

          for ( const auto& edge : face.edges() )
          {
            // Edge is an inner edge
            if ( isInner[ edge.nodeId(0) ] && isInner[ edge.nodeId(1) ] )
            {
              edges.push_back( {{ p[ edge.nodeId(0) ], p[ edge.nodeId(1) ] }} );
              newFace.emplace_back( edges.size() - 1 );
            }
            // Edge is intersected by halfspace boundary
            else if ( isInner[ edge.nodeId(0) ] ^ isInner[ edge.nodeId(1) ] )
            {

              // Intersection point is corner 0
              if ( isBoundary[ edge.nodeId(0) ] )
              {
                if ( lastIsId != std::size_t(-1) && lastIsId != p[ edge.nodeId(0) ] )
                {
                  edges.push_back( {{ p[ edge.nodeId(0) ], lastIsId }} );
                  newFace.emplace_back( edges.size() - 1 );

                  edges.push_back( {{ lastIsId, p[ edge.nodeId(0) ] }} );
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIsId = p[ edge.nodeId(0) ];
              }
              // Intersection point is corner 1
              else if ( isBoundary[ edge.nodeId(1) ] )
              {
                if ( lastIsId != std::size_t(-1) && lastIsId != p[ edge.nodeId(1) ] )
                {
                  edges.push_back( {{ lastIsId, p[ edge.nodeId(1) ] }} );
                  newFace.emplace_back( edges.size() - 1 );

                  edges.push_back( {{ p[ edge.nodeId(1) ], lastIsId }});
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIsId = p[ edge.nodeId(1) ];
              }
              else
              {
                // Intersection point is new point
                std::size_t isNodeId;

                ReverseEdge revEdge ( edge.nodeId(1), edge.nodeId(0) );
                auto it = newNodesMap.find( revEdge );

                if ( it == newNodesMap.end() )
                {
                  const Coordinate isNode = edge.intersection( halfSpace.boundary() );
                  newNodes.push_back( isNode );
                  isNodeId = n - 1 + newNodes.size();
                  newNodesMap.insert( MapValueType( { edge.nodeId(0), edge.nodeId(1) }, isNodeId ) );
                }
                else
                  isNodeId = it->second;


                if ( isInner[ edge.nodeId(0) ] )
                {
                  edges.push_back( {{ p[ edge.nodeId(0) ], isNodeId }} );
                  newFace.emplace_back( edges.size() - 1 );

                  if ( lastIsId != std::size_t(-1) && lastIsId != isNodeId )
                  {
                    edges.push_back( {{ isNodeId, lastIsId }} );
                    newFace.emplace_back( edges.size() - 1 );

                    edges.push_back( {{ lastIsId, isNodeId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                }
                else
                {
                  if ( lastIsId != std::size_t(-1) && lastIsId != isNodeId )
                  {
                    edges.push_back( {{ lastIsId, isNodeId }} );
                    newFace.emplace_back( edges.size() - 1 );

                    edges.push_back( {{ isNodeId, lastIsId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                  edges.push_back( {{ isNodeId, p[ edge.nodeId(1) ] }} );
                  newFace.emplace_back( edges.size() - 1 );
                }
                lastIsId = isNodeId;
              }
            }
          }
          if ( newFace.size() > 2 )
            faces.push_back( newFace );
        }

        // Sort subentities and nodes of intersection face
        // ===============================================

        // Erase null-edges
        for ( std::size_t i = 0; i < intersectionFace.size(); ++i )
        {
          const Edge& edge = edges[ intersectionFace[ i ] ];
          if ( edge[ 0 ] == edge[ 1 ] )
            intersectionFace.erase( intersectionFace.begin() + i );
        }

        // Sort edges of intersection face
        for ( std::size_t i = 0; i + 1 < intersectionFace.size(); ++i )
        {
          if ( edges[ intersectionFace[ i ] ][ 1 ] == edges[ intersectionFace[ i+1 ] ][ 0 ] )
            continue;

          const std::size_t index = edges[ intersectionFace[ i ] ][ 1 ];
          const auto pos = std::find_if(
            intersectionFace.begin() + i,
            intersectionFace.end(),
            [ &edges, index ]( const std::size_t& edgeId ) -> bool { return edges[ edgeId ][ 0 ] == index; }
          );

          assert ( pos != intersectionFace.end() );
          if ( pos != intersectionFace.end() )
            std::swap( intersectionFace[ i+1 ], *pos );
        }

        if ( intersectionFace.size() > 2 )
          faces.push_back( intersectionFace );


        // Build new polytop
        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( isInner[ i ] )
            nodes.push_back( polyhedron.node( i ) );

        nodes.insert( nodes.end(), newNodes.begin(), newNodes.end() );

        auto intersection = Polyhedron< Coord > ( { faces, edges, nodes } );
        assert ( !std::isnan( intersection.volume() ) );
        return intersection;
      }


      /**
       * \ingroup geo3d
       * \brief implementation for an intersection between a polyhedron and a hyperplane
       *
       * \tparam  Coord  type of the global coordinate
       * \return the intersection line
       */
      template< class Coord >
      auto intersect ( const Polyhedron< Coord >& polyhedron, const HyperPlane< Coord >& plane ) -> Dune::VoF::Face< Coord >
      {
        if ( !plane )
          return Dune::VoF::Face< Coord >();

        using Coordinate = Coord;
        using E = typename std::array< Coord, 2 >;

        const double eps = 1e-12;


        std::vector< bool > isInner ( polyhedron.nodes().size(), false );
        std::vector< E > edges;

        std::size_t n = 0;
        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( plane.levelSet( polyhedron.node( i ) ) <= eps )
          {
            isInner[ i ] = true;
            n++;
          }

        // no or all nodes are inner
        if ( n == 0 || n == polyhedron.nodes().size() )
        {
          for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
            if ( plane.levelSet( polyhedron.node( i ) ) >= -eps )
              return Dune::VoF::Face< Coord > ( { polyhedron.node( i ) } );
          return Dune::VoF::Face< Coord >();
        }



        for ( const auto& face : polyhedron.faces() )
        {
          Coord lastIntersectionPoint ( 0.0 );

          for ( const auto& edge : face.edges() )
            if ( isInner[ edge.nodeId(0) ] ^ isInner[ edge.nodeId(1) ] )
            {
              const Coordinate isPoint = edge.intersection( plane );

              if ( lastIntersectionPoint != Coord( 0.0 ) )
                edges.push_back( E( { lastIntersectionPoint, isPoint } ) );
              else
                lastIntersectionPoint = isPoint;
            }
        }

        // Sort edges of intersection face
        for ( std::size_t i = 0; i+1 < edges.size(); ++i )
        {
          if ( ( edges[ i ][ 1 ] - edges[ i+1 ][ 0 ] ).two_norm() < eps )
            continue;

          const Coord point = edges[ i ][ 1 ];

          const auto pos = std::find_if( edges.begin() + i + 1, edges.end(), [ &point, eps ]( const E& other ) -> bool { return ( point - other[ 0 ] ).two_norm() < eps; } );
          if ( pos != edges.end() )
          {
            std::swap( edges[ i+1 ], *pos );
            continue;
          }

          const auto pos2 = std::find_if( edges.begin() + i + 1, edges.end(), [ &point, eps ]( const E& other ) -> bool { return ( point - other[ 1 ] ).two_norm() < eps; } );
          if ( pos2 != edges.end() )
          {
            std::swap( edges[ i+1 ], *pos2 );
            std::swap( edges[ i+1 ][ 0 ], edges[ i+1 ][ 1 ] );
            continue;
          }

        }

        std::vector< Coordinate > nodes;
        for ( std::size_t i = 0; i < edges.size(); ++i )
          nodes.push_back( edges[ i ][ 0 ] );

        assert( nodes.size() > 0 );

        return typename Dune::VoF::Face< Coordinate >( nodes );
      }

    } // namespace __impl

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_3D_INTERSECT_HH
