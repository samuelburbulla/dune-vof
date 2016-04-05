#ifndef DUNE_VOF_GEOMETRY_3D_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_3D_INTERSECT_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/2d/polygon.hh>

namespace Dune {

  namespace VoF {

    /**
     * \brief namespace containing specific implementations
     */
    namespace __impl {


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

        const double eps = 1e-8;

        std::vector< Coordinate > nodes;
        std::vector< Edge > edges;
        std::vector< Face > faces;

        std::vector< bool > isInner ( polyhedron.nodes().size(), false );
        std::vector< Coordinate > newNodes;
        Face intersectionFace;

        std::size_t n = 0;
        std::vector< std::size_t > p;

        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( halfSpace.levelSet( polyhedron.node( i ) ) >= eps )
          {
            isInner[ i ] = true;
            p.push_back( n );
            n++;
          }
          else
            p.push_back( -1 );


        for ( const auto& face : polyhedron.faces() )
        {
          Face newFace;
          std::size_t lastIntersectionPointId = -1;

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
              const Coordinate isNode = edge.intersection( halfSpace.boundary() );

              // Intersection point is corner 0
              if ( ( isNode - edge.node(0) ).two_norm() <= eps )
              {
                if ( lastIntersectionPointId != std::size_t(-1) )
                {
                  edges.push_back( {{ p[ edge.nodeId(0) ], lastIntersectionPointId }} );
                  newFace.emplace_back( edges.size() - 1 );

                  edges.push_back( {{ lastIntersectionPointId, p[ edge.nodeId(0) ] }} );
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIntersectionPointId = p[ edge.nodeId(0) ];
              }
              // Intersection point is corner 1
              else if ( ( isNode - edge.node(1) ).two_norm() <= eps )
              {
                if ( lastIntersectionPointId != std::size_t(-1) )
                {
                  edges.push_back( {{ lastIntersectionPointId, p[ edge.nodeId(1) ] }} );
                  newFace.emplace_back( edges.size() - 1 );

                  edges.push_back( {{ p[ edge.nodeId(1) ], lastIntersectionPointId }});
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIntersectionPointId = p[ edge.nodeId(1) ];
              }
              else
              {
              // Intersection point is new point
                std::size_t isNodeId;
                const auto pos = std::find_if( newNodes.begin(), newNodes.end(),
                  [ isNode, eps ]( const auto& other ) -> bool { return (isNode - other).two_norm() < eps; }
                );

                if ( pos == newNodes.end() )
                {
                  newNodes.push_back( isNode );
                  isNodeId = n + newNodes.size() - 1;
                }
                else
                  isNodeId = n + pos - newNodes.begin();


                if ( isInner[ edge.nodeId(0) ] )
                {
                  edges.push_back( {{ p[ edge.nodeId(0) ], isNodeId }} );
                  newFace.emplace_back( edges.size() - 1 );

                  if ( lastIntersectionPointId != std::size_t(-1) )
                  {
                    edges.push_back( {{ isNodeId, lastIntersectionPointId }} );
                    newFace.emplace_back( edges.size() - 1 );

                    edges.push_back( {{ lastIntersectionPointId, isNodeId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                }
                else
                {
                  if ( lastIntersectionPointId != std::size_t(-1) )
                  {
                    edges.push_back( {{ lastIntersectionPointId, isNodeId }} );
                    newFace.emplace_back( edges.size() - 1 );

                    edges.push_back( {{ isNodeId, lastIntersectionPointId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                  edges.push_back( {{ isNodeId, p[ edge.nodeId(1) ] }} );
                  newFace.emplace_back( edges.size() - 1 );
                }
                lastIntersectionPointId = isNodeId;
              }
            }
          }
          if ( newFace.size() > 2 )
            faces.push_back( newFace );
        }



        // sort subentities and nodes of intersection face
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
          std::swap( intersectionFace[ i+1 ], *pos );
        }

        if ( intersectionFace.size() > 2 )
          faces.push_back( intersectionFace );


        // Build new polytop
        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( isInner[ i ] )
            nodes.push_back( polyhedron.node( i ) );

        nodes.insert( nodes.end(), newNodes.begin(), newNodes.end() );


        return Polyhedron< Coord > ( { faces, edges, nodes } );
      }


      /**
       * \ingroup geo3d
       * \brief implementation for an intersection between a polyhedron and a hyperplane
       *
       * \tparam  Coord  type of the global coordinate
       * \return the intersection line
       */
      template< class Coord >
      auto intersect ( const Polyhedron< Coord >& polyhedron, const HyperPlane< Coord >& plane ) -> Polygon< Coord >
      {
        if ( !plane )
          return Polygon< Coord >();

        using Coordinate = Coord;
        using Face = typename std::vector< size_t >;
        using Edge = typename std::array< size_t, 2 >;

        const double eps = 1e-8;

        std::vector< Edge > edges;

        std::vector< bool > isInner ( polyhedron.nodes().size(), false );
        std::vector< Coord > newNodes;
        Face intersectionFace;

        std::size_t n = 0;
        std::vector< std::size_t > p;


        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( plane.levelSet( polyhedron.node( i ) ) >= eps )
          {
            isInner[ i ] = true;
            p.push_back( n );
            n++;
          }
          else
            p.push_back( -1 );


        for ( const auto& face : polyhedron.faces() )
        {
          std::size_t lastIntersectionPointId = -1;

          for ( const auto& edge : face.edges() )
          {
            // Edge is intersected by plane
            if ( isInner[ edge.nodeId(0) ] ^ isInner[ edge.nodeId(1) ] )
            {
              const Coordinate isNode = edge.intersection( plane );

              // Intersection point is corner 0
              if ( ( isNode - edge.node(0) ).two_norm() <= eps )
              {
                if ( lastIntersectionPointId != std::size_t(-1) )
                {
                  edges.push_back( {{ lastIntersectionPointId, p[ edge.nodeId(0) ] }} );
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIntersectionPointId = p[ edge.nodeId(0) ];
              }
              // Intersection point is corner 1
              else if ( ( isNode - edge.node(1) ).two_norm() <= eps )
              {
                if ( lastIntersectionPointId != std::size_t(-1) )
                {
                  edges.push_back( {{ p[ edge.nodeId(1) ], lastIntersectionPointId }});
                  intersectionFace.emplace_back( edges.size() - 1 );
                }
                lastIntersectionPointId = p[ edge.nodeId(1) ];
              }
              else
              {
              // Intersection point is new point
                std::size_t isNodeId;
                const auto pos = std::find_if( newNodes.begin(), newNodes.end(),
                  [ isNode, eps ]( const auto& other ) -> bool { return (isNode - other).two_norm() < eps; }
                );

                if ( pos == newNodes.end() )
                {
                  newNodes.push_back( isNode );
                  isNodeId = n + newNodes.size() - 1;
                }
                else
                  isNodeId = n + pos - newNodes.begin();


                if ( isInner[ edge.nodeId(0) ] )
                {
                  if ( lastIntersectionPointId != std::size_t(-1) )
                  {
                    edges.push_back( {{ lastIntersectionPointId, isNodeId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                }
                else
                {
                  if ( lastIntersectionPointId != std::size_t(-1) )
                  {
                    edges.push_back( {{ isNodeId, lastIntersectionPointId }} );
                    intersectionFace.emplace_back( edges.size() - 1 );
                  }
                }
                lastIntersectionPointId = isNodeId;
              }
            }
          }
        }

        // sort subentities and nodes of intersection face
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
          std::swap( intersectionFace[ i+1 ], *pos );
        }

        std::vector< Coordinate > nodes;
        for ( std::size_t i = 0; i < polyhedron.nodes().size(); ++i )
          if ( isInner[ i ] )
            nodes.push_back( polyhedron.node( i ) );

        nodes.insert( nodes.end(), newNodes.begin(), newNodes.end() );


        if ( intersectionFace.size() < 2 )
          return Line< Coord >();
        else
        {
          std::vector< Coordinate > polygonNodes;
          for ( std::size_t i = 0; i < intersectionFace.size(); ++i )
            polygonNodes.push_back( nodes[ edges[ intersectionFace[ i ] ][ 0 ] ] );

          return make_polygon( polygonNodes );
        }
      }

    } // namespace __impl

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_3D_INTERSECT_HH
