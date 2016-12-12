#ifndef DUNE_VOF_GEOMETRY_POLYHEDRON_HH
#define DUNE_VOF_GEOMETRY_POLYHEDRON_HH

/* c++ includes */
#include <array>
#include <tuple>
#include <vector>
#include <cmath>

/* dune includes */
#include <dune/vof/geometry/utility.hh>

namespace Dune {

  namespace VoF {

    namespace __impl {

      template < class ParentType >
      struct Edge
      {
        using Coordinate = typename ParentType::Coordinate;
        template< class C > friend class Polyhedron;
        template< class P > friend class Face;

        Edge ( ParentType const *parent, const std::array< std::size_t, 2 >& nodeIds )
         : parent_ ( parent ), nodeIds_ ( nodeIds )
        {
          assert ( nodeIds.size() == 2 );
        };

        Edge ( const Edge& other, ParentType const *parent )
         : parent_ ( parent ), nodeIds_ ( other.nodeIds() ) {}


        const bool operator== ( const Edge& other ) const
        {
          return nodeIds_ == other.nodeIds(); // && &parent_ == &other.parent()
        }

        const ParentType *parent () const { return parent_; }
      private:
        void rebind( ParentType const *parent ) { parent_ = parent; }
      public:
        const std::array< std::size_t, 2 >& nodeIds() const { return nodeIds_; }

        const std::size_t& nodeId ( const std::size_t index ) const { return nodeIds_[ index ]; }

        const Coordinate& node ( const std::size_t index ) const { return parent_->node( nodeIds_[ index ] ); }

        const Coordinate outerNormal () const
        {
          Coordinate diff = node( 1 );
          diff -= node( 0 );

          Coordinate normal;
          normal[ 0 ] = diff[ 1 ];
          normal[ 1 ] = - diff[ 0 ];

          normal /= normal.two_norm();
          return normal;
        }

        const Coordinate outerNormal ( const Coordinate& planeNormal ) const
        {
          const Coordinate v1 = node( 1 ) - node( 0 );
          const Coordinate& v2 = planeNormal;

          //normal = generalizedCrossProduct( node( 1 ) - node( 0 ), planeNormal );

          Coordinate normal;
          normal[ 0 ] = v1[ 1 ] * v2[ 2 ] - v1[ 2 ] * v2[ 1 ];
          normal[ 1 ] = v1[ 2 ] * v2[ 0 ] - v1[ 0 ] * v2[ 2 ];
          normal[ 2 ] = v1[ 0 ] * v2[ 1 ] - v1[ 1 ] * v2[ 0 ];

          assert ( normal.two_norm() != 0 );
          normal /= normal.two_norm();
          return normal;
        }

        const Coordinate center () const
        {
          Coordinate center = node( 0 );
          center += node( 1 );
          center *= 0.5;
          return center;
        }

        const double volume () const
        {
          return ( node( 1 ) - node( 0 ) ).two_norm();
        }

        template < class Hyperplane >
        const Coordinate intersection ( const Hyperplane& hyperplane ) const
        {
          Coordinate p = node(0);
          p -= node(1);
          p *= ( ( hyperplane.normal() * node(0) ) + hyperplane.distance() ) / ( hyperplane.normal() * p * -1.0 );
          p += node(0);
          return p;
        }

      private:
        ParentType const *parent_;
        std::array< std::size_t, 2 > nodeIds_;
      };




      template < typename ParentType >
      struct Face
      {
        using Coordinate = typename ParentType::Coordinate;
        using Edge = typename ParentType::E;
        template< class C > friend class Polyhedron;

        Face ( ParentType const *parent, const std::vector< Edge >& edges )
         : parent_ ( parent ), edges_( edges )
        {
          for ( const auto& edge : edges )
            nodeIds_.push_back( edge.nodeId( 0 ) );
        };

        Face ( const Face& other, ParentType const *parent )
         : parent_ ( parent ), nodeIds_ ( other.nodeIds() )
        {
          for ( const auto& edge : other.edges() )
            edges_.emplace_back( Edge ( edge, parent ) );
        }


        const ParentType *parent () const { return parent_; }
      private:
        void rebind( ParentType const *parent )
        {
          parent_ = parent;
          for ( auto& edge : edges_ )
            edge.rebind( parent );
        }
      public:
        const std::vector< Edge >& edges () const { return edges_; }

        const Edge& edge ( const std::size_t index ) const { return edges_[ index ]; }

        const std::vector< std::size_t >& nodeIds () const { return nodeIds_; }

        const std::size_t& nodeId ( const std::size_t index ) const { return nodeIds_[ index ]; }

        const Coordinate& node ( const std::size_t index ) const { return parent_->node( nodeIds_[ index ] ); }


        const bool operator== ( const Face& other ) const
        {
          if ( &parent() != &other.parent() )
            return false;

          for ( const auto id : nodeIds() )
          {
            bool found = false;
            for ( const auto id2 : other.nodeIds() )
              if( id == id2 ) found = true;

            if ( !found )
              return false;
          }

          for ( const auto& edge : edges() )
          {
            bool found = false;
            for ( const auto& edge2 : other.edges() )
              if( edge == edge2 )
                found = true;

            if ( !found )
              return false;
          }

          return true;
        }

        const Coordinate outerNormal () const
        {
          assert ( nodeIds().size() >= 3 );

          //normal = generalizedCrossProduct( node( 0 ) - node( 1 ), node( 2 ) - node( 1 ) );

          const auto& v1 = node( 0 ) - node( 1 );
          const auto& v2 = node( 2 ) - node( 1 );

          Coordinate normal;
          normal[ 0 ] = v1[ 1 ] * v2[ 2 ] - v1[ 2 ] * v2[ 1 ];
          normal[ 1 ] = v1[ 2 ] * v2[ 0 ] - v1[ 0 ] * v2[ 2 ];
          normal[ 2 ] = v1[ 0 ] * v2[ 1 ] - v1[ 1 ] * v2[ 0 ];

          normal /= -1.0 * normal.two_norm();

          return normal;
        }

        const Coordinate center () const
        {
          Coordinate center ( 0.0 );
          for ( std::size_t i = 0; i < nodeIds_.size(); ++i )
            center += node( i );
          center *= 1.0 / nodeIds_.size();
          return center;
        }

        const double volume () const
        {
          double vol = 0;
          for ( const auto& edge : edges() )
            vol += edge.center() * edge.outerNormal( outerNormal() ) * edge.volume();
          return std::abs( vol / 2.0 );
        }

      private:
        ParentType const *parent_;
        std::vector< std::size_t > nodeIds_;
        std::vector< Edge > edges_;
      };

    } // namespace __impl



    template < class Coord >
    struct Polyhedron
    {
      using Coordinate = Coord;
      using F = typename __impl::Face< Polyhedron >;
      using E = typename __impl::Edge< Polyhedron >;

      static_assert( Coord::dimension == 3, "Dimension must be == 3." );

      Polyhedron () = default;

      Polyhedron (
        const std::vector< std::vector< std::size_t > >& faces,
        const std::vector< std::array< size_t, 2 > >& edges,
        const std::vector< Coordinate >& nodes
      ) : nodes_ ( std::move( nodes ) )
      {
        for ( const auto& edgeData : edges )
          edges_.emplace_back( this, edgeData );

        for ( const auto& faceData : faces )
        {
          std::vector< E > faceEdges;

          for ( const std::size_t id : faceData )
            faceEdges.push_back( edges_[ id ] );

          faces_.emplace_back( this, faceEdges );
        }
      };

    private:
      Polyhedron ( std::vector< F > faces, std::vector< E > edges, std::vector< Coordinate > nodes )
        : nodes_ ( std::move( nodes ) ), faces_ ( std::move( faces ) ), edges_ ( std::move( edges ) )
      {
        for ( auto& edge : edges_ )
          edge.rebind( this );

        for ( auto& face : faces_ )
          face.rebind( this );
      };

    public:
      Polyhedron ( const Polyhedron& other )
        : Polyhedron( other.faces_, other.edges_, other.nodes_ )
      {}

      Polyhedron ( Polyhedron&& other )
        : Polyhedron( std::move( other.faces_ ), std::move( other.edges_ ), std::move( other.nodes_ ) )
      {}

      Polyhedron& operator= ( Polyhedron&& ) = delete;
      Polyhedron& operator= ( const Polyhedron& ) = delete;


      Polyhedron ( const Polyhedron& other, const std::vector< Coord >& newNodes )
       : nodes_ ( newNodes )
      {
        for ( const F& face : other.faces() )
          faces_.emplace_back( F ( face, this ) );
        for ( const E& edge : other.edges() )
          edges_.emplace_back( E ( edge, this ) );
      }

      bool operator== ( const Polyhedron& other ) const
      {
        return nodes_ == other.nodes() && faces_ == other.faces() && edges_ == other.edges();
      }

      std::size_t size() const { return nodes().size(); }

      const std::vector< Coordinate >& nodes () const { return nodes_; }

      const Coordinate& node ( const std::size_t index ) const { return nodes_[ index ]; }
      const Coordinate& vertex ( const std::size_t index ) const { return node( index ); }

      const F& face ( const std::size_t index ) const { return faces_[ index ]; }

      const std::vector< F >& faces () const { return faces_; }

      const std::vector< E >& edges () const { return edges_; }


      const double volume () const
      {
        double vol = 0.0;
        for ( const auto& face : faces() )
          vol += face.center() * face.outerNormal() * face.volume();
        assert( !std::isnan( vol ) );
        return std::abs( vol / 3.0 );
      }

      const Coordinate center () const
      {
        Coordinate center ( 0.0 );
        for ( const auto& node : nodes_ )
          center += node;
        center *= 1.0 / nodes_.size();
        return center;
      }

      const F& attachedFaceToEdge ( const E& e ) const
      {
        for ( const F& face : faces() )
          for ( const E& edge : face.edges() )
            if ( edge == e )
              return face;

        DUNE_THROW( InvalidStateException, "Edge was not found in any face." );
        return face( 0 );
      }

      std::vector< E > attachedEdgesToNode ( const std::size_t index ) const
      {
        std::vector< E > attachedEdges;
        for ( const auto& edge : edges() )
        {
          if ( edge.nodeId( 0 ) == index )
            attachedEdges.push_back( edge );
        }

        return attachedEdges;
      }

    private:
      std::vector< Coordinate > nodes_;
      std::vector< F > faces_;
      std::vector< E > edges_;
    };



    // makePolyhedron
    // --------------
    template< class Geometry >
    static inline auto makePolyhedron ( const Geometry& geometry ) -> Polyhedron< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if ( type.isSimplex() )
      {
        const Container nodes { geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 2 ), geometry.corner( 3 ) };

        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 1 }}, {{ 1, 0 }}, {{ 0, 2 }}, {{ 2, 0 }}, {{ 1, 2 }}, {{ 2, 1 }}, {{ 0, 3 }}, {{ 3, 0 }}, {{ 1, 3 }}, {{ 3, 1 }}, {{ 3, 2 }}, {{ 2, 3 }}
        };

        const std::vector< std::vector< std::size_t > > faces { {{ 2, 5, 1 }}, {{ 0, 8, 7 }}, {{ 6, 10, 3 }}, {{ 4, 11, 9 }} };

        return Polyhedron< typename Geometry::GlobalCoordinate > ( faces, edges, nodes );
      }
      else if ( type.isCube() )
      {
        Container nodes;
        for ( std::size_t i = 0; i < 8; ++i )
         nodes.push_back( geometry.corner( i ) );

        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 4 }}, {{ 1, 5 }}, {{ 2, 6 }}, {{ 3, 7 }}, {{ 0, 2 }}, {{ 1, 3 }}, {{ 0, 1 }}, {{ 2, 3 }}, {{ 4, 6 }}, {{ 5, 7 }},
          {{ 4, 5 }}, {{ 6, 7 }}, {{ 4, 0 }}, {{ 5, 1 }}, {{ 6, 2 }}, {{ 7, 3 }}, {{ 2, 0 }}, {{ 3, 1 }}, {{ 1, 0 }}, {{ 3, 2 }},
          {{ 6, 4 }}, {{ 7, 5 }}, {{ 5, 4 }}, {{ 7, 6 }}
        };

        const std::vector< std::vector< std::size_t > > faces {
          {{ 0, 8, 14, 16 }}, {{ 5, 3, 21, 13 }}, {{ 6, 1, 22, 12 }}, {{ 19, 2, 11, 15 }}, {{ 4, 7, 17, 18 }}, {{ 10, 9, 23, 20 }}
        };
        return Polyhedron< typename Geometry::GlobalCoordinate >( faces, edges, nodes );
      }
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }

    template< class Geometry, class Map >
    static inline auto makePolyhedron ( const Geometry& geometry, Map&& map )
      -> Polyhedron< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if ( type.isSimplex() )
      {
        const Container nodes { map( geometry.corner( 0 ) ), map( geometry.corner( 1 ) ), map( geometry.corner( 2 ) ), map( geometry.corner( 3 ) ) };

        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 1 }}, {{ 1, 0 }}, {{ 0, 2 }}, {{ 2, 0 }}, {{ 1, 2 }}, {{ 2, 1 }}, {{ 0, 3 }}, {{ 3, 0 }}, {{ 1, 3 }}, {{ 3, 1 }}, {{ 3, 2 }}, {{ 2, 3 }}
        };

        const std::vector< std::vector< std::size_t > > faces { {{ 2, 5, 1 }}, {{ 0, 8, 7 }}, {{ 6, 10, 3 }}, {{ 4, 11, 9 }} };

        return Polyhedron< typename Geometry::GlobalCoordinate > ( faces, edges, nodes );
      }
      else if ( type.isCube() )
      {
        Container nodes;
        for ( std::size_t i = 0; i < 8; ++i )
         nodes.push_back( map( geometry.corner( i ) ) );

        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 4 }}, {{ 1, 5 }}, {{ 2, 6 }}, {{ 3, 7 }}, {{ 0, 2 }}, {{ 1, 3 }}, {{ 0, 1 }}, {{ 2, 3 }}, {{ 4, 6 }}, {{ 5, 7 }},
          {{ 4, 5 }}, {{ 6, 7 }}, {{ 4, 0 }}, {{ 5, 1 }}, {{ 6, 2 }}, {{ 7, 3 }}, {{ 2, 0 }}, {{ 3, 1 }}, {{ 1, 0 }}, {{ 3, 2 }},
          {{ 6, 4 }}, {{ 7, 5 }}, {{ 5, 4 }}, {{ 7, 6 }}
        };

        const std::vector< std::vector< std::size_t > > faces {
          {{ 0, 8, 14, 16 }}, {{ 5, 3, 21, 13 }}, {{ 6, 1, 22, 12 }}, {{ 19, 2, 11, 15 }}, {{ 4, 7, 17, 18 }}, {{ 10, 9, 23, 20 }}
        };
        return Polyhedron< typename Geometry::GlobalCoordinate >( faces, edges, nodes );
      }
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }



  } // end of namespace VoF

} // end of namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYHEDRON_HH
