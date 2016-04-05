#ifndef DUNE_VOF_GEOMETRY_POLYHEDRON_HH
#define DUNE_VOF_GEOMETRY_POLYHEDRON_HH

/* c++ includes */
#include <array>
#include <tuple>
#include <vector>


namespace Dune {

  namespace VoF {


    template < class Polyhedron >
    struct Edge
    {
      using Coordinate = typename Polyhedron::Coordinate;

      Edge ( const Polyhedron *parent, const std::array< std::size_t, 2 >& nodeIds )
       : parent_ ( parent ), nodeIds_ ( nodeIds )
      {
        assert ( nodeIds.size() == 2 );
      };

      Edge ( const Edge& other, const Polyhedron *parent )
       : parent_ ( parent ), nodeIds_ ( other.nodeIds() ) {}


      const bool operator== ( const Edge& other ) const
      {
        return nodeIds_ == other.nodeIds(); // && &parent_ == &other.parent()
      }

      const Polyhedron *parent () const { return parent_; }

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

        Coordinate normal;
        normal[ 0 ] = v1[ 1 ] * v2[ 2 ] - v1[ 2 ] * v2[ 1 ];
        normal[ 1 ] = v1[ 2 ] * v2[ 0 ] - v1[ 0 ] * v2[ 2 ];
        normal[ 2 ] = v1[ 0 ] * v2[ 1 ] - v1[ 1 ] * v2[ 0 ];

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
        assert( ( hyperplane.normal() * node(0) - hyperplane.distance() <= 0 ) ^  ( hyperplane.normal() * node(1) - hyperplane.distance() <= 0 ) );

        Coordinate p = node(0);
        p -= node(1);
        p *= ( ( hyperplane.normal() * node(0) ) - hyperplane.distance() ) / ( hyperplane.normal() * p * -1.0 );
        p += node(0);
        return p;
      }

    private:
      const Polyhedron *parent_;
      std::array< std::size_t, 2 > nodeIds_;
    };




    template < typename Polyhedron >
    struct Face
    {
      using Coordinate = typename Polyhedron::Coordinate;
      using Edge = typename Polyhedron::E;

      Face ( const Polyhedron *parent, const std::vector< Edge >& edges )
       : parent_ ( parent ), edges_( edges )
      {
        for ( const auto& edge : edges )
          nodeIds_.push_back( edge.nodeId( 0 ) );
      };

      Face ( const Face& other, const Polyhedron *parent )
       : parent_ ( parent ), nodeIds_ ( other.nodeIds() )
      {
        for ( const auto& edge : other.edges() )
          edges_.emplace_back( Edge ( edge, parent ) );
      }


      const Polyhedron *parent () const { return parent_; }

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
      const Polyhedron *parent_;
      std::vector< std::size_t > nodeIds_;
      std::vector< Edge > edges_;
    };




    template < class Coord >
    struct Polyhedron
    {
      using Coordinate = Coord;
      using F = Face< Polyhedron >;
      using E = Edge< Polyhedron >;
      static_assert(Coord::dimension == 3, "Dimension must be == 3." );

      Polyhedron (
        const std::vector< std::vector< std::size_t > >& faces,
        const std::vector< std::array< size_t, 2 > >& edges,
        const std::vector< Coordinate >& nodes
      ) : nodes_ ( nodes )
      {

        for ( const auto& edgeData : edges )
          edges_.emplace_back( E ( this, edgeData ) );

        for ( const auto& faceData : faces )
        {
          std::vector< E > faceEdges;

          for ( const std::size_t id : faceData )
            faceEdges.push_back( edges_[ id ] );

          faces_.emplace_back( F ( this, faceEdges ) );
        }

      };

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

      const std::vector< Coordinate >& nodes () const { return nodes_; }

      const Coordinate& node ( const std::size_t index ) const { return nodes_[ index ]; }

      const F& face ( const std::size_t index ) const { return faces_[ index ]; }

      const std::vector< F >& faces () const { return faces_; }

      const std::vector< E >& edges () const { return edges_; }


      const double volume () const
      {
        double vol = 0.0;
        for ( const auto& face : faces() )
          vol += face.center() * face.outerNormal() * face.volume();
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

        //Dune::Exception( "Warning: (attachedFaceToEdge) Edge was not found in any face." );
        std::cerr << "Edge was not found in any face." << std::endl;
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




    template< class Geometry >
    static inline auto make_polygon( const Geometry& geometry ) -> Polyhedron< typename Geometry::GlobalCoordinate >
    {
      using Container = std::vector< typename Geometry::GlobalCoordinate >;
      auto type = geometry.type();

      if ( type.isSimplex() )
      {
        const std::vector< DomainVector > nodes { geometry.corner( 0 ), geometry.corner( 1 ), geometry.corner( 2 ), geometry.corner( 3 ) };

        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 1 }}, {{ 1, 0 }}, {{ 0, 2 }}, {{ 2, 0 }}, {{ 1, 2 }}, {{ 2, 1 }}, {{ 0, 3 }}, {{ 3, 0 }}, {{ 1, 3 }}, {{ 3, 1 }}, {{ 3, 2 }}, {{ 2, 3 }}
        };

        const std::vector< std::vector< std::size_t > > faces { {{ 2, 5, 1 }}, {{ 0, 8, 7 }}, {{ 6, 10, 3 }}, {{ 4, 11, 9 }} };

        return Dune::VoF::Polyhedron< DomainVector > ( faces, edges, nodes );
      }
      else if ( type.isQube() )
      {
        const std::vector< DomainVector > nodes;
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
        return Dune::VoF::Polyhedron< DomainVector > ( faces, edges, nodes );
      }
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }



  } // end of namespace VoF

} // end of namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYHEDRON_HH
