#ifndef DUNE_VOF_RECONSTRUCTION_GENERALHEIGHTFUNCTION_HH
#define DUNE_VOF_RECONSTRUCTION_GENERALHEIGHTFUNCTION_HH

#include <cmath>

#include <algorithm>
#include <limits>
#include <vector>
#include <tuple>

#include <dune/common/fmatrix.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/polytope.hh>
#include <dune/vof/geometry/2d/intersect.hh>
#include <dune/vof/geometry/utility.hh>


namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Reconstruction
     * \brief   general height function reconstruction operator
     *
     * \tparam DF   discrete function type
     * \tparam RS   reconstruction set type
     * \tparam IR   initial reconstruction type
     */
    template< class GV, class VSt, class IR >
    struct GeneralHeightFunctionReconstruction
    {
      using GridView = GV;
      using VertexStencilSet = VSt;
      using InitialReconstruction = IR;

    private:
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using Reconstruction = HalfSpace< Coordinate >;

      static constexpr int dim = Coordinate::dimension;
      static constexpr int noc = static_cast< int > ( std::pow( 3, dim-1 ) );
      static constexpr double length = 3;
      static constexpr double width = 1;

      using Line = Dune::VoF::Line< Coordinate >;
      using Lines = std::array< std::array< Line, 2 >, noc >;
      using Segments = std::array< std::array< std::vector< std::array< double, 3 > >, 2 >, noc >;
      using Heights = Dune::FieldVector< double, noc >;

    public:
      explicit GeneralHeightFunctionReconstruction ( const VertexStencilSet &vertexStencilSet ) : initializer_( vertexStencilSet ) {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  ColorFunction
       * \tparam  ReconstructionSet
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class ColorFunction, class ReconstructionSet, class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags, bool communicate = false )
      {
        initializer_( color, reconstructions, flags );

        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, color, flags, reconstructions );
        }

        if ( communicate )
          reconstructions.communicate();
      }

      /**
       * \brief   (local) operator application
       *
       * \tparam  ColorFunction
       * \tparam  Flags
       * \tparam  ReconstructionSet
       * \param   entity          current element
       * \param   color           color functions
       * \param   flags           set of flags
       * \param   reconstructions reconstruction set
       */
      template< class ColorFunction, class Flags, class ReconstructionSet >
      void applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, ReconstructionSet &reconstructions )
      {
        const Coordinate &normal = reconstructions[ entity ].innerNormal();
        const Coordinate centroid = interface( entity, reconstructions ).centroid();

        const double h = std::pow( entity.geometry().volume(), 1.0 / dim );

        Coordinate shift = generalizedCrossProduct( normal );

        Lines lines;
        for( int i = 0; i < 3; ++i )
        {
          Coordinate mid = centroid;
          mid.axpy( ( i-1 ) * width * h, shift );

          Coordinate up = mid;
          Coordinate down = mid;

          up.axpy(  -length * h, normal );
          down.axpy( length * h, normal );

          lines[ i ][ 0 ] = Line( mid, up );
          lines[ i ][ 1 ] = Line( mid, down );
        }

        // Get every height function value as tuple of double: length, colorvalue, levelset of centroid
        Segments segments = computeSegments( color, lines, centroid, h, reconstructions[ entity ] );

        Heights heights ( 0.0 );
        double TOL = std::numeric_limits< double >::epsilon();

        for( std::size_t i = 0; i < noc; ++i )
          for( std::size_t r = 0; r < 2; ++r )
          {
            std::vector< std::array< double, 3 > > &segs = segments[ i ][ r ];
            double lastU = segs[ 0 ][ 1 ];
            heights[ i ] += lastU * segs[ 0 ][ 0 ];

            for ( std::size_t s = 1; s < segs.size(); ++s )
            {
              double u = segs[ s ][ 1 ];

              if ( r == 0 && u > lastU - TOL )
                u = 0.0;

              if ( r == 1 && u < lastU + TOL )
                u = 1.0;

              heights[ i ] += u * segs[ s ][ 0 ];
              lastU = u;
            }
          }

        Coordinate newNormal = computeNormal( heights, h, normal );

        reconstructions[ entity ] = locateHalfSpace( makePolytope( entity.geometry() ), newNormal, color[ entity ] );
      }

    private:

      template< class ColorFunction >
      Segments computeSegments( const ColorFunction &color, const Lines &lines, const Coordinate& center, const double h, const Reconstruction &reconstruction ) const
      {
        Segments segments;

        for ( const auto &entity : elements( color.gridView() ) )   // TODO: Do not run over the whole grid!
        {
          if ( ( entity.geometry().center() - center ).two_norm() > std::sqrt( ( length * h ) * ( length * h ) + h * h ) + h )
            continue;

         for ( int i = 0; i < noc; ++i )
          for ( int r = 0; r < 2; ++r )
          {
            Line segment = getSegment( lines[ i ][ r ], entity, color.gridView() );

            double segmentVolume = segment.volume();
            if ( segmentVolume < std::numeric_limits< double >::epsilon() )
              continue;

            segments[ i ][ r ].push_back( {{ segmentVolume, color[ entity ], reconstruction.levelSet( segment.centroid() ) }} );
          }
        }

        auto compDown = []( const std::array< double, 3 > &a, const std::array< double, 3 > &b ) -> bool { return a[ 2 ] < b[ 2 ]; };
        for ( int i = 0; i < noc; ++i )
          std::sort( segments[ i ][ 1 ].begin(), segments[ i ][ 1 ].end(), compDown );

        auto compUp = []( const std::array< double, 3 > &a, const std::array< double, 3 > &b ) -> bool { return a[ 2 ] > b[ 2 ]; };
        for ( int i = 0; i < noc; ++i )
          std::sort( segments[ i ][ 0 ].begin(), segments[ i ][ 0 ].end(), compUp );

        return segments;
      };

    #if GRIDDIM == 2
      Coordinate computeNormal ( const Heights &heights, const double h, const Coordinate &normal ) const
      {
        double Hx;
        if ( heights[ 0 ] == 0.0 )
          Hx = ( heights[ 2 ] - heights[ 1 ] ) / ( width * h );
        else if ( heights[ 2 ] == 0.0 )
          Hx = ( heights[ 1 ] - heights[ 0 ] ) / ( width * h );
        else
          Hx = ( heights[ 2 ] - heights[ 0 ] ) / ( width * 2.0 * h );

        Coordinate n ( { Hx, -1.0 } );
        normalize( n );

        // Rotate to reference frame
        Coordinate newNormal ( 0 );
        newNormal[ 0 ] = - n[ 0 ] * normal[ 1 ] - n[ 1 ] * normal[ 0 ];
        newNormal[ 1 ] = n[ 0 ] * normal[ 0 ] - n[ 1 ] * normal[ 1 ];

        normalize( newNormal );
        return newNormal;
      }

    #elif GRIDDIM == 3
      Coordinate computeNormal ( const Heights &heights, const Coordinate &normal ) const
      {
        double Hx;
        if ( heights[ 5 ] == 0.0 )
          Hx = ( heights[ 4 ] - heights[ 3 ] );
        else if ( heights[ 3 ] == 0.0 )
          Hx = ( heights[ 5 ] - heights[ 4 ] );
        else
          Hx = ( heights[ 5 ] - heights[ 3 ] ) / 2.0;

        double Hy;
        if ( heights[ 7 ] == 0.0 )
          Hy = ( heights[ 4 ] - heights[ 1 ] );
        else if ( heights[ 1 ] == 0.0 )
          Hy = ( heights[ 7 ] - heights[ 4 ] );
        else
          Hy = ( heights[ 7 ] - heights[ 1 ] ) / 2.0;


        int i = std::get< 0 >( orientation );
        int j = std::get< 1 >( orientation );

        Coordinate newNormal;
        newNormal[ i ] = -j;
        newNormal[ (i+1)%3 ] = Hx * -j;
        newNormal[ (i+2)%3 ] = Hy;

        normalize( newNormal );
        return newNormal;
      }
    #endif

      template< class Entity, class GridView >
      Line getSegment ( const Line& line, const Entity& entity, const GridView &gridView ) const
      {
        Line segment( line );

        for ( const auto &intersection : intersections( gridView, entity ) )
        {
          auto normal = intersection.centerUnitOuterNormal();
          normal *= -1.0;
          const auto center = intersection.geometry().center();

          const HalfSpace< Coordinate > halfspace( normal, center );

          segment = intersect( segment, halfspace );

          if ( segment == Line() )
            return Line();
        }
        return segment;
      }

      InitialReconstruction initializer_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif
