#ifndef DUNE_VOF_SWARTZCURVATURE_HH
#define DUNE_VOF_SWARTZCURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-grid includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Method
     * \brief curvature
     *
     * \tparam  GV  grid view
     * \tparam  ST  stencils
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class ST >
    struct SwartzCurvature
    {
      using GridView = GV;
      using StencilSet = ST;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit SwartzCurvature ( const StencilSet &stencils )
       : gridView_( stencils.gridView() ), stencils_( stencils )
      {}

      template< class DF, class ReconstructionSet, class Flags, class CurvatureSet >
      void operator() ( const DF &color, const ReconstructionSet &reconstructions, const Flags &flags,
                         CurvatureSet &curvatureSet, bool communicate = false ) const
      {
        for ( const auto& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          curvatureSet[ entity ] = 0.0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, reconstructions, flags, curvatureSet );
        }

        curvatureSet.communicate();

        CurvatureSet newCurvature ( curvatureSet );
        for ( const auto& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          if ( curvatureSet[ entity ] == 0.0 )
            averageCurvature( entity, curvatureSet, flags, newCurvature );

        }
        curvatureSet = newCurvature;

        if ( communicate )
          curvatureSet.communicate();
      }

    private:
      template< class ReconstructionSet, class Flags, class CurvatureSet >
      void applyLocal ( const Entity &entity, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvatureSet ) const
      {
        double weights = 0.0;
        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate centerEn = entity.geometry().center();

        //second moment: mu = |I|^2/12
        double muEn = interfaceEn.volume();
        muEn *= muEn / 12.0;

        for( const auto& neighbor1 : stencil( entity ) )
        {
          if ( !flags.isMixed( neighbor1 ) )
            continue;

          auto interfaceNb1 = interface( neighbor1, reconstructions );
          Coordinate centroidNb1 = interfaceNb1.centroid();
          Coordinate centerNb1 = neighbor1.geometry().center();

          Coordinate midwayNb1 = centroidNb1 + centroidEn;
          midwayNb1 *= 0.5;

          double triangleNb1 = ( centroidNb1 - centroidEn ).two_norm();

          double muNb1 = interfaceNb1.volume();
          muNb1 *= muNb1 / 12.0;

          Coordinate tanNb1 = centroidEn - centroidNb1;
          tanNb1 /= tanNb1.two_norm();

          Coordinate nu12 = generalizedCrossProduct( tanNb1 );
          if ( nu12 * reconstructions[ entity ].innerNormal() < 0 )
            continue;

          for( const auto& neighbor2 : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor2 ) )
              continue;

            if ( neighbor1 == neighbor2 )
              continue;

            auto interfaceNb2 = interface( neighbor2, reconstructions );
            Coordinate centroidNb2 = interfaceNb2.centroid();
            Coordinate centerNb2 = neighbor2.geometry().center();

            Coordinate midwayNb2 = centroidNb2 + centroidEn;
            midwayNb2 *= 0.5;

            double triangleNb2 = ( centroidNb2 - centroidEn ).two_norm();

            double muNb2 = interfaceNb2.volume();
            muNb2 *= muNb2 / 12.0;

            Coordinate midwayDiff = midwayNb2 - midwayNb1;
            double deltaS = midwayDiff.two_norm();

            Coordinate tanNb2 = centroidNb2 - centroidEn;
            tanNb2 /= tanNb2.two_norm();

            Coordinate nu23 = generalizedCrossProduct( tanNb2 );
            if ( nu23 * reconstructions[ entity ].innerNormal() < 0 )
              continue;

            if ( midwayDiff * ( tanNb1 + tanNb2 ) < 0 )
              continue;

            static_assert( GridView::dimension == 2, "Only implemented in 2D yet!" );
            auto times = []( const Coordinate &a, const Coordinate &b ) -> double { return a[0] * b[1] - a[1] * b[0]; };

            // taking all together
            double M = ( ( muNb2 - muEn ) / triangleNb2 - ( muEn - muNb1 ) / triangleNb1 ) / ( 2.0 * deltaS );

            double weight = 1.0; //interfaceNb1.volume() * interfaceNb2.volume();
            curvatureSet[ entity ] += weight * times( nu12, nu23 ) / ( ( 1.0 + M ) * deltaS );
            weights += weight;
          }
        }

        if ( weights > 0 )
          curvatureSet[ entity ] /= weights;
      }

      template< class CurvatureSet, class Flags >
      void averageCurvature( const Entity &entity, const CurvatureSet &curvatureSet, const Flags &flags,
                             CurvatureSet &newCurvature ) const
      {
        int n = 0;

        for ( const auto& neighbor : stencil( entity ) )
        {
          if ( !flags.isMixed( neighbor ) )
            continue;

          newCurvature[ entity ] += curvatureSet[ neighbor ];
          n++;
        }

        if ( n > 0 )
          newCurvature[ entity ] /= static_cast< double >( n );
      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      GridView gridView_;
      const StencilSet &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_SWARTZCURVATURE_HH
