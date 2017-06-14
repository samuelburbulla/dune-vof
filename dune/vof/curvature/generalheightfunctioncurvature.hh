#ifndef DUNE_VOF_GENERALHEIGHTFUNCTIONCURVATURE_HH
#define DUNE_VOF_GENERALHEIGHTFUNCTIONCURVATURE_HH

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

    // Calculate the curvature of the inferface in each mixed cell by the generalized height function method.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV  grid view
     * \tparam  ST  stencils
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class ST, class DF, class RS, class FL >
    struct GeneralHeightFunctionCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit GeneralHeightFunctionCurvature ( GridView gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils )
      {}

      template< class CurvatureSet >
      void operator() ( const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvatureSet ) const
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvatureSet[ entity ] = 0.0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, reconstructions, flags, curvatureSet );
        }

        curvatureSet.communicate();
      }

    private:
      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const ReconstructionSet &reconstructions, const Flags &flags,
                        CurvatureSet &curvatureSet ) const
      {
        double dx2uSum = 0.0;
        double dx2uw = 0.0;
        double dxuSum = 0.0;
        double dxuw = 0.0;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();
        Coordinate normalOrth = generalizedCrossProduct( normalEn );

        double xe = normalOrth * centroidEn;
        double ue = normalEn * centroidEn;

        for( const auto& neighbor1 : stencil( entity ) )
        {
          if ( !flags.isMixed( neighbor1 ) )
            continue;

          auto interfaceNb1 = interface( neighbor1, reconstructions );
          Coordinate centroidNb1 = interfaceNb1.centroid();

          double xn1 = normalOrth * centroidNb1;
          double un1 = normalEn * centroidNb1;

          double dxu = ( un1 - ue ) / ( xn1 - xe );
          const double weight = interfaceNb1.volume();
          dxuSum += dxu * weight;
          dxuw += weight;

          for( const auto& neighbor2 : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor2 ) )
              continue;

            if ( neighbor1 == neighbor2 )
              continue;

            auto interfaceNb2 = interface( neighbor2, reconstructions );
            Coordinate centroidNb2 = interfaceNb2.centroid();

            double xn2 = normalOrth * centroidNb2;
            double un2 = normalEn * centroidNb2;

            double tmp = 0.5 * ( xn2 - xn1 ) * ( xn2 - xe ) * ( xn1 - xe );
            const double dx2u = ( ( un2 - ue ) * ( xn1 - xe ) - ( un1 - ue ) * ( xn2 - xe ) ) / tmp;
            const double weight = interfaceNb2.volume() * interfaceNb1.volume();
            dx2uSum += dx2u * weight;
            dx2uw += weight;
          }
        }

        double dxu = dxuSum / dxuw;
        double dx2u = dx2uSum / dx2uw;

        curvatureSet[ entity ] = dx2u / std::pow( 1.0 + dxu * dxu, 3.0 / 2.0 );

      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      GridView gridView_;
      const Stencils &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GENERALHEIGHTFUNCTIONCURVATURE_HH
