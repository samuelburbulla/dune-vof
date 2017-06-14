#ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
#define DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-vof includes
#include <dune/vof/dataset.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/stencil/heightfunctionstencil.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>
#include <dune/vof/reconstruction/heightfunction.hh>
#include <dune/vof/utility.hh>


namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell by the height function method.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV    grid view
     * \tparam  VNST  vertex neighbor stencils
     */
    template< class GV, class VNST >
    struct CartesianHeightFunctionCurvature
     : public HeightFunctionReconstruction< GV, VNST, ModifiedYoungsReconstruction< GV, VNST > >
    {
    private:
      using BaseType = HeightFunctionReconstruction< GV, VNST, ModifiedYoungsReconstruction< GV, VNST > >;
    public:
      using GridView = typename BaseType::GridView;
      using StencilSet = typename BaseType::StencilSet;
      using Stencil = Dune::VoF::HeightFunctionStencil< GridView >;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      static constexpr int dim = GridView::dimension;

      using Heights = Dune::FieldVector< double, Stencil::noc >;
      using Orientation = std::tuple< int, int >;

    public:
      explicit CartesianHeightFunctionCurvature ( const StencilSet &vertexNeighborStencils )
        : BaseType( vertexNeighborStencils ) {}

      template< class ColorFunction, class ReconstructionSet, class Flags, class CurvatureSet >
      void operator() ( const ColorFunction &color, const ReconstructionSet &reconstructions, const Flags &flags,
                        CurvatureSet &curvature, bool communicate = false ) const
      {
        for ( const auto& entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          curvature[ entity ] = 0.0;
          satisfiesConstraint_[ entity ] = 0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, color, reconstructions, curvature );
        }

        curvature.communicate();
        satisfiesConstraint_.communicate();

        for ( int i = 0; i < 3; ++i )
        {
          CurvatureSet newCurvature ( color.gridView() );

          for ( const auto& entity : elements( color.gridView(), Partitions::interiorBorder ) )
          {
            if ( !flags.isMixed( entity ) )
              continue;

            averageCurvature( entity, curvature, flags, newCurvature, (i > 0) );
          }
          curvature = newCurvature;
        }

        if ( communicate )
          curvature.communicate();
      }

    private:
      using BaseType::satisfiesConstraint_;
      using BaseType::vertexStencil;

      template< class ColorFunction, class ReconstructionSet, class CurvatureSet >
      void applyLocal ( const Entity &entity, const ColorFunction &color, const ReconstructionSet &reconstructions,
                        CurvatureSet &curvature ) const
      {
        const Coordinate &normal = reconstructions[ entity ].innerNormal();
        const Orientation orientation = this->getOrientation( normal );

        const auto entityInfo = GridView::Grid::getRealImplementation( entity ).entityInfo();
        const Stencil stencil ( color.gridView(), entityInfo, orientation, 1 );

        Heights heights = this->getHeightValues( color, stencil );

        // Constraint
        double uMid = heights[ ( heights.size() - 1 ) / 2 ];
        int effTdown = stencil.effectiveTdown();
        if ( uMid < effTdown || uMid > effTdown + 1 )
          return;

        double deltaX = std::pow( entity.geometry().volume(), 1.0 / dim );

        for( std::size_t i = 0; i < decltype( stencil )::noc; ++i )
          if ( heights[ i ] == 0.0 )
            return;

        satisfiesConstraint_[ entity ] = 1;

        curvature[ entity ] = kappa( heights, deltaX );
      }

    # if GRIDDIM == 1
      double kappa ( const Heights &heights, const double dx ) const
      {
        return 0;
      }
    # elif GRIDDIM == 2
      double kappa ( const Heights &heights, const double dx ) const
      {
        double Hx = ( heights[ 2 ] - heights[ 0 ] ) / 2.0;
        double Hxx = ( heights[ 2 ] - 2 * heights[ 1 ] + heights[ 0 ] ) / dx;

        return - Hxx / std::pow( 1.0 + Hx * Hx, 3.0 / 2.0 );
      }
    #elif GRIDDIM == 3
      double kappa ( const Heights &heights, const double dx ) const
      {
        double Hx = ( heights[ 5 ] - heights[ 3 ] ) / 2.0;
        double Hy = ( heights[ 7 ] - heights[ 1 ] ) / 2.0;
        double Hxx = ( heights[ 5 ] - 2.0 * heights[ 4 ] + heights[ 3 ] ) / dx;
        double Hyy = ( heights[ 7 ] - 2.0 * heights[ 4 ] + heights[ 1 ] ) / dx;
        double Hxy = ( heights[ 8 ] - heights[ 2 ] - heights[ 6 ] + heights[ 0 ] ) / ( 4.0 * dx );

        return - ( Hxx + Hyy + Hxx * Hy * Hy + Hyy * Hx * Hx - 2.0 * Hxy * Hx * Hy ) / ( std::pow( 1.0 + Hx * Hx + Hy * Hy, 3.0 / 2.0 ) );
      }
    #endif

      template< class CurvatureSet, class Flags >
      void averageCurvature( const Entity &entity, const CurvatureSet &curvature, const Flags &flags, CurvatureSet &newCurvature, bool ignore ) const
      {
        int n = 0;

        if ( satisfiesConstraint_[ entity ] || ignore )
        {
          newCurvature[ entity ] += curvature[ entity ];
          n++;
        }

        for( const auto& neighbor : vertexStencil( entity ) )
        {
          if ( !flags.isMixed( neighbor ) )
            continue;

          if ( satisfiesConstraint_[ neighbor ] || ignore )
          {
            newCurvature[ entity ] += curvature[ neighbor ];
            n++;
          }
        }

        if ( n > 0 )
          newCurvature[ entity ] /= static_cast< double >( n );
      }
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
