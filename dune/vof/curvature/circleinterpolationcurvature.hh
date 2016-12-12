#ifndef DUNE_VOF_CIRCLEINTERPOLATIONCURVATURE_HH
#define DUNE_VOF_CIRCLEINTERPOLATIONCURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell by an interpolation of a circle.

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
    struct CircleInterpolationCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit CircleInterpolationCurvature ( const GridView &gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils )
      {}

      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, DiscreteFunction &curvature )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature[ entity ] = 0.0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, flags, curvature );
        }
      }

    private:
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, DiscreteFunction &curvature )
      {
        Matrix AtA( 0.0 );
        Coordinate Atb( 0.0 );
        std::vector< Coordinate > points;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();

        for( const auto& neighbor : stencil( entity ) )
        {
          if ( !flags.isMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          points.push_back( centroidNb );
        }

        for ( const auto &point : points )
        {
          Coordinate p = point;
          p -= centroidEn;
          const ctype weight = p.two_norm2();
          p *= 2.0 * weight;
          AtA += outerProduct( p, p );
          Atb.axpy( weight * ( point.two_norm2() - centroidEn.two_norm2() ), p );
        }
        Coordinate centerOfCircle;
        AtA.solve( centerOfCircle, Atb );

        double radius = ( centerOfCircle - centroidEn ).two_norm();
        double sign = - 1 + 2 * ( 0 < ( centerOfCircle - centroidEn ) * reconstructions[ entity ].innerNormal() );
        curvature[ entity ] = - sign / radius;
      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const GridView &gridView_;
      const Stencils &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CIRCLEINTERPOLATIONCURVATURE_HH
