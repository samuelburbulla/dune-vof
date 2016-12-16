#ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
#define DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH

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

    // Calculate the curvature of the inferface in each mixed cell by the height function method.

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
    struct CartesianHeightFunctionCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit CartesianHeightFunctionCurvature ( const GridView &gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils )
      {}

      template< class CurvatureSet >
      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature[ entity ] = 0.0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, curvature );
        }

        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          if ( curvature[ entity ] != 0.0 )
            continue;

          averageCurvature( entity, curvature ); // TODO: only use cells with fulfill the necessary condition for curvature calculation
        }
      }

    private:

      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, CurvatureSet &curvature )
      {
        double tol = std::numeric_limits< double >::epsilon();

        Coordinate normal = reconstructions[ entity ].innerNormal();
        Coordinate center = entity.geometry().center();

        Coordinate direction ( 0.0 );
        Coordinate directionOrth ( 0.0 );
        if ( std::abs( normal[ 0 ] ) > std::abs( normal[ 1 ] ) )
        {
         direction[ 0 ] = 1.0;
         directionOrth[ 1 ] = 1.0;
        }
        else
        {
         direction[ 1 ] = 1.0;
         directionOrth[ 0 ] = 1.0;
        }

        double left = 0.0, right = 0.0;
        double middle = uh[ entity ];
        double dx = std::numeric_limits< double >::max();
        for( const auto& neighbor : stencil( entity ) )
        {
          Coordinate d = center - neighbor.geometry().center();

          if ( std::abs( d * directionOrth ) < tol )
            middle += uh[ neighbor ];
          else if ( d * directionOrth > 0 )
            left += uh[ neighbor ];
          else if ( d * directionOrth < 0 )
            right += uh[ neighbor ];

          dx = std::min( dx, d.two_norm() );
        }

        right *= dx;
        left *= dx;
        middle *= dx;

        if ( dx < middle && middle < 2 * dx )
        {
          // Workaround for boundary: linear interpolation
          if ( right < tol )
            right = middle + ( middle - left );
          if ( left < tol )
            left = middle + ( middle - right );

          double Hx = ( right - left ) / ( 2 * dx );
          double Hxx = ( right - 2 * middle + left ) / ( dx * dx );

          curvature[ entity ] = Hxx / std::pow( 1.0 + Hx * Hx, 3.0 / 2.0 );
        }
      }

      template< class CurvatureSet >
      void averageCurvature( const Entity &entity, CurvatureSet &curvature )
      {
        int n = 0;
        for( const auto& neighbor : stencil( entity ) )
        {
          curvature[ entity ] += curvature[ neighbor ];
          n++;
        }
        curvature[ entity ] /= n;
      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const GridView &gridView_;
      const Stencils &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
