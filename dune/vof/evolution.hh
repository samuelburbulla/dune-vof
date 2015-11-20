#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometricutility.hh>

namespace Dune
{
  namespace VoF
  {

    // Evolution
    // ---------

    template< class RS, class DF >
    struct Evolution
    {
      using ReconstructionSet = RS;
      using ColorFunction = DF;

      using GridView = typename ColorFunction::GridView;

    private:
      using Reconstruction = typename ReconstructionSet::Reconstruction;
      using Entity = typename ReconstructionSet::Entity;

      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using Polygon = Polygon2D< Coordinate >;

      using ctype = typename ColorFunction::ctype;
    public:
      explicit Evolution ( double eps )
       : eps_( eps )
      {}

      template< class Velocity, class Flags >
      void operator() ( const ColorFunction &color, const ReconstructionSet &reconstructions, const Velocity& velocity, const double dt,
                        ColorFunction &update, const Flags &flags ) const
      {
        update.clear();

        for( const auto &entity : elements( color.gridView() ) )
        {
          if( !flags.isMixed( entity ) && !flags.isActive( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, flags, dt, color, reconstructions, velocity, update);
        }
      }

    private:
      template< class Velocity, class Flags >
      void applyLocal ( const Entity &entity, const Flags &flags, const double dt, const ColorFunction &color, const ReconstructionSet &reconstructions,
                        const Velocity& velocity, ColorFunction &update ) const
      {
        const auto geoEn = entity.geometry();

        for ( const auto &intersection : intersections( color.gridView(), entity ) )
        {
          if ( !intersection.neighbor() )
            continue;

          const auto geoIs = intersection.geometry();

          const Coordinate outerNormal = intersection.centerUnitOuterNormal();
          Coordinate v = velocity( geoIs.center() );
          v *= dt;

          const auto &neighbor = intersection.outside();

          upwind_.clear();
          flux_.clear();

          ctype flux = 0.0;
          if ( v * outerNormal > 0 ) // outflow
          {
            // insert in cw/ccw order
            upwind_.addVertex( geoIs.corner( 0 ) );
            upwind_.addVertex( geoIs.corner( 1 ) );
            upwind_.addVertex( geoIs.corner( 1 ) - v );
            upwind_.addVertex( geoIs.corner( 0 ) - v );

            if ( flags.isMixed( entity ) || flags.isFullAndMixed( entity ) )
              flux = truncVolume( entity, reconstructions );
            else if ( color[ entity ] >= (1 -eps_) )
              flux = upwind_.volume();
          }
          else if ( v * outerNormal < 0 ) // inflow
          {
            // insert in cw/ccw order
            upwind_.addVertex( geoIs.corner( 0 ) );
            upwind_.addVertex( geoIs.corner( 0 ) - v );
            upwind_.addVertex( geoIs.corner( 1 ) - v );
            upwind_.addVertex( geoIs.corner( 1 ) );

            if ( flags.isMixed( neighbor ) || flags.isFullAndMixed( neighbor ) )
              flux = -truncVolume( neighbor, reconstructions );
            else if ( color[ neighbor ] >= (1 -eps_) )
              flux = -upwind_.volume();
          }

          update[ entity ] -= flux / geoEn.volume();
        }
      }

      ctype truncVolume ( const Entity &entity, const ReconstructionSet &reconstructions ) const
      {
        polygonLineIntersection( upwind_, reconstructions[ entity ], flux_ );
        polyAddInnerVertices( upwind_, reconstructions[ entity ], flux_ );

        return flux_.volume();
      }

      double eps_;
      mutable Polygon flux_, upwind_;
    };


    // evolution
    // --------

    template< class ReconstructionSet, class ColorFunction >
    static inline auto evolution ( const ReconstructionSet&, const ColorFunction&, double eps ) -> decltype( Evolution< ReconstructionSet, ColorFunction >( eps ) )
    {
      return Evolution< ReconstructionSet, ColorFunction >( eps );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
