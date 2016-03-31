#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

#include <functional>

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>


namespace Dune
{
  namespace VoF
  {

    // Evolution
    // ---------

    /**
     * \ingroup Method
     * \brief operator for time evoution
     * \details Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 24ff
     *
     * \tparam  RS  reconstructions set type
     * \tparam  DF  discrete function type
     */
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
      using Polygon_ = Polygon< Coordinate >;

      using ctype = typename ColorFunction::ctype;
    public:
      explicit Evolution ( double eps )
       : eps_( eps )
      {}

      /**
       * \brief (gobal) operator application
       *
       * \tparam  Velocity        velocity field type
       * \tparam  Flags
       * \param   color           discrete function
       * \param   reconstructions set of reconstructions
       * \param   velocity        velocity field
       * \param   dt              time step
       * \param   update          discrete function of flow
       * \param   flags           set of flags
       */
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
      /**
       * \brief (local) operator application
       *
       * \tparam  Velocity        velocity field type
       * \tparam  Flags
       * \param   entity          current element
       * \param   flags           set of flags
       * \param   dt              time step
       * \param   color           discrete function
       * \param   reconstructions set of reconstructions
       * \param   velocity        velocity field
       * \param   update          discrete function of flow
       */
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

          auto upwind = upwindPolygon( geoIs, v );

          ctype flux = 0.0;
          if ( v * outerNormal > 0 ) // outflow
          {
            if ( flags.isMixed( entity ) || flags.isFullAndMixed( entity ) )
              flux = truncVolume( upwind, reconstructions[ entity ] );
            else if ( color[ entity ] >= (1 -eps_) )
              flux = upwind.volume();
          }
          else if ( v * outerNormal < 0 ) // inflow
          {
            if ( flags.isMixed( neighbor ) || flags.isFullAndMixed( neighbor ) )
              flux = -truncVolume( upwind, reconstructions[ neighbor ] );
            else if ( color[ neighbor ] >= (1 -eps_) )
              flux = -upwind.volume();
          }

          update[ entity ] -= flux / geoEn.volume();
        }
      }

      /**
       * \brief generate upwind polygon
       *
       * \tparam  IntersectionGeometry
       * \param   iGeometry intersection geometry
       * \param   v         upwind shift
       */
      template< class IntersectionGeometry >
      inline Polygon_ upwindPolygon ( const IntersectionGeometry& iGeometry, const Coordinate& v ) const
      {
        if ( ( generalizedCrossProduct( iGeometry.corner( 1 ) - iGeometry.corner( 0 ) ) * v ) < 0.0 )
          return Polygon_( { iGeometry.corner( 0 ), iGeometry.corner( 1 ), iGeometry.corner( 1 ) - v, iGeometry.corner( 0 ) - v } );
        else
          return Polygon_( { iGeometry.corner( 1 ), iGeometry.corner( 0 ), iGeometry.corner( 0 ) - v, iGeometry.corner( 1 ) - v } );
      }

      /**
       * \brief volume of truncated upwind polygon
       */
      inline ctype truncVolume ( const Polygon_& upwind, const Reconstruction& halfSpace ) const
      {
        Polygon_ intersection = intersect( std::cref( upwind ), std::cref( halfSpace ) );
        return intersection.volume();
      }

      double eps_;
    };


    // evolution
    // --------

    /**
     * \ingroup Method
     * \brief generate time evolution operator
     *
     * \tparam  ReconstructionSet
     * \tparam  ColorFunction
     * \param   rs                  reconstruction set
     * \param   eps                 marker tolerance
     * \return [description]
     */
    template< class ReconstructionSet, class ColorFunction >
    static inline auto evolution ( const ReconstructionSet&, const ColorFunction&, double eps ) -> decltype( Evolution< ReconstructionSet, ColorFunction >( eps ) )
    {
      return Evolution< ReconstructionSet, ColorFunction >( eps );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
