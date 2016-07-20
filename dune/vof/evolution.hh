#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

#include <functional>
#include <type_traits>

//- dune-common includes
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

//- local includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/upwindpolygon.hh>
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
      using Polytope = typename std::conditional< Coordinate::dimension == 2, Polygon< Coordinate >, Polyhedron< Coordinate > >::type;

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
       * \param   timeProvider    time provider
       * \param   update          discrete function of flow
       * \param   flags           set of flags
       */
      template< class Velocity, class Flags, class TimeProvider >
      double operator() ( const ColorFunction &color, const ReconstructionSet &reconstructions, const Velocity& velocity, TimeProvider &timeProvider,
                        ColorFunction &update, const Flags &flags ) const
      {
        double elapsedTime = - MPI_Wtime();

        update.clear();

        for( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if( !flags.isMixed( entity ) && !flags.isActive( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, flags, timeProvider, color, reconstructions, velocity, update);
        }

        elapsedTime += MPI_Wtime();
        return elapsedTime;
      }

    private:
      /**
       * \brief (local) operator application
       *
       * \tparam  Velocity        velocity field type
       * \tparam  Flags
       * \param   entity          current element
       * \param   flags           set of flags
       * \param   timeProvider    time provider
       * \param   color           discrete function
       * \param   reconstructions set of reconstructions
       * \param   velocity        velocity field
       * \param   update          discrete function of flow
       */
      template< class Velocity, class Flags, class TimeProvider >
      void applyLocal ( const Entity &entity, const Flags &flags, TimeProvider &timeProvider, const ColorFunction &color, const ReconstructionSet &reconstructions,
                        const Velocity& velocity, ColorFunction &update ) const
      {
        const auto geoEn = entity.geometry();

        for ( const auto &intersection : intersections( color.gridView(), entity ) )
        {
          const auto geoIs = intersection.geometry();

          const Coordinate outerNormal = intersection.centerUnitOuterNormal();
          Coordinate v = velocity( geoIs.center() );

          double dtEst = geoEn.volume() / std::abs( intersection.integrationOuterNormal( typename decltype( geoIs )::LocalCoordinate( 0 ) ) * v );
          timeProvider.provideTimeStepEstimate( dtEst );

          v *= timeProvider.deltaT();

          auto upwind = upwindPolygon( geoIs, v );

          ctype flux = 0.0;
          if ( v * outerNormal > 0 ) // outflow
          {
            if ( flags.isMixed( entity ) || flags.isFullAndMixed( entity ) )
              flux = truncVolume( upwind, reconstructions[ entity ] );
            else if ( color[ entity ] >= ( 1 - eps_ ) )
              flux = upwind.volume();
          }
          else if ( v * outerNormal < 0 ) // inflow
          {
            if ( !intersection.neighbor() )
              continue;

            const auto &neighbor = intersection.outside();

            if ( flags.isMixed( neighbor ) || flags.isFullAndMixed( neighbor ) )
              flux = -truncVolume( upwind, reconstructions[ neighbor ] );
            else if ( color[ neighbor ] >= ( 1 - eps_ ) )
              flux = -upwind.volume();
          }
          assert( flux == flux );

          update[ entity ] -= flux / geoEn.volume();
        }
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
