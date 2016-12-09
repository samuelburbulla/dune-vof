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
     * \tparam  DF  discrete function type
     * \tparam  RS  reconstructions set type
     * \tparam  FL  flags set type
     * \tparam  VE  velocity field type
     */
    template< class DF, class RS, class FL, class VE >
    struct Evolution
    {
      using ReconstructionSet = RS;
      using Flags = FL;
      using DiscreteFunction = DF;
      using Velocity = VE;

      using GridView = typename DiscreteFunction::GridView;

    private:
      using Entity = typename ReconstructionSet::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using ctype = typename DiscreteFunction::ctype;

    public:
      explicit Evolution ( const DiscreteFunction& update,
                           const ReconstructionSet& reconstructions,
                           const Flags& flags,
                           Velocity& velocity )
        : gridView_( update.gridView() ),
          reconstructions_( reconstructions ),
          flags_( flags ),
          velocity_( velocity )
      {}

      /**
       * \brief (gobal) operator application
       *
       * \param   t               time
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      void operator() ( double t, double deltaT, DiscreteFunction &update ) const
      {
        dtEst_ = std::numeric_limits< double >::max();

        update.clear();

        for( const auto &entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          if( !flags_.isActive( entity ) )
            continue;

          applyLocal( entity, reconstructions_, flags_, velocity_, t, deltaT, update );
        }
      }

      double provideTimeStepEstimate () const
      {
        return dtEst_;
      }

    private:
      /**
       * \brief (local) operator application
       *
       * \param   entity          current element
       * \param   reconstructions set of reconstructions
       * \param   flags           set of flags
       * \param   velocity        velocity field
       * \param   t               time
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      void applyLocal ( const Entity &entity,
                        const ReconstructionSet &reconstructions,
                        const Flags &flags,
                        Velocity& velocity,
                        double t,
                        double deltaT,
                        DiscreteFunction &update ) const
      {
        double volume = entity.geometry().volume();

        for ( const auto &intersection : intersections( gridView(), entity ) )
        {
          if ( !intersection.neighbor() )
            continue;

          const Coordinate &outerNormal = intersection.centerUnitOuterNormal();

          velocity.bind( intersection );
          Coordinate v = velocity( intersection.geometry().center(), t + 0.5 * deltaT );

          dtEst_ = std::min( dtEst_, volume / std::abs( intersection.integrationOuterNormal( 0.0 ) * v ) );

          v *= deltaT;

          ctype flux = 0.0;

          // outflow
          if ( v * outerNormal > 0 )
          {
            geometricFlux( intersection.inside(), intersection.geometry(), reconstructions, flags, v, flux );
            flux *= -1.0;
          }
          // inflow
          else if ( v * outerNormal < 0 )
            geometricFlux( intersection.outside(), intersection.geometry(), reconstructions, flags, v, flux );

          update[ entity ] += flux / volume;

        }
      }

      /**
       * \brief Computes and returns the geometric flux
       * \details [long description]
       *
       * \tparam  Geometry
       * \param   entity
       * \param   geometry
       * \param   reconstructions
       * \param   flags
       * \param   v
       * \param   flux
       */
      template< class Geometry >
      void geometricFlux ( const Entity& entity,
                           const Geometry& geometry,
                           const ReconstructionSet& reconstructions,
                           const Flags& flags,
                           const Coordinate& v,
                           ctype& flux ) const
      {
        flux = 0.0;
        auto upwind = upwindPolygon( geometry, v );

        if ( flags.isMixed( entity ) )
          flux = truncVolume( upwind, reconstructions[ entity ] );
        else if ( flags.isFull( entity ) )
          flux = upwind.volume();
      }

      const GridView& gridView() const { return gridView_; }

      const GridView& gridView_;
      const ReconstructionSet& reconstructions_;
      const Flags& flags_;
      Velocity& velocity_;

      mutable double dtEst_;
    };

    // evolution
    // --------

    /**
     * \ingroup Method
     * \brief generate time evolution operator
     *
     * \tparam  DiscreteFunction
     * \tparam  ReconstructionSet
     * \tparam  Flags
     * \tparam  Velocity
     * \return [description]
     */
    template< class DiscreteFunction, class ReconstructionSet, class Flags, class Velocity >
    static inline auto evolution ( const DiscreteFunction& df, const ReconstructionSet& rs, const Flags& fl, Velocity& ve )
     -> decltype( Evolution< DiscreteFunction, ReconstructionSet, Flags, Velocity >( df, rs, fl, ve ) )
    {
      return Evolution< DiscreteFunction, ReconstructionSet, Flags, Velocity >( df, rs, fl, ve );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
