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
     * \tparam  FL  flags set type
     */
    template< class RS, class FL >
    struct Evolution
    {
      using ReconstructionSet = RS;
      using Flags = FL;

    private:
      using Entity = typename ReconstructionSet::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using ctype = typename Entity::Geometry::ctype;

    public:
      explicit Evolution ( const ReconstructionSet& reconstructions,
                           const Flags& flags )
        : reconstructions_( reconstructions ),
          flags_( flags )
      {}

      /**
       * \brief (gobal) operator application
       *
       * \param   velocity        velocity
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      template< class Velocity, class DiscreteFunction >
      double operator() ( Velocity& velocity, double deltaT, DiscreteFunction &update ) const
      {
        double dtEst = std::numeric_limits< double >::max();

        update.clear();

        for( const auto &entity : elements( update.gridView(), Partitions::interiorBorder ) )
        {
          if( !flags_.isActive( entity ) )
            continue;

          using std::min;
          dtEst = min( dtEst, applyLocal( entity, reconstructions_, flags_, velocity, deltaT, update ) );
        }

        return dtEst;
      }

    private:
      /**
       * \brief (local) operator application
       *
       * \param   entity          current element
       * \param   reconstructions set of reconstructions
       * \param   flags           set of flags
       * \param   velocity        velocity field
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      template< class Velocity, class DiscreteFunction >
      double applyLocal ( const Entity &entity,
                          const ReconstructionSet &reconstructions,
                          const Flags &flags,
                          Velocity& velocity,
                          double deltaT,
                          DiscreteFunction &update ) const
      {
        double dtEst = std::numeric_limits< double >::max();
        double volume = entity.geometry().volume();

        for ( const auto &intersection : intersections( update.gridView(), entity ) )
        {
          if ( !intersection.neighbor() )
            continue;

          const Coordinate &outerNormal = intersection.centerUnitOuterNormal();

          velocity.bind( intersection );
          const auto& refElement = ReferenceElements< ctype, std::decay_t< decltype( intersection ) >::mydimension >::general( intersection.type() );
          Coordinate v = velocity( refElement.position( 0, 0 ) );

          using std::abs;
          using std::min;
          dtEst = min( dtEst, volume / ( intersection.geometry().volume() * abs( outerNormal * v ) ) );

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

        return dtEst;
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

      const ReconstructionSet& reconstructions_;
      const Flags& flags_;

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
     * \return [description]
     */
    template< class ReconstructionSet, class Flags >
    static inline auto evolution ( const ReconstructionSet& rs, const Flags& fl )
     -> decltype( Evolution< ReconstructionSet, Flags >( rs, fl ) )
    {
      return Evolution< ReconstructionSet, Flags >( rs, fl );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
