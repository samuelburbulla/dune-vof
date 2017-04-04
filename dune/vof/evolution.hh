#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

#include <functional>
#include <type_traits>

//- dune-common includes
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

//- dune-grid includes
#include <dune/grid/common/partitionset.hh>

//- local includes
#include <dune/vof/common/commoperation.hh>
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
     * \tparam  GV  grid view type
     */
    template< class GV >
    struct Evolution
    {
      using GridView = GV;

    private:
      using Entity = typename GridView::template Codim< 0 >::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using ctype = typename Entity::Geometry::ctype;
      static constexpr std::size_t dim = GridView::dimension;

    public:
      explicit Evolution ( GridView gridView ) : gridView_( gridView ) {}

      /**
       * \brief (gobal) operator application
       *
       * \param   velocity        velocity
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      template< class ReconstructionSet, class Flags, class Velocity, class DiscreteFunction >
      double operator() ( const ReconstructionSet& reconstructions, const Flags& flags, Velocity& velocity, double deltaT, DiscreteFunction &update ) const
      {
        double dtEst = std::numeric_limits< double >::max();

        update.clear();

        for( const auto &entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          if( !flags.isMixed( entity ) )
            continue;

          using std::min;
          dtEst = min( dtEst, applyLocal( entity, reconstructions, flags, velocity, deltaT, update ) );
        }

        update.communicate( Dune::All_All_Interface, CommOperation::Add() );

        return gridView().comm().min( dtEst );
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
      template< class ReconstructionSet, class Flags, class Velocity, class DiscreteFunction >
      double applyLocal ( const Entity &entity,
                          const ReconstructionSet &reconstructions,
                          const Flags &flags,
                          Velocity& velocity,
                          double deltaT,
                          DiscreteFunction &update ) const
      {
        double sumFluxes = 0.0;
        double volume = entity.geometry().volume();

        for ( const auto &intersection : intersections( gridView(), entity ) )
        {
          if ( !intersection.neighbor() )
            continue;

          const Coordinate &outerNormal = intersection.centerUnitOuterNormal();

          velocity.bind( intersection );
          const auto& refElement = ReferenceElements< ctype, dim-1 >::general( intersection.type() );
          Coordinate v = velocity( refElement.position( 0, 0 ) );

          using std::abs;
          using std::min;
          sumFluxes += intersection.geometry().volume() * abs( outerNormal * v );

          v *= deltaT;

          ctype flux = 0.0;

          const auto &neighbor = intersection.outside();

          // outflow
          if ( v * outerNormal > 0 )
          {
            geometricFlux( intersection.inside(), intersection.geometry(), reconstructions, flags, v, flux );

            update[ entity ] -= flux / volume;

            if( flags.isFull( neighbor ) )
              update[ neighbor ] -= ( ( v * intersection.integrationOuterNormal( refElement.position( 0, 0 ) ) ) - flux ) / neighbor.geometry().volume();

            if ( flags.isEmpty( neighbor ) )
              update[ neighbor ] += flux / neighbor.geometry().volume();

            if ( flags.isMixed( neighbor ) )
              update[ neighbor ] += flux / neighbor.geometry().volume();
          }

          // inflow from full neighbors
          if ( flags.isFull( neighbor ) )
            if ( v * outerNormal < 0 )
              update[ entity ] -= ( v * intersection.integrationOuterNormal( refElement.position( 0, 0 ) ) ) / volume;
        }

        return volume / sumFluxes;
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
      template< class Geometry, class ReconstructionSet, class Flags >
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

      GridView gridView_;
    };

    // evolution
    // --------

    /**
     * \ingroup Method
     * \brief generate time evolution operator
     *
     * \tparam  GridView
     * \return [description]
     */
    template< class GridView >
    static inline auto evolution ( const GridView& gv )
     -> decltype( Evolution< GridView >( gv ) )
    {
      return Evolution< GridView >( gv );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
