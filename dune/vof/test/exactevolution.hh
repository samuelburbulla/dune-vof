#ifndef DUNE_VOF_EXACTEVOLUTION_HH
#define DUNE_VOF_EXACTEVOLUTION_HH

#include <functional>
#include <type_traits>

//- dune-common includes
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

//- local includes
#include "utility.hh"
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
     *
     * \tparam  DF  discrete function type
     */
    template< class DF >
    struct ExactEvolution
    {
      using ColorFunction = DF;

      using GridView = typename ColorFunction::GridView;

    private:
      using Entity = typename GridView::template Codim< 0 >::Entity;

      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using Polytope = typename std::conditional< Coordinate::dimension == 2, Polygon< Coordinate >, Polyhedron< Coordinate > >::type;

      using ctype = typename ColorFunction::ctype;
    public:
      explicit ExactEvolution () {}

      /**
       * \brief (gobal) operator application
       *
       * \tparam  Velocity        velocity field type
       * \tparam  Flags
       * \param   color           discrete function
       * \param   circle          circle from exact solution
       * \param   velocity        velocity field
       * \param   timeProvider    time provider
       * \param   update          discrete function of flow
       * \param   flags           set of flags
       */
      template< class Circle, class Velocity, class Flags, class TimeProvider >
      double operator() ( const ColorFunction &color, const Circle& circle, const Velocity& velocity, TimeProvider &timeProvider,
                        ColorFunction &update, const Flags &flags ) const
      {
        double elapsedTime = - MPI_Wtime();

        update.clear();

        for( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
          applyLocal( entity, flags, timeProvider, color, circle, velocity, update);

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
       * \param   circle          circle from exact solution
       * \param   velocity        velocity field
       * \param   update          discrete function of flow
       */
      template< class Velocity, class Flags, class TimeProvider, class Circle >
      void applyLocal ( const Entity &entity, const Flags &flags, const TimeProvider &timeProvider, const ColorFunction &color, const Circle &circle,
                        const Velocity& velocity, ColorFunction &update ) const
      {
        const auto geoEn = entity.geometry();

        for ( const auto &intersection : intersections( color.gridView(), entity ) )
        {
          const auto geoIs = intersection.geometry();

          const Coordinate outerNormal = intersection.centerUnitOuterNormal();
          Coordinate v = velocity( geoIs.center(), timeProvider.time() + 0.5 * timeProvider.deltaT() );
          v *= timeProvider.deltaT();

          auto upwind = upwindPolygon( geoIs, v );

          const double t = timeProvider.time();
          const Coordinate c = circle.center( t );
          const double r = circle.radius( t );

          ctype flux = 0.0;
          if ( v * outerNormal > 0 ) // outflow
            flux = -intersectionVolume( upwind, c, r );
          else if ( v * outerNormal < 0 ) // inflow
            flux = intersectionVolume( upwind, c, r );

          update[ entity ] += flux / geoEn.volume();
        }
      }
    };

    // exactEvolution
    // --------------

    /**
     * \ingroup Method
     * \brief generate time evolution operator
     *
     * \tparam  ColorFunction
     * \return [description]
     */
    template< class ColorFunction >
    static inline auto exactEvolution ( const ColorFunction& ) -> decltype( ExactEvolution< ColorFunction >() )
    {
      return ExactEvolution< ColorFunction >();
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EXACTEVOLUTION_HH
