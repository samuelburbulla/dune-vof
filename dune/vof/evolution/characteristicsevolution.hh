#ifndef DUNE_VOF_CHARACTERISTICSEVOLUTION_HH
#define DUNE_VOF_CHARACTERISTICSEVOLUTION_HH

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

    // CharacteristicsEvolution
    // ------------------------

    /**
     * \ingroup Method
     * \brief operator for time evoution
     *
     * \tparam  GV  grid view type
     */
    template< class GV >
    struct CharacteristicsEvolution
    {
      using GridView = GV;

    private:
      using Entity = typename GridView::template Codim< 0 >::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using ctype = typename Entity::Geometry::ctype;
      static constexpr std::size_t dim = GridView::dimension;

    public:
      explicit CharacteristicsEvolution ( GridView gridView ) : gridView_( gridView ) {}

      /**
       * \brief (gobal) operator application
       *
       * \param   velocity        velocity
       * \param   deltaT          delta t
       * \param   update          discrete function of flow
       */
      template< class ColorFunction, class ReconstructionSet, class Flags, class Velocity >
      double operator() ( ColorFunction &color,
                          const ReconstructionSet& reconstructions,
                          const Flags& flags,
                          const Velocity& velocity,
                          double time,
                          double deltaT ) const
      {
        double dtEst = std::numeric_limits< double >::max();

        color.clear();

        for( const auto &entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          if( !( flags.isMixed( entity ) || flags.isFull( entity ) ) )
            continue;

          using std::min;
          dtEst = min( dtEst, applyLocal( entity, reconstructions, flags, velocity, time, deltaT, color ) );
        }

        color.communicate( Dune::All_All_Interface, CommOperation::Add() );

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
      template< class ReconstructionSet, class Flags, class Velocity, class ColorFunction >
      double applyLocal ( const Entity &entity,
                          const ReconstructionSet &reconstructions,
                          const Flags &flags,
                          Velocity& velocity,
                          double time,
                          double deltaT,
                          ColorFunction &color ) const
      {
        auto polytope = makePolytope( entity.geometry() );

        if ( flags.isMixed( entity ) )
          polytope = intersect( polytope, reconstructions[ entity ] );

        for ( int i = 0; i < polytope.size(); ++i )
          getNewPosition( polytope.vertex( i ), velocity, time, deltaT );

        double controlVolume = polytope.volume();

        for ( const auto &elem : elements( gridView(), Partitions::all ) )
        {
          const auto &geo = elem.geometry();
          Polygon< Coordinate > cut = intersect( polytope, makePolytope( geo ) );
          double vol = cut.volume();
          color[ elem ] += vol / geo.volume();
          controlVolume -= vol;
        }

        assert( std::abs( controlVolume ) < 1e-14 );

        return deltaT;
      }

      /**
       * \brief Computes and returns the new postiton of a coordinate following the streamlines
       *
       * \tparam  Velocity
       * \param   corner
       * \param   velocity
       * \param   deltaT
       */
      template< class Velocity >
      void getNewPosition ( Coordinate& corner,
                            const Velocity& velocity,
                            double time,
                            double dt ) const
      {
        // Runge Kutta 4. Ordnung (Wikipedia: Klassisches Runge-Kutta-Verfahren)
        Coordinate k1 = velocity( corner, time );

        Coordinate x2 = corner;
        x2.axpy( dt / 2.0, k1 );
        Coordinate k2 = velocity( x2, time + dt / 2.0 );

        Coordinate x3 = corner;
        x3.axpy( dt / 2.0, k2 );
        Coordinate k3 = velocity( x3, time + dt / 2.0 );

        Coordinate x4 = corner;
        x4.axpy( dt, k3 );
        Coordinate k4 = velocity( x4, time + dt );

        corner.axpy( dt * 1.0 / 6.0, k1 );
        corner.axpy( dt * 1.0 / 3.0, k2 );
        corner.axpy( dt * 1.0 / 3.0, k3 );
        corner.axpy( dt * 1.0 / 6.0, k4 );
      }

      const GridView& gridView() const { return gridView_; }

      GridView gridView_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CHARACTERISTICSEVOLUTION_HH
