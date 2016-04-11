#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH

#include <cmath>
#include <numeric>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    // ModifiedSwartzReconstruction
    // ----------------------------


    /**
     * \ingroup   Reconstruction
     * \brief     modified Swartz reconstruction operator
     * \details   Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 15ff
     *
     * \tparam  DF  discrete function type
     * \tparam  RS  reconstruction set type
     * \tparam  StS stencils type
     * \tparam  IR  initial reconstruction type
     */
    template< class DF, class RS, class StS, class IR >
    struct ModifiedSwartzReconstruction
    {
      using ColorFunction = DF;
      using ReconstructionSet = RS;
      using StencilSet = StS;
      using InitialReconstruction = IR;

      using GridView = typename ColorFunction::GridView;

    private:
      using Reconstruction = typename ReconstructionSet::Reconstruction;
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename ColorFunction::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      ModifiedSwartzReconstruction ( StencilSet &stencils, InitialReconstruction initializer,
                                     const std::size_t maxIterations = 50 )
       : stencils_( stencils ), initializer_( initializer ), maxIterations_( maxIterations )
      {}

      explicit ModifiedSwartzReconstruction ( StencilSet &stencils, const std::size_t maxIterations = 50 )
       : stencils_( stencils ), initializer_( stencils ), maxIterations_( maxIterations )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  Flags           set
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class Flags >
      double operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        double elapsedTime = initializer()( color, reconstructions, flags );
        elapsedTime = - MPI_Wtime();

        for ( const auto &entity : elements( color.gridView() ) )
        {
          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, flags, color, reconstructions );
        }
        elapsedTime += MPI_Wtime();

        auto exchange = typename ReconstructionSet::Exchange ( reconstructions );
        color.gridView().communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );

        return elapsedTime;
      }

    private:
      /**
       * \brief   (local) operator application
       *
       * \tparam  Flags
       * \param   entity          current element
       * \param   flags           set of flags
       * \param   color           color functions
       * \param   reconstructions set of reconstruction
       */
      template< class Flags >
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, ReconstructionSet &reconstructions ) const
      {

        std::size_t iterations = 0;
        Reconstruction &reconstruction = reconstructions[ entity ];
        Coordinate newNormal, normal = reconstruction.innerNormal();

        const auto geoEn = entity.geometry();
        const auto& stencilEn = stencil( entity );
        auto polygonEn = makePolytope( geoEn );

        do
        {
          auto it1 = intersect( std::cref( polygonEn ), reconstruction.boundary() );
          auto lineEn = static_cast< typename decltype( it1 )::Result > ( it1 );

          newNormal = Coordinate( 0 );
          for( const auto &neighbor : stencilEn )
          {
            if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
              continue;

            // disregard empty neighbors
            if ( reconstructions[ neighbor ].innerNormal() == Coordinate( 0 ) )
              continue;

            auto polygonNb = makePolytope( neighbor.geometry() );
            auto it2 = intersect( std::cref( polygonNb ), locateHalfSpace( polygonNb, normal, color[ neighbor ] ).boundary() );
            auto lineNb = static_cast< typename decltype( it2 )::Result > ( it2 );

            Coordinate direction = lineNb.centroid() - lineEn.centroid();
            Coordinate centerNormal = normal;
            centerNormal.axpy( -(normal * direction) / direction.two_norm2(), direction );

            // Coordinate centerNormal = generalizedCrossProduct( lineNb.centroid() - lineEn.centroid() );
            // assert( centerNormal.two_norm2() > 0.0 );
            normalize( centerNormal );

            // if ( ( centerNormal * normal ) < 0.0 )
            //   centerNormal *= -1.0;

            // if ( centerNormal * normal < std::cos( M_PI / 3.0 ) )
            //   continue;


            newNormal += centerNormal;
          }

          if ( newNormal == Coordinate( 0 ) )
            break;

          normalize( newNormal );

          reconstruction = locateHalfSpace( polygonEn, newNormal, color[ entity ] );

          std::swap( newNormal, normal );
          ++iterations;
        }
        while ( (normal - newNormal).two_norm2() > 1e-8 && iterations < maxIterations_ );
      }

      Stencil stencil ( const Entity &entity ) const { return stencils_[ entity ]; } // rework stencil
      const InitialReconstruction &initializer () const { return initializer_; }

      StencilSet &stencils_;
      InitialReconstruction initializer_;
      const std::size_t maxIterations_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
