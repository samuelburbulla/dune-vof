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
     * \tparam  GV  grid view
     * \tparam  StS stencils type
     * \tparam  IR  initial reconstruction type
     */
    template< class GV, class StS, class IR >
    struct ModifiedSwartzReconstruction
    {
      using GridView =GV;
      using StencilSet = StS;
      using InitialReconstruction = IR;

    private:
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      ModifiedSwartzReconstruction ( const StencilSet &stencils, InitialReconstruction initializer,
                                     const std::size_t maxIterations = 50 )
       : stencils_( stencils ), initializer_( initializer ), maxIterations_( maxIterations )
      {}

      explicit ModifiedSwartzReconstruction ( const StencilSet &stencils, const std::size_t maxIterations = 50 )
       : stencils_( stencils ), initializer_( stencils ), maxIterations_( maxIterations )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  ColorFunction
       * \tparam  ReconstructionSet
       * \tparam  Flags           set
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class ColorFunction, class ReconstructionSet, class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        initializer()( color, reconstructions, flags );

        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, color, flags, reconstructions );
        }

        reconstructions.communicate();
      }

    private:
      /**
       * \brief   (local) operator application
       *
       * \tparam  ColorFunction
       * \tparam  Flags
       * \tparam  ReconstructionSet
       * \param   entity          current element
       * \param   flags           set of flags
       * \param   color           color functions
       * \param   reconstructions  set of reconstruction
       */
      template< class ColorFunction, class Flags, class ReconstructionSet >
      void applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, ReconstructionSet &reconstructions ) const
      {
        std::size_t iterations = 0;
        auto &reconstruction = reconstructions[ entity ];
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
            if ( !flags.isMixed( neighbor ) )
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

            if( centerNormal.two_norm2() >= std::numeric_limits< decltype( centerNormal.two_norm2() ) >::epsilon() )
              normalize( centerNormal );

            // if ( ( centerNormal * normal ) < 0.0 )
            //   centerNormal *= -1.0;

            newNormal += centerNormal;
            assert( !std::isnan( newNormal[0] ) );
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

      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }
      const InitialReconstruction &initializer () const { return initializer_; }

      const StencilSet &stencils_;
      InitialReconstruction initializer_;
      const std::size_t maxIterations_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
