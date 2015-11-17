#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH

#include <numeric>
#include <vector>

#include <dune/vof/geometricutility.hh>

namespace Dune
{
  namespace VoF
  {

    // maybe further "dunified" with trait classes ...

    // ModifiedSwartzReconstruction
    // ----------------------------

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
                                     const std::size_t maxIterations = 30 )
       : stencils_( stencils ), initializer_( initializer ), maxIterations_( maxIterations )
      {}

      explicit ModifiedSwartzReconstruction ( StencilSet &stencils, const std::size_t maxIterations = 30 )
       : stencils_( stencils ), initializer_( stencils ), maxIterations_( maxIterations )
      {}


      template< class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        initializer()( color, reconstructions, flags );

        for ( const auto &entity : elements( color.gridView() ) )
        {
          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, flags, color, reconstructions );
          reconstructions.intersections( entity ) = intersectionsEn_;
        }
      }

    private:
      template< class Flags >
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, ReconstructionSet &reconstructions ) const
      {
        std::size_t iterations = 0;
        Reconstruction &reconstruction = reconstructions[ entity ];
        Coordinate &normal = reconstruction.normal();
        intersectionsEn_ = reconstructions.intersections( entity );
        Coordinate newNormal;

        const auto geoEn = entity.geometry();
        const auto &stencilEn = stencil( entity );
        do
        {
          assert( intersectionsEn_.size() != 0 );
          Coordinate centerEn = std::accumulate( intersectionsEn_.begin(), intersectionsEn_.end(), Coordinate( 0.0 ) );
          centerEn *= ( 1.0 / static_cast< typename Coordinate::value_type >( intersectionsEn_.size() ) );

          newNormal = Coordinate( 0.0 );

          std::size_t count = 0;

          for( const auto &neighbor : stencilEn )
          {
            if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
              continue;

            Reconstruction reconstructionNb = reconstructions[ neighbor ];

            if ( ( reconstructionNb.normal() * normal ) <= 0.0 )
              continue;

            reconstructionNb.normal() = normal;
            const auto geoNb = neighbor.geometry();

            computeInterfaceLinePosition( geoNb, color[ neighbor ], reconstructionNb, intersectionsNb_ );

            assert( intersectionsNb_.size() != 0 );
            Coordinate centerNb = std::accumulate( intersectionsNb_.begin(), intersectionsNb_.end(), Coordinate( 0.0 ) );
            centerNb *= ( 1.0 / static_cast< typename Coordinate::value_type >( intersectionsNb_.size() ) );

            Coordinate centerNormal = rotate90degreesCounterClockwise( centerNb - centerEn );
            assert( centerNormal.two_norm2() > 0.0 ); // this assertion is fishy !

            if ( ( centerNormal * normal ) < 0.0 )
              centerNormal *= -1.0;

            newNormal.axpy( 1.0 / centerNormal.two_norm() , centerNormal );
            ++count;
          }

          if( count != 0 )
          {
            newNormal *= 1.0 / static_cast< typename Coordinate::value_type >( count );
            std::swap( newNormal, normal );
            computeInterfaceLinePosition( geoEn, color[ entity ], reconstruction, intersectionsEn_ );
          }
          else
            break;

          ++iterations;
        }
        while ( (normal - newNormal).two_norm2() > 1e-8 && iterations < maxIterations_ );
      }

      void normalize ( Coordinate &normal ) const
      {
        const auto length2 = normal.two_norm2();
        assert( length2 > std::numeric_limits< decltype( length2 ) >::epsilon() );
        normal *= 1.0 / std::sqrt( length2 );
      }

      Stencil stencil ( const Entity &entity ) const { return stencils_[ entity ]; } // rework stencil
      const InitialReconstruction &initializer () const { return initializer_; }

      StencilSet &stencils_;
      InitialReconstruction initializer_;
      const std::size_t maxIterations_;

      mutable std::vector< Coordinate > intersectionsEn_;
      mutable std::vector< Coordinate > intersectionsNb_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
