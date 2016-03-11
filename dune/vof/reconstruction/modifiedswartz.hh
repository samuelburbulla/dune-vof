#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH

#include <cmath>
#include <numeric>
#include <vector>

#include <dune/vof/geometry/algorithm.hh>

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

        auto exchange = typename ReconstructionSet::Exchange ( reconstructions );
        color.gridView().communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );
      }

    private:
      template< class Flags >
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, ReconstructionSet &reconstructions ) const
      {
        std::size_t iterations = 0;
        Reconstruction &reconstruction = reconstructions[ entity ];
        Coordinate newNormal, &normal = reconstruction.normal();
        intersectionsEn_ = reconstructions.intersections( entity );

        const auto geoEn = entity.geometry();
        const auto &stencilEn = stencil( entity );
        do
        {
          assert( intersectionsEn_.size() != 0 );
          Coordinate centerEn = std::accumulate( intersectionsEn_.begin(), intersectionsEn_.end(), Coordinate( 0.0 ) );
          centerEn *= ( 1.0 / static_cast< typename Coordinate::value_type >( intersectionsEn_.size() ) );

          newNormal = Coordinate( 0.0 );
          for( const auto &neighbor : stencilEn )
          {
            if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
              continue;

            //if ( ( reconstructions[ neighbor ].normal() * normal ) <= 0.0 )
              //continue;

            if ( reconstructions[ neighbor ].normal().two_norm() < std::numeric_limits< double >::epsilon() )
              continue;

            Reconstruction reconstructionNb( normal, 0.0 );
            computeInterfaceLinePosition( neighbor.geometry(), color[ neighbor ], reconstructionNb, intersectionsNb_ );

            assert( intersectionsNb_.size() != 0 );
            Coordinate centerNb = std::accumulate( intersectionsNb_.begin(), intersectionsNb_.end(), Coordinate( 0.0 ) );
            centerNb *= ( 1.0 / static_cast< typename Coordinate::value_type >( intersectionsNb_.size() ) );

            Coordinate centerNormal = rotateCCW( centerNb - centerEn );
            assert( centerNormal.two_norm2() > 0.0 );
            normalize( centerNormal );

            if ( ( centerNormal * normal ) < 0.0 )
              centerNormal *= -1.0;

            if ( centerNormal * normal < std::cos( M_PI / 3.0 ) )
              continue;

            newNormal += centerNormal;
          }

          if ( newNormal == Coordinate( 0 ) )
            break;

          normalize( newNormal );

          std::swap( newNormal, normal );
          computeInterfaceLinePosition( geoEn, color[ entity ], reconstruction, intersectionsEn_ );

          ++iterations;
        }
        while ( (normal - newNormal).two_norm2() > 1e-8 && iterations < maxIterations_ );
      }

      void normalize ( Coordinate &normal ) const
      {
        normal /= normal.two_norm();
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
