#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH

#include <cmath>
#include <numeric>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/polygon.hh>
#include <dune/vof/geometry/polytope.hh>
#include <dune/vof/geometry/utility.hh>

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
        Coordinate newNormal, normal = reconstruction.innerNormal();

        const auto geoEn = entity.geometry();
        const auto& stencilEn = stencil( entity );
        auto polygonEn = make_polygon( geoEn );

        do
        {
          Line< Coordinate > lineEn = intersect( std::cref( polygonEn ), reconstruction.boundary() );

          newNormal = Coordinate( 0 );
          for( const auto &neighbor : stencilEn )
          {
            if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
              continue;

            //if ( ( reconstructions[ neighbor ].normal() * normal ) <= 0.0 )
              //continue;

            if ( reconstructions[ neighbor ].innerNormal().two_norm() < std::numeric_limits< double >::epsilon() )
              continue;

            const auto& polygonNb = make_polygon( neighbor.geometry() );
            Line< Coordinate > lineNb = intersect( std::cref( polygonNb ), locateHalfSpace( polygonNb, normal, color[ neighbor ] ).boundary() );

            Coordinate centerNormal = generalizedCrossProduct( lineNb.centroid() - lineEn.centroid() );
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
