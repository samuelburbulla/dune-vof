#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDSWARTZ_HH

#include <cmath>
#include <numeric>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/partitionset.hh>

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
      using GridView = GV;
      using StencilSet = StS;
      using InitialReconstruction = IR;

    private:
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      ModifiedSwartzReconstruction ( const StencilSet &stencils, InitialReconstruction initializer,
                                     const std::size_t maxIterations = 3 )
       : stencils_( stencils ), initializer_( initializer ), maxIterations_( maxIterations )
      {}

      explicit ModifiedSwartzReconstruction ( const StencilSet &stencils, const std::size_t maxIterations = 3 )
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
        auto &reconstruction = reconstructions[ entity ];

        Coordinate normal = reconstruction.innerNormal();
        std::size_t iterations = 0;
        double residuum;

        do
        {
          Coordinate oldNormal = normal;
          normal = 0;

          auto interfaceEn = interface( entity, reconstructions );

          for( const auto &intersection : intersections( color.gridView(), entity ) )
          {
            if ( !intersection.neighbor() )
              continue;

            const auto &neighbor = intersection.outside();

            if ( !flags.isMixed( neighbor ) )
              continue;

            auto interfaceNb = interface( neighbor, reconstructions );

            // TODO: 3D
          #if GRIDDIM == 2

            Coordinate centerNormal = generalizedCrossProduct( interfaceNb.centroid() - interfaceEn.centroid() );

            if ( centerNormal * oldNormal < 0 )
              centerNormal *= -1.0;
            normalize( centerNormal );

            double weight = interfaceNb.volume();
            normal.axpy( weight, centerNormal );


          #elif GRIDDIM == 3

            Coordinate centerNormal( 0.0 );

            for( const auto &intersection2 : intersections( color.gridView(), entity ) )
            {
              if ( !intersection2.neighbor() )
                continue;

              const auto &neighbor2 = intersection2.outside();

              if ( neighbor == neighbor2 )
                continue;

              if ( !flags.isMixed( neighbor2 ) )
                continue;

              auto interfaceNb2 = interface( neighbor2, reconstructions );

              centerNormal = generalizedCrossProduct(
                interfaceNb.centroid()  - interfaceEn.centroid(),
                interfaceNb2.centroid() - interfaceEn.centroid()
              );

              if ( centerNormal * oldNormal < 0 )
                centerNormal *= -1.0;
              normalize( centerNormal );

              double weight = interfaceNb.volume() * interfaceNb2.volume();
              normal.axpy( weight, centerNormal );
            }
          #endif

          }

          normalize( normal );

          assert( normal.two_norm() > 0 );

          reconstruction = locateHalfSpace( makePolytope( entity.geometry() ), normal, color[ entity ] );

          residuum = 1.0 - ( normal * oldNormal );
          ++iterations;
        }
        while ( residuum > 1e-8 && iterations < maxIterations_ );
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
