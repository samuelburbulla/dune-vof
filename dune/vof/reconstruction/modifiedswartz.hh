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

      static constexpr int dim = GridView::dimension;

    public:
      ModifiedSwartzReconstruction ( const StencilSet &stencils, InitialReconstruction initializer,
                                     const std::size_t maxIterations = 10 )
       : stencils_( stencils ), initializer_( initializer ), maxIterations_( maxIterations )
      {}

      explicit ModifiedSwartzReconstruction ( const StencilSet &stencils, const std::size_t maxIterations = 10 )
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

        #if GRIDDIM == 1
          return;
        #endif

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
       * \brief   (local) operator application for 2d
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
      auto applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, ReconstructionSet &reconstructions ) const
       -> std::enable_if_t< Coordinate::dimension == 2, void >
      {
        auto &reconstruction = reconstructions[ entity ];

        Coordinate normal = reconstruction.innerNormal();
        std::size_t iterations = 0;
        double residuum = 0.0;

        const auto polytope = makePolytope( entity.geometry() );

        do
        {
          Coordinate oldNormal = normal;
          normal = 0;

          const auto hs = locateHalfSpace( polytope, oldNormal, color[ entity ] );
          auto interfaceEn = intersect( polytope, hs.boundary(), eager );

          for( const auto &neighbor : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor ) )
              continue;

            const auto polytopeNb = makePolytope( neighbor.geometry() );
            const auto hsNb = locateHalfSpace( polytopeNb, oldNormal, color[ neighbor ] );
            auto interfaceNb = intersect( polytopeNb, hsNb.boundary(), eager );

            Coordinate centerNormal = generalizedCrossProduct( interfaceNb.centroid() - interfaceEn.centroid() );

            if ( centerNormal * oldNormal < 0 )
              centerNormal *= -1.0;
            normalize( centerNormal );

            double weight = 1.0;
            normal.axpy( weight, centerNormal );
          }

          if ( normal.two_norm() < 0.5 )
            return;

          normalize( normal );

          assert( normal.two_norm() > 0 );

          reconstruction = locateHalfSpace( polytope, normal, color[ entity ] );

          residuum = 1.0 - ( normal * oldNormal );
          ++iterations;
        }
        while ( residuum > 1e-12 && iterations < maxIterations_ ); // residuum is always positive
      }

      /**
       * \brief   (local) operator application for 3d
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
      auto applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, ReconstructionSet &reconstructions ) const
       -> std::enable_if_t< Coordinate::dimension == 3, void >
      {
        auto &reconstruction = reconstructions[ entity ];

        Coordinate normal = reconstruction.innerNormal();
        std::size_t iterations = 0;
        double residuum = 0.0;

        const auto polytope = makePolytope( entity.geometry() );

        do
        {
          Coordinate oldNormal = normal;
          normal = 0;

          const auto hs = locateHalfSpace( polytope, oldNormal, color[ entity ] );
          auto interfaceEn = intersect( polytope, hs.boundary(), eager );

          for( const auto &neighbor : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor ) )
              continue;

            const auto polytopeNb = makePolytope( neighbor.geometry() );
            const auto hsNb = locateHalfSpace( polytopeNb, oldNormal, color[ neighbor ] );
            auto interfaceNb = intersect( polytopeNb, hsNb.boundary(), eager );

            Coordinate centerNormal( 0.0 );

            for( const auto &neighbor2 : stencil( entity ) )
            {
              if ( neighbor == neighbor2 )
                continue;

              if ( !flags.isMixed( neighbor2 ) )
                continue;

              if ( color.gridView().indexSet().index( neighbor2 ) > color.gridView().indexSet().index( neighbor ) )
                continue;

              const auto polytopeNb2 = makePolytope( neighbor2.geometry() );
              const auto hsNb2 = locateHalfSpace( polytopeNb2, oldNormal, color[ neighbor2 ] );
              auto interfaceNb2 = intersect( polytopeNb2, hsNb2.boundary(), eager );


              centerNormal = generalizedCrossProduct(
                interfaceNb.centroid()  - interfaceEn.centroid(),
                interfaceNb2.centroid() - interfaceEn.centroid()
              );

              if ( centerNormal.two_norm() < 0.5 )
                continue;

              if ( centerNormal * oldNormal < 0 )
                centerNormal *= -1.0;
              normalize( centerNormal );
            }

            double weight = 1.0;
            normal.axpy( weight, centerNormal );
          }

          if ( normal.two_norm() < 0.5 )
            return;

          normalize( normal );

          assert( normal.two_norm() > 0 );

          reconstruction = locateHalfSpace( polytope, normal, color[ entity ] );

          residuum = 1.0 - ( normal * oldNormal );
          ++iterations;
        }
        while ( residuum > 1e-12 && iterations < maxIterations_ ); // residuum is always positive
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
