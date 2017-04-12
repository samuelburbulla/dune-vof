#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH

#include <cmath>

#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/grid/common/partitionset.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/utility.hh>

namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Reconstruction
     * \brief   modified Youngs reconstruction operator
     * \details Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 15ff
     *
     * \tparam GV   grid view
     * \tparam StS  stencils type
     */
    template< class GV, class StS >
    struct ModifiedYoungsReconstruction
    {
      using GridView = GV;
      using StencilSet = StS;

    private:
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;
      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, dim, dim >;
      using Vector = FieldVector< ctype, dim >;

    public:
      explicit ModifiedYoungsReconstruction ( const StencilSet &stencils )
       : stencils_( stencils )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  ColorFunction
       * \tparam  ReconstructionSet
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class ColorFunction, class ReconstructionSet, class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags, bool communicate = false ) const
      {
        reconstructions.clear();
        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, color, flags, reconstructions[ entity ] );
        }

        if ( communicate )
          reconstructions.communicate();
      }

      /**
       * \brief   (local) operator application
       *
       * \tparam  ColorFunction
       * \tparam  Flags
       * \tparam  Reconstruction
       * \param   entity          current element
       * \param   color           color functions
       * \param   flags           set of flags
       * \param   reconstruction  single reconstruction
       */
      template< class ColorFunction, class Flags, class Reconstruction >
      void applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, Reconstruction &reconstruction ) const
      {
        Coordinate normal( 0.0 );

        const auto geometry = entity.geometry();

        const Coordinate center = geometry.center();
        const double colorEn = color[ entity ];

        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );
        for ( const auto &neighbor : stencil( entity ) )
        {
          Vector d = neighbor.geometry().center() - center;
          const ctype weight = 1.0 / d.two_norm2();
          d *= weight;
          AtA += outerProduct( d, d );
          Atb.axpy( weight * ( clamp( color[ neighbor ], 0.0, 1.0 ) - colorEn ), d );
        }
        AtA.solve( normal, Atb );

        if( normal.two_norm2() < std::numeric_limits< ctype >::epsilon() )
          return;

        normalize( normal );

        if ( flags.isMixed( entity ) )
        {
          auto polytope = makePolytope( geometry );
          reconstruction = locateHalfSpace( polytope, normal, colorEn );
        }
        else
          reconstruction = Reconstruction( normal, 0.0 );
      }

    private:
      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const StencilSet &stencils_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
