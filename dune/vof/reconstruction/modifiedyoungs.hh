#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH

#include <cmath>

#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Reconstruction
     * \brief   modified Youngs reconstruction operator
     * \details Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 15ff
     *
     * \tparam DF   discrete function type
     * \tparam RS   reconstruction set type
     * \tparam StS  stencils type
     */
    template< class DF, class RS, class StS >
    struct ModifiedYoungsReconstruction
    {
      using ColorFunction = DF;
      using ReconstructionSet = RS;
      using StencilSet = StS;

      using GridView = typename ColorFunction::GridView;

    private:
      using Reconstruction = typename ReconstructionSet::Reconstruction;
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename ColorFunction::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;

      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, dim, dim >;
      using Vector = Coordinate;

    public:
      explicit ModifiedYoungsReconstruction ( StencilSet &stencils )
       : stencils_( stencils )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        reconstructions.clear();
        for ( const auto &entity : elements( color.gridView() ) )
        {
          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, flags, color, reconstructions[ entity ] );
        }

        auto exchange = typename ReconstructionSet::Exchange ( reconstructions );
        color.gridView().communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );
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
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, Reconstruction &reconstruction ) const
      {
        const auto geometry = entity.geometry();

        Coordinate normal;
        const Coordinate center = geometry.center();
        const auto colorEn = color[ entity ];
        assert( colorEn == colorEn );

        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );
        for ( const auto &neighbor : stencil( entity ) )
        {
          Coordinate d = neighbor.geometry().center() - center;
          const ctype wk = 1.0 / d.two_norm2();
          d *= wk;
          AtA += outerProduct( d, d );
          assert( color[ neighbor ] == color[ neighbor ] );
          Atb.axpy( wk * ( color[ neighbor ] - colorEn ), d );
        }
        AtA.solve( normal, Atb );

        if( normal.two_norm2() < std::numeric_limits< ctype >::epsilon() )
        {
          reconstruction = Reconstruction();
          return;
        }

        normalize( normal );

        auto polytope = makePolytope( geometry );
        reconstruction = locateHalfSpace( polytope, normal, colorEn );
      }

      Matrix outerProduct ( const Vector &a, const Vector &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < dim; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      Stencil stencil ( const Entity &entity ) const { return stencils_[ entity ]; } // rework stencils

      StencilSet &stencils_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
