#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH

#include <cmath>

#include <limits>
#include <vector>

//- dune-common includes
#include <dune/common/fmatrix.hh>

#include <dune/vof/geometricutility.hh>

namespace Dune
{
  namespace VoF
  {

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

      template< class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        reconstructions.clear();
        for ( const auto &entity : elements( color.gridView() ) )
          if ( applyLocal( entity, flags, color, reconstructions[ entity ] ) )
            reconstructions.intersections( entity ) = intersections_;
      }

    private:
      template< class Flags >
      bool applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, Reconstruction &reconstruction ) const
      {
        if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
          return false;

        const auto geometry = entity.geometry();

        Coordinate &normal = reconstruction.normal();
        const Coordinate center = geometry.center();
        const auto colorEn = color[ entity ];

        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );
        for ( const auto &neighbor : stencil( entity ) )
        {
          Coordinate d = neighbor.geometry().center() - center;
          const ctype wk = 1.0 / d.two_norm2();
          d *= wk;
          AtA += outerProduct( d, d );
          Atb.axpy( wk * ( color[ neighbor ] - colorEn ), d );
        }
        AtA.solve( normal, Atb );

        if( !normalize( normal ) )
          return false;

        computeInterfaceLinePosition( geometry, colorEn, reconstruction, intersections_ );

        return true;
      }

      Matrix outerProduct ( const Vector &a, const Vector &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < dim; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      bool normalize ( Coordinate &normal ) const
      {
        const ctype length2 = normal.two_norm2();
        if ( length2 < std::numeric_limits< ctype >::epsilon() )
        {
          normal = Coordinate( 0.0 );

          return false;
        }
        normal *= 1.0 / std::sqrt( length2 );

        return true;
      }

      Stencil stencil ( const Entity &entity ) const { return stencils_[ entity ]; } // rework stencils

      StencilSet &stencils_;

      mutable std::vector< Coordinate > intersections_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
