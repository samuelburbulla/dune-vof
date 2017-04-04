#ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGSSECONDORDER_HH
#define DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGSSECONDORDER_HH

#include <cmath>

#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>

namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Reconstruction
     * \brief   modified Youngs second order reconstruction operator
     * \details Denner F., Van Wachem B. G. M., Balanced-Force VoF Framework For Arbitrary Meshes, p. 231f
     *
     * \tparam DF   discrete function type
     * \tparam RS   reconstruction set type
     * \tparam StS  stencils type
     */
    template< class GV, class StS >
    struct ModifiedYoungsSecondOrderReconstruction
    {
      using GridView = GV;
      using StencilSet = StS;

    private:
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;
      static constexpr int derivatives = dim + dim * ( dim + 1 ) / 2;
      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, derivatives, derivatives >;
      using Vector = FieldVector< ctype, derivatives >;

      using FirstOrderReconstruction = ModifiedYoungsReconstruction< GridView, StencilSet >;

    public:
      explicit ModifiedYoungsSecondOrderReconstruction ( const StencilSet &stencils )
       : stencils_( stencils ), firstOrderReconstruction( stencils )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class ColorFunction, class ReconstructionSet, class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        reconstructions.clear();
        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
          applyLocal( entity, flags, color, reconstructions[ entity ] );

        reconstructions.communicate();
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
      template< class Flags, class ColorFunction, class Reconstruction >
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, Reconstruction &reconstruction ) const
      {
        Coordinate normal( 0.0 );

        const auto geometry = entity.geometry();

        const Coordinate center = geometry.center();
        const double colorEn = color[ entity ];

        // second order
        if ( stencil( entity ).size() >= std::pow( 3, dim ) - 1  )
        {
          Matrix AtA( 0.0 );
          Vector Atb( 0.0 );
          for ( const auto &neighbor : stencil( entity ) )
          {
            Coordinate d = neighbor.geometry().center() - center;

            Vector v; // v = ( dx,  dy, dxx, dyy, dxy )
            for ( std::size_t i = 0; i < dim; ++i )
              v[ i ] = d[ i ];

            for ( std::size_t i = 0; i < dim; ++i )
              v[ i + dim ] = 0.5 * d[ i ] * d[ i ];

            int n = 0;
            for ( std::size_t i = 0; i < dim; ++i )
              for ( std::size_t j = i+1; j < dim; ++j )
                if ( i == j )
                  continue;
                else
                {
                  v[ n + 2 * dim ] = d[ i ] * d[ j ];
                  n++;
                }

            const ctype weight = 1.0 / d.two_norm2();
            v *= weight;
            AtA += outerProduct( v, v );
            Atb.axpy( weight * ( color[ neighbor ] - colorEn ), v );
          }

          Vector solution( 0.0 );
          AtA.solve( solution, Atb );

          for ( std::size_t i = 0; i < dim; ++i )
            normal[ i ] = solution[ i ];


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
        else
        {
          firstOrderReconstruction.applyLocal( entity, flags, color, reconstruction );
        }
      }

      Matrix outerProduct ( const Vector &a, const Vector &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < derivatives; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const StencilSet &stencils_;
      const FirstOrderReconstruction firstOrderReconstruction;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGSSECONDORDER_HH
