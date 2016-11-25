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
      static constexpr int derivatives = dim + dim * ( dim + 1 ) / 2;
      using ctype = typename Coordinate::value_type;
      using MatrixD = FieldMatrix< ctype, derivatives, derivatives >;
      using VectorD = FieldVector< ctype, derivatives >;
      using Matrix = FieldMatrix< ctype, dim, dim >;
      using Vector = FieldVector< ctype, dim >;

    public:
      explicit ModifiedYoungsReconstruction ( const StencilSet &stencils )
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
      double operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags ) const
      {
        double elapsedTime = - MPI_Wtime();

        reconstructions.clear();
        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          //if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) && !flags.isActive( entity ) )
            //continue;

          applyLocal( entity, flags, color, reconstructions[ entity ] );
        }
        elapsedTime += MPI_Wtime();

        auto exchange = typename ReconstructionSet::Exchange ( reconstructions );
        color.gridView().communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );

        return elapsedTime;
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
        Coordinate normal( 0.0 );

        const auto geometry = entity.geometry();

        const Coordinate center = geometry.center();
        const auto colorEn = color[ entity ];

        // second order
        if ( stencil( entity ).size() >= 8 )
        {
          MatrixD AtA( 0.0 );
          VectorD Atb( 0.0 );
          for ( const auto &neighbor : stencil( entity ) )
          {
            Coordinate d = neighbor.geometry().center() - center;

            VectorD v; // v = ( dx,  dy, dxx, dyy, dxy )
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

          VectorD solution;
          AtA.solve( solution, Atb );

          for ( std::size_t i = 0; i < dim; ++i )
            normal[ i ] = solution[ i ];
        }
        else // first order
        {
          Matrix AtA( 0.0 );
          Vector Atb( 0.0 );
          for ( const auto &neighbor : stencil( entity ) )
          {
            Coordinate d = neighbor.geometry().center() - center;
            const ctype weight = 1.0 / d.two_norm2();
            d *= weight;
            AtA += outerProduct( d, d );
            Atb.axpy( weight * ( color[ neighbor ] - colorEn ), d );
          }
          AtA.solve( normal, Atb );
        }

        if( normal.two_norm2() < std::numeric_limits< ctype >::epsilon() )
          return;

        normalize( normal );

        if ( flags.isMixed( entity ) || flags.isFullAndMixed( entity ) )
        {
          auto polytope = makePolytope( geometry );
          reconstruction = locateHalfSpace( polytope, normal, colorEn );
        }
      }

      Matrix outerProduct ( const Vector &a, const Vector &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < dim; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      MatrixD outerProduct ( const VectorD &a, const VectorD &b ) const
      {
        MatrixD m( 0.0 );
        for ( std::size_t i = 0; i < derivatives; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const StencilSet &stencils_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
