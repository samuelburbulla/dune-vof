#ifndef DUNE_VOF_NORMALDIVERGENCECURVATURE_HH
#define DUNE_VOF_NORMALDIVERGENCECURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-grid includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell by the divergence of the normal vector field.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV  grid view
     * \tparam  ST  stencils
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class ST, class DF, class RS, class FL >
    struct NormalDivergenceCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;
      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, dim, dim >;
      using Vector = FieldVector< ctype, dim >;

    public:
      explicit NormalDivergenceCurvature ( const GridView &gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils )
      {}

      template< class CurvatureSet >
      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature[ entity ] = 0.0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, flags, curvature );
        }

        CurvatureSet tmpCurvature( curvature );
        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applySmoothing1st( entity, uh, tmpCurvature, curvature );
        }

        CurvatureSet tmpCurvature2( curvature );
        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applySmoothing2nd( entity, uh, reconstructions, tmpCurvature2, curvature );
        }

      }

    private:
      double smoothingWeight( const double color )
      {
        return std::pow( 1.0 - 2.0 * std::abs( 0.5 - color ), 8.0 );
      }

      template< class CurvatureSet >
      void applySmoothing1st ( const Entity &entity, const DiscreteFunction &uh, const CurvatureSet &tmpCurvature, CurvatureSet &curvature )
      {
        double kappa = tmpCurvature[ entity ] * smoothingWeight( uh[ entity ] );
        double sumWeights = smoothingWeight( uh[ entity ] );

        for ( const auto &neighbor : stencil( entity ) )
        {
          kappa += tmpCurvature[ neighbor ] * smoothingWeight( uh[ neighbor ] );
          sumWeights += smoothingWeight( uh[ neighbor ] );
        }

        kappa /= sumWeights;
        curvature[ entity ] = kappa;
      }

      double additionalSmoothingWeight( const Coordinate &normal, const Coordinate &delta )
      {
        return std::pow( std::abs( normal * delta ), 8.0 );
      }

      template< class CurvatureSet >
      void applySmoothing2nd ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const CurvatureSet &tmpCurvature, CurvatureSet &curvature )
      {
        Coordinate normal = reconstructions[ entity ].innerNormal();

        double kappa = tmpCurvature[ entity ] * smoothingWeight( uh[ entity ] );
        double sumWeights = smoothingWeight( uh[ entity ] );

        for ( const auto &neighbor : stencil( entity ) )
        {
          Coordinate delta = neighbor.geometry().center() - entity.geometry().center();
          delta /= delta.two_norm();

          kappa += tmpCurvature[ neighbor ] * smoothingWeight( uh[ neighbor ] ) * additionalSmoothingWeight( normal, delta );
          sumWeights += smoothingWeight( uh[ neighbor ] ) * additionalSmoothingWeight( normal, delta );
        }

        kappa /= sumWeights;
        curvature[ entity ] = kappa;
      }

      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature )
      {
        Coordinate center = entity.geometry().center();

        Matrix nablaN( 0.0 );

        for ( std::size_t k = 0; k < dim; ++k )
        {
          Matrix AtA( 0.0 );
          Coordinate Atb( 0.0 );

          for( const auto &neighbor : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor ) )
              continue;

            Coordinate centerNb = neighbor.geometry().center();
            Coordinate d = centerNb - center;
            const ctype weight = 1.0 / d.two_norm2();
            d *= weight;
            AtA += outerProduct( d, d );
            Atb.axpy( weight * ( reconstructions[ neighbor ].innerNormal()[ k ] - reconstructions[ entity ].innerNormal()[ k ] ), d );
          }

          Coordinate dNk ( 0.0 );
          AtA.solve( dNk, Atb );
          curvature[ entity ] += dNk[ k ];
        }
      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const GridView &gridView_;
      const Stencils &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_NORMALDIVERGENCECURVATURE_HH
