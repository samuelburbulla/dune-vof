#ifndef DUNE_VOF_PARABOLAINTERPOLATIONCURVATURE_HH
#define DUNE_VOF_PARABOLAINTERPOLATIONCURVATURE_HH

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

    // Calculate the curvature of the inferface in each mixed cell by an interpolation of a parabola.

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
    struct ParabolaInterpolationCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, 2, 2 >;
      using Vector = FieldVector< ctype, 2 >;

    public:
      explicit ParabolaInterpolationCurvature ( GridView gridView, const Stencils &stencils )
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

        curvature.communicate();
      }

    private:
      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature )
      {
        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();

        /*
        Vector p ( { 0, 0 } );
        const ctype weight = 1.0;
        p *= weight;
        AtA += outerProduct( p, p );
        Atb.axpy( 0.0, p );
        */

        // with normal
        Coordinate n = reconstructions[ entity ].innerNormal();
        Vector p2 ( { 2.0 * 0.0, 1.0 } );
        const ctype weight = interfaceEn.volume();
        AtA.axpy( weight, outerProduct( p2, p2 ) );
        Atb.axpy( weight * ( - ( n * generalizedCrossProduct( normalEn ) ) / ( n * normalEn ) ), p2 );

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) )
           continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffC = centroidNb - centroidEn;

          double dx = generalizedCrossProduct( normalEn ) * diffC;
          double dy = normalEn * diffC;

          /*
          // fit centroids
          Vector p ( { dx * dx, dx } );
          const ctype weight = 1.0 / diffC.two_norm();
          AtA.axpy( weight, outerProduct( p, p ) );
          Atb.axpy( weight * dy, p );
          */

          // with normal
          Coordinate n = reconstructions[ neighbor ].innerNormal();
          Vector p2 ( { 2.0 * dx, 1.0 } );
          const ctype weight = interfaceNb.volume();;
          AtA.axpy( weight, outerProduct( p2, p2 ) );
          Atb.axpy( weight * ( - ( n * generalizedCrossProduct( normalEn ) ) / ( n * normalEn ) ), p2 );

        }
        Vector abc;
        AtA.solve( abc, Atb );

        curvature[ entity ] = - 2.0 * abc[0] / std::pow( 1.0 + abc[1] * abc[1], 3.0 / 2.0 );
      }

      const GridView &gridView () const { return gridView_; }
      const auto &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      GridView gridView_;
      const Stencils &stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_PARABOLAINTERPOLATIONCURVATURE_HH
