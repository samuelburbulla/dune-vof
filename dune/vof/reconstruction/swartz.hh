#ifndef DUNE_VOF_RECONSTRUCTION_SWARTZ_HH
#define DUNE_VOF_RECONSTRUCTION_SWARTZ_HH

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

    // SwartzReconstruction
    // --------------------


    /**
     * \ingroup   Reconstruction
     * \brief     modified Swartz reconstruction operator
     * \details   Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 15ff
     *
     * \tparam  GV  grid view
     * \tparam  StS stencils type
     */
    template< class GV, class StS >
    struct SwartzReconstruction
    {
      using GridView = GV;
      using StencilSet = StS;

    private:
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;
      using ctype = typename Coordinate::value_type;


    public:
      explicit SwartzReconstruction ( const StencilSet &stencils, const std::size_t maxIterations = 50 )
       : stencils_( stencils ), maxIterations_( maxIterations )
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
        initializer_( color, reconstructions, flags );

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
        Coordinate normalOld = normal;

        const auto geometry = entity.geometry();
        Coordinate center = geometry.center();
        const double colorEn = clamp( color[ entity ], 0.0, 1.0 );

        auto fInverse = [ &normalOld ] ( const auto &geometry, const auto color )
        {
          const auto polytope = makePolytope( geometry );
          const auto hs = locateHalfSpace( polytope, normalOld, color );
          return hs.levelSet( geometry.center() );
        };

        auto computeNormal = [ normalOld ] ( const Coordinate &x1, const Coordinate &x2, const double a1, const double a2, Coordinate &normal ) -> int
        {
          using Matrix = FieldMatrix< ctype, dim, dim >;
          using Vector = FieldVector< ctype, dim >;

          const Matrix X { x1 - Coordinate( 0.5 ), x2 - Coordinate( 0.5 ) };

          Vector a { a1, a2 };
          Vector one ( -1.0 );

          Vector b, q;
          X.solve( b, a );
          X.solve( q, one );

          Vector r1, r2;
          X.mv( b, r1 );
          X.mv( q, r2 );

          double a_ = q * q;
          double b_ = 2.0 * ( q * b );
          double c_ = b * b - 1.0;
          double diskr = b_ * b_ - 4.0 * a_ * c_;
          if ( diskr < 1e-14 )
            return 0;
          else if ( std::abs( diskr ) < 1e-14 )
          {
            double s = - b_ / ( 2.0 * a_ );
            normal += b;
            normal.axpy( s, q );
            return 1;
          }
          else
          {
            double s1 = ( - b_ + std::sqrt( diskr ) ) / ( 2.0 * a_ );
            double s2 = ( - b_ - std::sqrt( diskr ) ) / ( 2.0 * a_ );

            Coordinate guess = b;
            guess.axpy( s1, q );

            if ( normalOld * guess > 0 )
            {
              normal += guess;
              return 1;
            }
            else
            {
              guess.axpy( s2 - s1, q );
              normal += guess;
              return 1;
            }
          }
        };


        HalfSpace< Coordinate > hsEn;
        do
        {
          int n = 0;
          normalOld = normal;
          normal = Coordinate ( 0.0 );

          double sigmaEn = fInverse( geometry, colorEn );

          for ( const auto &neighbor : stencil( entity ) )
          {
            if ( !flags.isMixed( neighbor ) )
              continue;

            double colorNb = clamp( color[ neighbor ], 0.0, 1.0 );
            const auto geometryNb = neighbor.geometry();
            auto centerNb = geometryNb.center();

            double sigmaNb = fInverse( geometryNb, colorNb );

            n += computeNormal( center, centerNb, sigmaEn, sigmaNb, normal );
          }
          if( n == 0 )
            return;
          normalize( normal );
    break; // makes second order too!
        }
        while ( std::acos( normal * normalOld ) > 1e-6 );

        reconstruction = locateHalfSpace( makePolytope( geometry ), normal, colorEn );
      }

      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const StencilSet &stencils_;
      const std::size_t maxIterations_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_SWARTZ_HH
