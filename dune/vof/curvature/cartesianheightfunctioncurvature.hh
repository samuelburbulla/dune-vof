#ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
#define DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-vof includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/stencil/heightfunctionstencil.hh>


namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell by the height function method.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV  grid view
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class VNST, class DF, class RS, class FL >
    struct CartesianHeightFunctionCurvature
    {
      using GridView = GV;
      using Stencils = Dune::VoF::HeightFunctionStencils< GridView >;
      using VertexNeighborStencils = VNST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      using Heights = Dune::FieldVector< double, Stencils::Stencil::noc >;

    public:
      explicit CartesianHeightFunctionCurvature ( GridView gridView, const VertexNeighborStencils &vertexNeighborStencils )
       : gridView_( gridView ),
         stencils_( gridView ),
         vertexNeighborStencils_( vertexNeighborStencils ),
         indexSet_( gridView.indexSet() ),
         satisfiesConstraint_( gridView )
      {}

      template< class CurvatureSet >
      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature, bool communicate = false )
      {
        for ( const auto& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          curvature[ entity ] = 0.0;
          satisfiesConstraint( entity ) = 0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, curvature );
        }

        curvature.communicate();
        satisfiesConstraint_.communicate();

        for ( const auto& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          if ( satisfiesConstraint( entity ) )
            continue;

          averageCurvature( entity, curvature );
        }

        if ( communicate )
          curvature.communicate();
      }

    private:
      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, CurvatureSet &curvature )
      {
        const Coordinate &normal = reconstructions[ entity ].innerNormal();
        auto stencil = heightFunctionStencil( normal, entity );

        Heights heights ( 0.0 );

        for( std::size_t i = 0; i < stencil.columns(); ++i )
          for( int t = stencil.tdown(); t <= stencil.tup(); ++t )
          {
            if ( !stencil.valid( i, t ) )
              continue;

            double u = uh[ stencil( i, t ) ];

            // local monotonic variation
            if ( t < 0 )
            {
              if ( u < uh[ stencil( i, t+1 ) ] - 1e-8 )
                u = 1.0;
            }
            else if ( t > 0 )
            {
              if ( u > uh[ stencil( i, t-1 ) ] + 1e-8 )
                u = 0.0;
            }

            heights[ i ] += u;
          }

        // Constraint
        double uMid = heights[ ( heights.size() - 1 ) / 2 ];
        int effTdown = stencil.effectiveTdown();
        if ( uMid < effTdown || uMid > effTdown + 1 )
          return;

        for( std::size_t i = 0; i < decltype( stencil )::noc; ++i )
          if ( heights[ i ] == 0.0 )
            return;

        satisfiesConstraint( entity ) = 1;

        curvature[ entity ] = kappa( heights, getOrientation( normal ), stencils_.deltaX() );
      }


      template< class Coordinate >
      auto kappa ( const Heights &heights, const Coordinate &orientation, const double dx ) const
      -> typename std::enable_if< Coordinate::dimension == 2, double >::type
      {
        double Hx = ( heights[ 2 ] - heights[ 0 ] ) / 2.0;
        double Hxx = ( heights[ 2 ] - 2 * heights[ 1 ] + heights[ 0 ] ) / dx;

        return - Hxx / std::pow( 1.0 + Hx * Hx, 3.0 / 2.0 );
      }

      template< class Coordinate >
      auto kappa ( const Heights &heights, const Coordinate &orientation, const double dx ) const
      -> typename std::enable_if< Coordinate::dimension == 3, double >::type
      {
        double Hx = ( heights[ 5 ] - heights[ 3 ] ) / 2.0;
        double Hy = ( heights[ 7 ] - heights[ 1 ] ) / 2.0;
        double Hxx = ( heights[ 5 ] - 2.0 * heights[ 4 ] + heights[ 3 ] ) / dx;
        double Hyy = ( heights[ 7 ] - 2.0 * heights[ 4 ] + heights[ 1 ] ) / dx;
        double Hxy = ( heights[ 8 ] - heights[ 2 ] - heights[ 6 ] + heights[ 0 ] ) / ( 4.0 * dx );

        return - ( Hxx + Hyy + Hxx * Hy * Hy + Hyy * Hx * Hx - 2.0 * Hxy * Hx * Hy ) / ( std::pow( 1.0 + Hx * Hx + Hy * Hy, 3.0 / 2.0 ) );
      }

      static inline Coordinate getOrientation( const Coordinate &normal )
      {
        std::size_t dir = 0;
        double max = std::numeric_limits< double >::min();
        for ( std::size_t i = 0; i < Coordinate::dimension; ++i )
        if ( std::abs( normal[ i ] ) > max )
        {
          dir = i;
          max = std::abs( normal[ i ] );
        }
        Coordinate orientation;
        orientation[ dir ] = ( normal[ dir ] > 0 ) ? -1.0 : 1.0;

        return orientation;
      }

      template< class CurvatureSet >
      void averageCurvature( const Entity &entity, CurvatureSet &curvature ) const
      {
        int n = 0;
        for( const auto& neighbor : vertexNeighborStencil( entity ) )
        {
          if ( satisfiesConstraint( neighbor ) )
          {
            curvature[ entity ] += curvature[ neighbor ];
            n++;
          }
        }
        if ( n > 0 )
          curvature[ entity ] /= n;
      }

      const GridView &gridView () const { return gridView_; }

      const IndexSet &indexSet () const { return indexSet_; }

      const auto index ( const Entity &entity ) const { return indexSet_.index( entity ); }

      const auto &heightFunctionStencil ( const Coordinate &normal, const Entity &entity ) const { return stencils_( normal, entity ); }

      const auto &vertexNeighborStencil ( const Entity &entity ) const { return vertexNeighborStencils_[ entity ]; }

      const std::size_t satisfiesConstraint ( const Entity &entity ) const { return satisfiesConstraint_[ entity ]; }
      std::size_t &satisfiesConstraint ( const Entity &entity ) { return satisfiesConstraint_[ entity ]; }

      GridView gridView_;
      Stencils stencils_;
      const VertexNeighborStencils &vertexNeighborStencils_;
      IndexSet indexSet_;
      Dune::VoF::DataSet< GridView, std::size_t > satisfiesConstraint_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
