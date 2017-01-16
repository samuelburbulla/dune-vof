#ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
#define DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH

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

    // Calculate the curvature of the inferface in each mixed cell by the height function method.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV  grid view
     * \tparam  ST  stencils
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class ST, class VNST, class DF, class RS, class FL >
    struct CartesianHeightFunctionCurvature
    {
      using GridView = GV;
      using Stencils = ST;
      using VertexNeighborStencils = VNST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit CartesianHeightFunctionCurvature ( GridView gridView, const Stencils &stencils, const VertexNeighborStencils &vertexNeighborStencils )
       : gridView_( gridView ),
         stencils_( stencils ),
         vertexNeighborStencils_( vertexNeighborStencils ),
         indexSet_( gridView.indexSet() ),
         satisfiesConstraint_( indexSet_.size( 0 ), 0 )
      {}

      template< class CurvatureSet >
      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags, CurvatureSet &curvature )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature[ entity ] = 0.0;
          satisfiesConstraint( entity ) = 0;

          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, curvature );
        }

        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          if ( satisfiesConstraint( entity ) )
            continue;

          averageCurvature( entity, curvature );
        }

      }

    private:

      template< class CurvatureSet >
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, CurvatureSet &curvature )
      {
        const Coordinate &normal = reconstructions[ entity ].innerNormal();
        auto stencil = heightFunctionStencil( normal, entity );

        Dune::FieldVector< double, decltype( stencil )::noc > heights ( 0.0 );

        for( int i = stencil.cMin(); i < stencil.cMax(); ++i )
          for( int t = -stencil.tdown(); t <= stencil.tup(); ++t )
          {
            double u = uh[ stencil( i, t ) ];

            // local monotonic variation
            if ( t < 0 )
            {
              if ( u < uh[ stencil( i, t+1 ) ] )
                u = 1.0;
            }
            else if ( t > 0 )
            {
              if ( u > uh[ stencil( i, t-1 ) ] )
                u = 0.0;
            }

            heights[ i - stencil.cMin() ] += u;
          }

        if ( stencil.tdown() < heights[ 1 ] && heights[ 1 ] < stencil.tdown() + 1 )
        {
          satisfiesConstraint( entity ) = 1;

          for( std::size_t i = 0; i < stencil.columns(); ++i )
            if ( heights[ i ] == 0.0 )
              return;

          double Hx = ( heights[ 2 ] - heights[ 0 ] ) / 2.0;
          double Hxx = ( heights[ 2 ] - 2 * heights[ 1 ] + heights[ 0 ] ) / stencils_.deltaX();

          curvature[ entity ] = - Hxx / std::pow( 1.0 + Hx * Hx, 3.0 / 2.0 );
        }
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

      const std::size_t satisfiesConstraint ( const Entity &entity ) const { return satisfiesConstraint_[ index( entity ) ]; }
      std::size_t &satisfiesConstraint ( const Entity &entity ) { return satisfiesConstraint_[ index( entity ) ]; }

      GridView gridView_;
      const Stencils &stencils_;
      const VertexNeighborStencils &vertexNeighborStencils_;
      IndexSet indexSet_;
      std::vector< std::size_t > satisfiesConstraint_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CARTESIANHEIGHTFUNCTIONCURVATURE_HH
