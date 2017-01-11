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
        double tol = std::numeric_limits< double >::epsilon();

        Coordinate normal = reconstructions[ entity ].innerNormal();
        Coordinate center = entity.geometry().center();

        std::size_t directionAxis = ( std::abs( normal[ 0 ] ) > std::abs( normal[ 1 ] ) ) ? 0 : 1;
        Coordinate direction ( 0.0 );
        direction[ directionAxis ] = 1.0;

        double left = 0.0, right = 0.0;
        double middle = uh[ entity ];
        double dx = std::numeric_limits< double >::max();
        for( const auto& neighbor : stencil( directionAxis, entity ) )
        {
          Coordinate d = center - neighbor.geometry().center();

          Coordinate delta = d;
          delta.axpy( - ( d * direction ), direction );

          double distance = delta[ 1 - directionAxis ];

          if ( std::abs( distance ) < tol )
            middle += uh[ neighbor ];
          else if ( distance > 0 )
            left += uh[ neighbor ];
          else if ( distance < 0 )
            right += uh[ neighbor ];

          if ( 0 < std::abs( distance ) && std::abs( distance ) < dx )
            dx = std::abs( distance );
        }

        right *= dx;
        left *= dx;
        middle *= dx;

        if ( 3.0 * dx < middle && middle < 4.0 * dx )
        {
          satisfiesConstraint( entity ) = 1;

          if ( right < tol || left < tol )
            return;

          double Hx = ( right - left ) / ( 2 * dx );
          double Hxx = ( right - 2 * middle + left ) / ( dx * dx );

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

      const auto &stencil ( const std::size_t d, const Entity &entity ) const { return stencils_( d, entity ); }

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
