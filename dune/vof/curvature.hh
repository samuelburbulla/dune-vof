#ifndef DUNE_VOF_CURVATURE_HH
#define DUNE_VOF_CURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell.

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
    struct Curvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;
      using ctype = typename Coordinate::value_type;
      static constexpr int dim = Coordinate::dimension;
      using Matrix = FieldMatrix< ctype, dim, dim >;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );

    public:
      explicit Curvature ( const GridView &gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils ), curvature_( indexSet().size( 0 ) )
      {}

      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature_[ index( entity ) ] = 0.0;

          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, flags );
        }
      }

      double operator[] ( const Entity& entity ) const
      {
        return curvature_[ index( entity ) ];
      }

      double operator[] ( const int index ) const
      {
        return curvature_[ index ];
      }

      std::size_t size () const { return curvature_.size(); }

      const GridView &gridView () const { return gridView_; }

    private:
      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags )
      {
        /*
        double uxx = 0.0, uxy = 0.0, uyy = 0.0;

        Coordinate cEn = entity.geometry().center();
        Coordinate ex ( { 1.0, 0.0 } );
        Coordinate ey ( { 1.0, 0.0 } );
        Coordinate exy = ex + ey;
        exy /= exy.two_norm();


        for( const auto& neighbor : stencils_[ entity ] )
        {
          Coordinate cNb = neighbor.geometry().center();
          Coordinate dX = cNb - cEn;

          uxx += ( uh[ neighbor ] - uh[ entity ] ) * std::abs( ex * dX ) / dX.two_norm2();
          uyy += ( uh[ neighbor ] - uh[ entity ] ) * std::abs( ey * dX ) / dX.two_norm2();
          uxy += ( uh[ neighbor ] - uh[ entity ] ) * std::abs( exy * dX ) / dX.two_norm2();
        }

        double ux = reconstructions[ entity ].innerNormal()[ 0 ];
        double uy = reconstructions[ entity ].innerNormal()[ 1 ];

        double tmp = std::sqrt( 1.0 + ux * ux + uy * uy );
        curvature_[ index( entity ) ] = uxx; //( uxx + uyy + uxx * uy * uy + uyy * ux * ux - 2.0 * uxy * ux * uy ) / ( tmp * tmp * tmp );
        */
        int n = 0;
        double h = 0.0;
        double divN ( 0.0 );

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffC = centroidEn - centroidNb;

          divN += ( normalEn * diffC ) / diffC.two_norm();
          h += std::abs( diffC * generalizedCrossProduct( normalEn ) );
          n++;
        }
        h /= n;
        curvature_[ index( entity ) ] = - ( divN / h );

        /*
        // Interpolate with circle through points
        Matrix AtA( 0.0 );
        Coordinate Atb( 0.0 );
        std::vector< Coordinate > points;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          points.push_back( centroidNb );
        }

        for ( const auto &point : points )
        {
          Coordinate p = point;
          p -= centroidEn;
          const ctype weight = p.two_norm2();
          p *= 2.0 * weight;
          AtA += outerProduct( p, p );
          Atb.axpy( weight * ( point.two_norm2() - centroidEn.two_norm2() ), p );
        }
        Coordinate centerOfCircle;
        AtA.solve( centerOfCircle, Atb );

        double radius = ( centerOfCircle - centroidEn ).two_norm();
        double sign = - 1 + 2 * ( 0 < ( centerOfCircle - centroidEn ) * reconstructions[ entity ].innerNormal() );
        curvature_[ index( entity ) ] = sign / radius;
        */

      }

      auto interface( const Entity &entity, const ReconstructionSet &reconstructions ) const
      {
        auto polygon = makePolytope( entity.geometry() );
        auto it = intersect( polygon, reconstructions[ entity ].boundary() );
        auto interface = static_cast< typename decltype( it )::Result > ( it );
        return interface;
      }

      Matrix outerProduct ( const Coordinate &a, const Coordinate &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < dim; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      const IndexSet& indexSet () const { return gridView().indexSet(); }
      Index index ( const Entity &entity ) const { return indexSet().index( entity ); }

      GridView gridView_;
      const Stencils &stencils_;
      std::vector< double > curvature_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURE_HH
