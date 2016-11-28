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
        // Height function
        Coordinate normal = reconstructions[ entity ].innerNormal();
        Coordinate center = entity.geometry().center();

        Coordinate direction ( 0.0 );
        Coordinate directionOrth ( 0.0 );
        if ( std::abs( normal[ 0 ] ) > std::abs( normal[ 1 ] ) )
        {
         direction[ 0 ] = 1.0;
         directionOrth[ 1 ] = 1.0;
        }
        else
        {
         direction[ 1 ] = 1.0;
         directionOrth[ 0 ] = 1.0;
        }

        double left = 0.0, right = 0.0;
        double middle = uh[ entity ];
        double dx = std::numeric_limits< double >::max();
        for( const auto& neighbor : stencils_[ entity ] )
        {
          Coordinate d = neighbor.geometry().center() - center;

          if ( std::abs( d * directionOrth ) < std::numeric_limits< double >::epsilon() )
            middle += uh[ neighbor ];
          else if ( d * directionOrth > 0 )
            right += uh[ neighbor ];
          else if ( d * directionOrth < 0 )
            left += uh[ neighbor ];


          dx = std::min( dx, d.two_norm() );
        }
        right *= dx;
        left *= dx;
        middle *= dx;

        if ( dx < middle && middle < 2 * dx )
        {
          double Hx = ( right - left ) / ( 2 * dx );
          double Hxx = ( right - 2 * middle + left ) / ( dx * dx );

          curvature_[ index( entity ) ] = - Hxx / std::pow( 1.0 + Hx * Hx, 2.0 / 3.0 );
        }
        */
        /*
        // Least squares for gradients
        Coordinate center = entity.geometry().center();

        for ( std::size_t k = 0; k < dim; ++k )
        {
          Matrix AtA( 0.0 );
          Coordinate Atb( 0.0 );

          for( const auto& neighbor : stencils_[ entity ] )
          {
            Coordinate d = neighbor.geometry().center() - center;
            const ctype weight = 1.0 / d.two_norm2();
            d *= weight;
            AtA += outerProduct( d, d );
            Atb.axpy( weight * ( reconstructions[ neighbor ].innerNormal()[ k ] - reconstructions[ entity ].innerNormal()[ k ] ), d );
          }

          Coordinate dNk;
          AtA.solve( dNk, Atb );

          curvature_[ index( entity ) ] -= dNk[ k ];
        }
        */
        // Finite differences
        int n = 0;
        double h = 0.0;
        double divN ( 0.0 );

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffC = centroidEn - centroidNb;

          double dx = std::abs( generalizedCrossProduct( normalEn ) * diffC );
          divN += ( normalEn * diffC ) / ( dx * dx );
          n++;
        }
        curvature_[ index( entity ) ] = - divN * 2.0 / n;
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
