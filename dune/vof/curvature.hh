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
    template< class GV, class ST, class RS, class FL >
    struct Curvature
    {
      using GridView = GV;
      using Stencils = ST;
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

      void operator() ( const ReconstructionSet &reconstructions, const Flags &flags )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature_[ index( entity ) ] = 0.0;

          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, reconstructions, flags );
        }
      }

      double operator[] ( const Entity& entity ) const
      {
        return curvature_[ index( entity ) ];
      }

    private:
      void applyLocal ( const Entity &entity, const ReconstructionSet &reconstructions, const Flags &flags )
      {
        double sumWeights = 0.0;
        double divN ( 0.0 );

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();

        sumWeights = 0.0;
        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffN = reconstructions[ entity ].innerNormal() - reconstructions[ neighbor ].innerNormal();
          Coordinate diffC = centroidEn - centroidNb;

          double weight = diffC.two_norm2();

          divN += weight * ( diffN * diffC ) / diffC.two_norm2();
          sumWeights += weight;
        }

        curvature_[ index( entity ) ] = - divN / sumWeights;

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

      const GridView &gridView () const { return gridView_; }
      const IndexSet& indexSet () const { return gridView().indexSet(); }
      Index index ( const Entity &entity ) const { return indexSet().index( entity ); }

      GridView gridView_;
      const Stencils &stencils_;
      std::vector< double > curvature_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURE_HH
