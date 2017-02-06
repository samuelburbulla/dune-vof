#ifndef DUNE_VOF_RECURSIVEINTERPOLATION_HH
#define DUNE_VOF_RECURSIVEINTERPOLATION_HH

// dune-geometry includes
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/virtualrefinement.hh>

namespace Dune
{
  namespace VoF
  {
    // Recursive interpolation of arbitrary initial data
    // =================================================

    template< class GridView >
    struct RecursiveInterpolation
    {
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using ctype = typename GridView::ctype;
      static constexpr int dim = GridView::dimension;
      using LocalCoordinate = FieldVector< ctype, dim >;
      using Refinement = VirtualRefinement< dim, ctype >;

    public:
      RecursiveInterpolation( const GridView &gridView, std::size_t depth ) : gridView_( gridView ), depth_( depth ) {}

      template < class Function, class DiscreteFunction >
      void operator() ( const Function &f, DiscreteFunction &uh )
      {
        for ( const auto &entity : elements( gridView_ ) )
        {
          const auto geometry = entity.geometry();
          auto localfunction = [ &f, geometry ]( const LocalCoordinate &x ) -> int { return (int) ( f( geometry.global( x ) ) + 0.3 ); };
          const auto refElementGeometry = ReferenceElements< ctype, dim >::general( geometry.type() ).template geometry< 0 >( 0 );
          uh[ entity ] = recursiveEvaluation( refElementGeometry, localfunction, depth_ );
        }
      }

    private:
      template< class Embedding, class LocalFunction >
      static double recursiveEvaluation( const Embedding &embedding, const LocalFunction &lf, std::size_t level )
      {
        int sum = 0;
        for ( int i = 0; i < embedding.corners(); ++i )
          sum += lf( embedding.corner( i ) );

        if ( level == 0 )
          return embedding.volume() * ( static_cast<double> ( sum ) / static_cast< double >( embedding.corners() ) );

        if ( sum == embedding.corners() )
          return embedding.volume();
        if ( sum == 0 )
          return 0.0;

        std::vector< LocalCoordinate > subCoordinates;

        Refinement &refinement = buildRefinement< dim, ctype >( embedding.type(), embedding.type() );
        for( auto sit = refinement.vBegin( 1 ), send = refinement.vEnd( 1 ); sit != send; ++sit )
          subCoordinates.push_back( embedding.global( sit.coords() ) );

        double value = 0.0;
        for( auto sit = refinement.eBegin( 1 ), send = refinement.eEnd( 1 ); sit != send; ++sit )
        {
          std::vector< LocalCoordinate > subCorners;
          for( int i : sit.vertexIndices() )
            subCorners.push_back( subCoordinates[ i ] );

          MultiLinearGeometry< ctype, dim, dim > subEmbedding ( embedding.type(), subCorners );

          value += recursiveEvaluation( subEmbedding, lf, level - 1 );
        }

        return value;
      }

      const GridView &gridView_;
      const std::size_t depth_;
    };

  } // namespace VoF

}  // namespace Dune

#endif

