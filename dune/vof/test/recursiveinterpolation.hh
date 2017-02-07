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
      using Matrix = FieldMatrix< ctype, dim+1, dim+1 >;
      using Vector = FieldVector< ctype, dim+1 >;

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


      template< class Embedding, class LocalFunction >
      static double recursiveEvaluation( const Embedding &embedding, const LocalFunction &lf, std::size_t level, bool allLevels = false )
      {
        if ( !allLevels || level == 0 )
        {
          int sum = 0;
          for ( int i = 0; i < embedding.corners(); ++i )
            sum += lf( embedding.corner( i ) );

          if ( sum == embedding.corners() )
            return embedding.volume();
          if ( sum == 0 )
            return 0.0;

          // Fit halfspace to the phase change point in all edges
          if ( level == 0 )
          {
            Matrix AtA( 0.0 );
            Vector Atb( 0.0 );

            LocalCoordinate someInnerPoint;
            std::vector< LocalCoordinate > points;

            const auto &ref = ReferenceElements< ctype, dim >::general( embedding.type() );

            for ( int s = 0; s < ref.size( dim-1 ); ++s )
            {
              auto edge = ref.template geometry< dim-1 >( s );
              LocalCoordinate p0 = embedding.global( edge.corner( 0 ) );
              LocalCoordinate p1 = embedding.global( edge.corner( 1 ) );

              if ( !( lf( p0 ) ^ lf( p1 ) ) )
                continue;

              if ( lf( p0 ) == 1 )
                someInnerPoint = p0;
              if ( lf( p1 ) == 1 )
                someInnerPoint = p1;

              auto convexCombination = [ p0, p1, lf ]( ctype x ) -> int
              {
                LocalCoordinate pos;
                pos.axpy( 1.0 - x, p0 );
                pos.axpy( x, p1 );
                return lf( pos );
              };
              ctype result = bisection( convexCombination );
              LocalCoordinate intersection;
              intersection.axpy( 1.0 - result, p0 );
              intersection.axpy( result, p1 );

              points.push_back( intersection );
            }

            HalfSpace< LocalCoordinate > hs;

            if ( points.size() < dim )
              return 0.0;
            else if ( points.size() == dim )
              hs = HalfSpace< LocalCoordinate >( points );
            /*
            else
            {
              for ( const auto point : points )
              {
                Vector d;
                for ( std::size_t i = 0; i < dim; ++i )
                  d[ i ] = point[ i ];
                d[ dim ] = 1.0;

                const ctype weight = 1.0;
                d *= weight;
                AtA += outerProduct( d, d );
                Atb.axpy( 0.0, d );
              }
              Vector solution;
              AtA.solve( solution, Atb );

              LocalCoordinate normal;
              for ( std::size_t i = 0; i < dim; ++i )
                normal[ i ] = solution[ i ];

              hs = HalfSpace< LocalCoordinate >( normal, solution[ dim ] );
            }
            */

            if ( hs.levelSet( someInnerPoint ) < 0 )
              hs = hs.flip();

            auto polytope = makePolytope( embedding );
            auto part = intersect( polytope, hs, eager );
            return part.volume();
          }
        }
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

          value += recursiveEvaluation( subEmbedding, lf, level-1, allLevels );
        }

        return value;
      }

    private:
      template< class Function >
      static ctype bisection( const Function &f )
      {
        assert ( f( 0.0 ) ^ f( 1.0 ) );
        return bisectionRecursion( f, 0.0, 1.0, 45 );
      }

      template< class Function >
      static ctype bisectionRecursion( const Function &f, ctype a, ctype b, int level )
      {
        ctype m = 0.5 * ( a + b );

        if ( level == 0 )
          return m;

        if ( f( a ) ^ f( m ) )
          return bisectionRecursion( f, a, m, level-1 );
        else
          return bisectionRecursion( f, m, b, level-1 );
      }

      const GridView &gridView_;
      const std::size_t depth_;
    };

  } // namespace VoF

}  // namespace Dune

#endif

