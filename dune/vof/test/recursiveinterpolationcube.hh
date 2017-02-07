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
          uh[ entity ] = recursiveEvaluation( LocalCoordinate( 0 ), 1.0, localfunction, depth_ );
        }
      }

    private:
      template< class LocalFunction >
      static double recursiveEvaluation( const LocalCoordinate &origin, ctype h, const LocalFunction &lf, std::size_t level )
      {
        const int corners = 1 << dim;

        int sum = 0;
        for ( int i = 0; i < corners; ++i )
        {
          LocalCoordinate corner ( origin );
          for ( int j = 0; j < dim; ++j )
            corner[ j ] += h * static_cast< ctype >( ( i >> j ) & 1 );
          sum += lf( corner );
        }

        if ( level == 0 )
          return volume( h ) * ( static_cast< double > ( sum ) / static_cast< double >( corners ) );

        if ( sum == corners )
          return volume( h );
        if ( sum == 0 )
          return 0.0;

        h /= 2;
        double value = 0.0;
        for ( int i = 0; i < corners; ++i )
        {
          LocalCoordinate childOrigin( origin );
          for ( int j = 0; j < dim; ++j )
            childOrigin[ j ] += h * static_cast< ctype >( ( i >> j ) & 1 );
          value += recursiveEvaluation( childOrigin, h, lf, level - 1 );
        }

        return value;
      }

      static ctype volume ( double h )
      {
        double volume = 1.0;
        for ( std::size_t i = 0; i < dim; ++i )
          volume *= h;
        return volume;
      }

      const GridView &gridView_;
      const std::size_t depth_;
    };

  } // namespace VoF

}  // namespace Dune

#endif

