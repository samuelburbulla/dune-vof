#ifndef DUNE_VOF_SECONGYOUNGSNORMALGUESSING_HH
#define DUNE_VOF_SECONGYOUNGSNORMALGUESSING_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>


namespace Dune
{
  namespace VoF
  {

    template< class Grid, class DomainVector >
    class SecondYoungsNormalGuessing
    {

    public:
      template< class C, class D >
      SecondYoungsNormalGuessing ( const Grid &grid, C &concentration, std::vector< bool > &cellIsMixed, const D &domain, const double eps )
      {
        auto gridView = grid.leafGridView();

        normals.resize( gridView.indexSet().size( 0 ) );

        for( auto&& entity : elements( gridView ) )
        {
          int entityIndex = gridView.indexSet().index( entity );

          if( cellIsMixed[ entityIndex ] && concentration[ entityIndex ] <= 1 - eps )
          {

            const int n = domain.cellsInDomain[ entityIndex ].size();

            Dune::DynamicMatrix< double > A( n, 2 );
            Dune::DynamicVector< double > b( n );

            DomainVector xi = entity.geometry().center();

            for( int i = 0; i < n; ++i )
            {
              int neighborIndex = domain.cellsInDomain[ entityIndex ][ i ];
              const auto &neighbor = grid.entity( domain.seeds[ neighborIndex ] );


              DomainVector xk = neighbor.geometry().center();
              double wk = 1.0 / ( (xk - xi).one_norm()  );
              wk *= wk;

              A[ i ][ 0 ] = wk * ( xk[ 0 ] - xi[ 0 ] );
              A[ i ][ 1 ] = wk * ( xk[ 1 ] - xi[ 1 ] );

              b[ i ] = wk * ( concentration[ neighborIndex ] - concentration[ entityIndex ] );
            }

            Dune::DynamicMatrix< double > AtA( 2, 2 );

            for( int i = 0; i < 2; ++i )
              for( int j = 0; j < 2; ++j )
              {
                AtA[ i ][ j ] = 0;

                for( int k = 0; k < n; ++k )
                  AtA[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
              }

            Dune::DynamicVector< double > x( 2 );
            Dune::DynamicVector< double > Atb( 2 );

            A.mtv( b, Atb );

            AtA.solve( x, Atb );

            if( x.one_norm() > 0 )
            {
              x *= 1.0 / x.two_norm();
              normals[ entityIndex ] = x;
            }
            else                                     // isolated cell, unmix it!
            {
              concentration[ entityIndex ] = ( concentration[ entityIndex ] < 0.5 ) ? 0.0 : 1.0;
              cellIsMixed[ entityIndex ] = false;
            }

          }
        }

      }

      ~SecondYoungsNormalGuessing ()
      {
        normals.clear();
      }

      void getNormal ( const int entityIndex, Line2D< DomainVector > &g )
      {
        g.n = normals[ entityIndex ];
      }

    private:

      std::vector< DomainVector > normals;

    };



  }       // end of namespace VoF
} // end of namespace Dune



#endif