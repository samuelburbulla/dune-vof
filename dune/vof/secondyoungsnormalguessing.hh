#ifndef DUNE_VOF_SECONGYOUNGSNORMALGUESSING_HH
#define DUNE_VOF_SECONGYOUNGSNORMALGUESSING_HH

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>


namespace Dune
{
  namespace VoF
  {

      template< class GridView, class ColorFunction, class Flags, class Domain, class ReconstructionSet >
      void SecondYoungsNormalGuessing ( const GridView &gridView, const ColorFunction &colorFunction, const Flags &flags, const Domain &domain,
                                        ReconstructionSet &reconstructionSet )
      {

        const int dimworld = GridView::dimensionworld;
        typedef typename GridView::ctype ctype;
        typedef typename Dune::FieldVector< ctype, dimworld > DomainVector;


        for( auto&& entity : elements( gridView ) )
        {
          if( flags.isMixed( entity ) )
          {

            const int n = domain[ entity ].size();

            Dune::DynamicMatrix< double > A( n, 2 );
            Dune::DynamicVector< double > b( n );

            DomainVector xi = entity.geometry().center();

            for( std::size_t i = 0; i < domain[ entity ].size(); ++i )
            {
              auto neighbor = domain[ entity ][ i ];
              DomainVector xk = neighbor.geometry().center();
              double wk = 1.0 / ( (xk - xi).two_norm2()  );

              A[ i ]  = xk;
              A[ i ] -= xi;
              A[ i ] *= wk;

              b[ i ] = wk * ( colorFunction[ neighbor ] - colorFunction[ entity ] );
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

            assert( x.two_norm() > 0 );

            x *= 1.0 / x.two_norm();
            reconstructionSet[ entity ].normal() = x;

          }
        }

      }

  }       // end of namespace VoF
} // end of namespace Dune



#endif