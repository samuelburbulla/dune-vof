#ifndef __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__
#define __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__

#include <stdexcept>
#include <stdio.h>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>

#include "geometricaltoolbox.hh"
#include "secondyoungsnormalguessing.hh"


namespace Dune
{
  namespace VoF
  {

    template< class GV, class E, class V >
    double computeInterfaceLinePosition( const GV &gridView, const E &entity, const V &n, double concentration );

    template< class GV, class E, class V >
    V interfaceLineCentroid( const GV &gridView, const E &entity, const V &normal, const double pEntity );

    double brentsMethod ( double a, double b, double (*f) ( double p ) );



    template< class R >
    void clearReconstruction ( R &reconstruction )
    {
      typedef typename Dune::FieldVector< double, 2 > fvector;
      fvector null; null = 0;
      for( std::size_t i = 0; i < reconstruction.size(); ++i )
        reconstruction[ i ] = std::array< fvector, 3 >( { null, null, null } );
    }



    template< class Grid, class C, class R, class F, class D >
    void reconstruct ( const Grid &grid, C &concentration, R &reconstruction, F &cellIsMixed, const D &domain, const double eps )
    {
      const int dimworld = Grid::dimensionworld;
      typedef typename Dune::FieldVector< double, dimworld > fvector;
      typedef typename Grid::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;


      Dune::VoF::SecondYoungsNormalGuessing< Grid, fvector > guessedNormals( grid, concentration, cellIsMixed, domain, eps );


      int count = 0;

      const GridView gridView = grid.leafGridView();

      fvector normalOld;
      Line2D< fvector > g;

      for( auto&& entity : elements( gridView ) )
      {
        count = 0;

        int entityIndex = gridView.indexSet().index( entity );

        if( cellIsMixed[ entityIndex ] && concentration[ entityIndex ] <= 1 - eps )
        {
          auto geoEntity = entity.geometry();


          guessedNormals.getNormal( entityIndex, g );



          double sumCount;
          fvector centroidLineNormal, sumNormals;

          do
          {
            sumCount = 0;
            sumNormals = 0;


            auto entityIntersections = computeInterfaceLinePosition( gridView, entity, geoEntity, concentration[ entityIndex ], g );


            fvector centroid1;
            // reconstruction doesn't intersect
            if( entityIntersections.size() == 0 )
            {
              cellIsMixed[ entityIndex ] = false;
              break;
            }
            // get middle of intersections
            else
            {
              int c = 0;
              for( std::size_t i = 0; i < entityIntersections.size(); ++i )
              {
                centroid1 += entityIntersections[ i ];
                c++;
              }
              centroid1 *= 1.0 / c;
            }



            for( int neighborIndex : domain.cellsInDomain[ entityIndex ] )

              if( cellIsMixed[ neighborIndex ] && concentration[ neighborIndex ] <= 1 - eps )
              {

                const Entity &neighbor = grid.entity( domain.seeds[ neighborIndex ] );

                Line2D< fvector > h;

                // nimm diese Zelle nur, wenn die geratene Normal schon in eine aehnliche Richtung zeigt
                guessedNormals.getNormal( neighborIndex, h );

                if( g.n * h.n < 0 )
                  continue;



                auto geoNeighbor = neighbor.geometry();

                h.n = g.n;

                auto neighborIntersections = computeInterfaceLinePosition( gridView, neighbor, geoNeighbor, concentration[ neighborIndex ], h );

                fvector centroid2;
                // reconstruction doesn't intersect
                if( neighborIntersections.size() == 0 )
                {
                  cellIsMixed[ neighborIndex ] = false;
                  continue;
                }
                // get middle of intersections
                else
                {
                  int c = 0;
                  for( std::size_t i = 0; i < neighborIntersections.size(); ++i )
                  {
                    centroid2 += neighborIntersections[ i ];
                    c++;
                  }
                  centroid2 *= 1.0 / c;
                }

                centroidLineNormal = centroid2 - centroid1;
                rotate90degreesCounterClockwise( centroidLineNormal );


                assert( centroidLineNormal.two_norm() > 0 );



                // errechnete Normale muss in eine aehnliche Richtung wie geschaetze Normale auf jeder Zelle zeigen
                if( centroidLineNormal * g.n < 0 )
                  centroidLineNormal *= -1.0;


                sumNormals.axpy( 1.0 / (geoEntity.center() - centroid2).two_norm(), centroidLineNormal );
                sumCount++;
              }



            normalOld = g.n;

            if( sumCount > 0 )
            {
              g.n = sumNormals;
              g.n *= 1.0 / sumCount;
            }

            assert( g.n.two_norm() > 1e-14 );

            count++;

          } while( (g.n - normalOld).two_norm2() > 1e-8 && count < 30 );            // limit number of loops
          

          auto entityIntersections
            = computeInterfaceLinePosition( gridView, entity, geoEntity, concentration[ entityIndex ], g );

          if( entityIntersections.size() != 0 )
            reconstruction[ entityIndex ] = std::array< fvector, 3 >( { entityIntersections[ 0 ], entityIntersections[ 1 ], g.n } );

        }
      }

    }




    template< class GV, class E, class Geo, class V >
    std::vector< V > computeInterfaceLinePosition ( const GV &gridView, const E &entity, const Geo &geo, double concentration, Line2D< V > &g )
    {
      double pMin = 0, pMax = 0, volume;
      //use bigger range than [0,1] initially
      double volMin = -1;
      double volMax = 2;

      // Initial guess for p
      for( int i = 0; i < geo.corners(); ++i )
      {
        g.p = geo.corner( i ) * g.n;
        g.p = -g.p;

        volume = getVolumeFraction( gridView, entity, geo, g );


        if( volume <= volMax && volume >= concentration )
        {
          pMax = g.p;
          volMax = volume;
        }

        if( volume >= volMin && volume <= concentration )
        {
          pMin = g.p;
          volMin = volume;
        }
      }


      g.p = brentsMethod( pMin, pMax,
                          [ &gridView, &entity, &geo, &concentration, &g ] ( double p ) -> double
                          { Line2D< V > h( g.n, p ); return getVolumeFraction( gridView, entity, geo, h ) - concentration; } );

      return lineCellIntersections( gridView, entity, geo, g );
    }



    //Quelle: Wikipedia (Brent-Verfahren)
    template< class F >
    double brentsMethod ( double a, double b, F f )
    {
      const double TOL = 1e-14;

      double fa, fb, fc, c, d, e, p, q, m, s, tol, r;


      fa = f( a );
      fb = f( b );

      if( std::abs( fa ) < TOL )
        return a;
      if( std::abs( fb ) < TOL )
        return b;

      assert( fa * fb < 0 );

      c = a; fc = fa;
      d = b - a; e = d;

      int iter = 0;
      int maxiter = 1000;

      while( iter < maxiter )
      {

        iter++;

        if( fb * fc > 0 )
        {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }

        if( std::abs( fc ) < std::abs( fb ) )
        {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }

        tol = 2.0 * 1e-8 * std::abs( b ) + TOL;
        m = ( c - b ) / 2.0;

        if( std::abs( m ) > tol && std::abs( fb ) > 0 )
        {

          if( std::abs( e ) < tol || std::abs( fa ) <= std::abs( fb ) )
          {
            d = m;
            e = m;
          }
          else
          {
            s = fb / fa;

            if( a == c )
            {
              p = 2.0 * m * s;
              q = 1 - s;
            }
            else
            {
              q = fa / fc;
              r = fb / fc;
              p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
              q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
            }

            if( p > 0 )
              q = -q;
            else
              p = -p;

            s = e;
            e = d;

            if( 2.0 * p < 3.0 * m * q - std::abs( tol * q ) && p < std::abs( s * q / 2.0 ) )
              d = p / q;
            else
            {
              d = m;
              e = m;
            }
          }

          a = b;
          fa = fb;

          if( std::abs( d ) > tol )
            b = b + d;
          else
          {
            if( m > 0 )
              b = b + tol;
            else
              b = b - tol;
          }
        }
        else
          break;

        fb = f( b );
      }

      //std::cerr << "Brent steps:" << iter << std::endl;
      return b;

    }


  }
}

#endif

