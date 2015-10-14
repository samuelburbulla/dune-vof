#ifndef __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__
#define __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__

#include <dune/common/fvector.hh>


#include "geometricaltoolbox.hh"
#include "reconstruct.hh"


namespace Dune
{
  namespace VoF
  {

    template< class G, class V, class R, class D, class fvector >
    void evolve ( const G &grid, V &c, R &reconstruction, D &domain, const int numberOfCells, double t, double dt,
                  const double eps, std::vector< bool > &cellIsMixed,  std::vector< bool > &cellIsActive, const std::vector<fvector> &velocityField,
      std::vector< double > &overundershoots )
    {
      typedef typename G::LeafGridView GridView;


      GridView gridView = grid.leafGridView();


      // Calculate discrete velocity divergence in loop
      //std::vector< double > divergence( gridView.size( 0 ), 0 );


      std::map< std::pair<int,int>, double > intersectionFluxes;


      // Calculate volume fluxes for each cell
      V update( c.size() );
      for( std::size_t i = 0; i < c.size(); ++i )
        update[ i ] = 0;


      for( auto &&entity : elements( gridView ) )
      {

        int entityIndex = gridView.indexSet().index( entity );


        if( cellIsMixed[ entityIndex ] || cellIsActive[ entityIndex ] )
        {

          auto entityGeo = entity.geometry();


          for( auto &&intersection : intersections( gridView, entity ) )
          {
            if( intersection.neighbor() )
            {

              const auto &neighbor = intersection.outside();
              int neighborIndex = gridView.indexSet().index( neighbor );


              const auto isGeo = intersection.geometry();

              fvector outerNormal = intersection.centerUnitOuterNormal();

              fvector velocity = velocityField[ entityIndex ];
              velocity += velocityField[ neighborIndex ];
              velocity *= 0.5;

              fvector ufEn = velocityField[ entityIndex ];
              ufEn *= c[ entityIndex ];

              fvector ufNe = velocityField[ neighborIndex ];
              ufNe *= c[ neighborIndex ];

              //divergence[ entityIndex ] += intersection.integrationOuterNormal( 0 ) * ( ufEn + ufNe ) * 0.5;

              std::pair<int,int> fluxIndex ( std::min( entityIndex, neighborIndex ), std::max( entityIndex, neighborIndex ) );

              // do not calculate fluxes twice
              //if ( intersectionFluxes.find( fluxIndex ) == intersectionFluxes.end() )
              //{
                // build time integration polygon
                velocity *= dt;


                Polygon2D< fvector > timeIntegrationPolygon;

                timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) );
                timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) );

                timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) - velocity );
                timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) - velocity );


                // outflows
                if( velocity * outerNormal > 0 )
                {

                  // build flux polygon with volume of material being fluxed
                  Polygon2D< fvector > fluxPolygon;
                  double fluxVolume = 0;

                  if( cellIsMixed[ entityIndex ] )
                  {
                    Line2D< fvector > reconstLine( reconstruction[ entityIndex ][ 2 ], reconstruction[ entityIndex ][ 0 ] );

                    polygonLineIntersection( timeIntegrationPolygon, reconstLine, fluxPolygon );
                    polyAddInnerVertices( timeIntegrationPolygon, reconstLine, fluxPolygon );

                    fluxVolume = fluxPolygon.volume();
                  }

                  else if( c[ entityIndex ] > 1 - eps )
                    fluxVolume = timeIntegrationPolygon.volume();


                  update[ entityIndex ] -= fluxVolume / entityGeo.volume();

                  intersectionFluxes.insert( std::pair<std::pair<int,int>,double>( fluxIndex , -fluxVolume ) );
                }

                //inflow
                else if( velocity * outerNormal < 0 )
                {
                  // build phase polygon of the neighbor
                  Polygon2D< fvector > fluxPolygon;
                  double fluxVolume = 0;

                  if( cellIsMixed[ neighborIndex ] )
                  {
                    Line2D< fvector > reconstLine( reconstruction[ neighborIndex ][ 2 ], reconstruction[ neighborIndex ][ 0 ] );

                    polygonLineIntersection( timeIntegrationPolygon, reconstLine, fluxPolygon );
                    polyAddInnerVertices( timeIntegrationPolygon, reconstLine, fluxPolygon );

                    fluxVolume = fluxPolygon.volume();
                  }
                  else if( c[ neighborIndex ] > 1 - eps )
                    fluxVolume = timeIntegrationPolygon.volume();

                  update[ entityIndex ] += fluxVolume / entityGeo.volume();

                  intersectionFluxes.insert( std::pair<std::pair<int,int>,double>( fluxIndex , fluxVolume ) );
                }
              //}
              //else // intersection was already calculated
              //{
                //update[ entityIndex ] -= intersectionFluxes.find( fluxIndex )->second / entityGeo.volume();
              //}
            } // end of intersection is neighbor
          }
        }
      }



      // Advance volume fractions f in time
      for( std::size_t i = 0; i < c.size(); ++i )
      {
        // discrete velocity divergence correction and advantage
        if ( cellIsMixed[ i ] || cellIsActive[ i ] )
        {
          //update[ i ] += divergence[ i ] * c[ i ] * dt * 0.5;
          c[ i ] += update[ i ];
          //c[ i ] /= 1.0 - divergence[ i ] * dt * 0.5;
        }
      }




      // Look if any volume fraction undershoots or overshoots
      for( std::size_t i = 0; i < c.size(); ++i )
      {

        if ( c [ i ] < 0 )
        {
          overundershoots[ i ] = c[ i ];
        }
        else if ( c[ i ] > 1 )
        {
          overundershoots[ i ] = c[ i ] - 1;
        }
        else
        {
          overundershoots[ i ] = 0;
        }

        // c[ i ] = std::max( 0.0, c[ i ] );
        // c[ i ] = std::min( 1.0, c[ i ] );
     }

    }


  }
}

#endif










