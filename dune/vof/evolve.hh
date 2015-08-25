#ifndef __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__
#define __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__

#include <dune/common/fvector.hh>


#include "geometricaltoolbox.hh"
#include "reconstruct.hh"


namespace Dune
{
  namespace VoF
  {

    template< class G, class V, class R, class D >
    void evolve ( const G &grid, V &c, R &reconstruction, D &domain, const int numberOfCells, double t, double dt,
                  const double eps, std::vector< bool > &cellIsMixed,  std::vector< bool > &cellIsActive )
    {
      const int dimworld = G::dimensionworld;
      typedef typename G::ctype ct;
      typedef typename Dune::FieldVector< ct, dimworld > fvector;
      typedef typename G::LeafGridView GridView;


      GridView gridView = grid.leafGridView();


      // Calculate discrete velocity divergence in loop 
      std::vector< double > divergence( gridView.size( 0 ) );


      std::map< std::pair<int,int>, double > intersectionFluxes;


      // Calculate volume fluxes for each cell
      V update( c.size() );
      for( std::size_t i = 0; i < c.size(); ++i )
        update[ i ] = 0;


      for( auto &&entity : elements( gridView ) )
      {

        int entityIndex = gridView.indexSet().index( entity );
        divergence[ entityIndex ] = 0;
	

        if( cellIsMixed[ entityIndex ] || cellIsActive[ entityIndex ] )
        {

          auto entityGeo = entity.geometry();


          for( auto &&intersection : intersections( gridView, entity ) )
          {

            const auto &neighbor = intersection.outside();
            int neighborIndex = gridView.indexSet().index( neighbor );


            const auto isGeo = intersection.geometry();

            fvector outerNormal = intersection.centerUnitOuterNormal();
            fvector velocity = psi( isGeo.center(), t );

            // discrete divergence
            //std::cout << "Calc.loop: " << divergence[ entityIndex ] << std::endl;

            std::pair<int,int> fluxIndex ( std::min( entityIndex, neighborIndex ), std::max( entityIndex, neighborIndex ) );
            
            // do not calculate fluxes twice
            if ( intersectionFluxes.find( fluxIndex ) == intersectionFluxes.end() )
            {

              auto neighborGeo = neighbor.geometry();
           


              //Compute edge volume fluxes
              double edgeVolumeFlux = 0.0;


              // build time integration polygon
              velocity *= dt;

              Polygon2D< fvector > timeIntegrationPolygon;

              timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) );
              timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) );

	      fvector flux = outerNormal;
	      flux *= ( velocity * outerNormal ); 
	      	
              timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) - flux );
              timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) - flux );

              // outflows
              if( velocity * outerNormal > 0 )
              {
 		//divergence[ entityIndex ] -= std::abs( velocity[0] * velocity[1] * 0.5 ) / entityGeo.volume();
		//divergence[ neighborIndex ] += std::abs( velocity[0] * velocity[1] * 0.5 ) / neighborGeo.volume();

                // build phase polygon of entity
                Polygon2D< fvector > phasePolygon;

                if( cellIsMixed[ entityIndex ] )
                {
                  phasePolygon.addVertex( reconstruction[ entityIndex ][ 0 ] );
                  phasePolygon.addVertex( reconstruction[ entityIndex ][ 1 ] );

                  Line2D< fvector > reconstLine( reconstruction[ entityIndex ][ 2 ], reconstruction[ entityIndex ][ 0 ] );
                  for( auto v : getInnerVertices( entityGeo, reconstLine ) )
                    phasePolygon.addVertex( v );
                }

                else if( c[ entityIndex ] > 1 - eps )
                  for( int i = 0; i < entityGeo.corners(); ++i )
                    phasePolygon.addVertex( entityGeo.corner( i ) );


                edgeVolumeFlux = polygonIntersectionVolume( timeIntegrationPolygon, phasePolygon );

                update[ entityIndex ] -= edgeVolumeFlux / entityGeo.volume();

                intersectionFluxes.insert( std::pair<std::pair<int,int>,double>( fluxIndex , -edgeVolumeFlux ) );
              }

              //inflow
              else if( velocity * outerNormal < 0 )
                if( intersection.neighbor() )
                {
		  //divergence[ entityIndex ] += std::abs( velocity[0] * velocity[1] * 0.5 ) / entityGeo.volume() ;
		  //divergence[ neighborIndex ] -= std::abs( velocity[0] * velocity[1] * 0.5 ) / neighborGeo.volume();
		  
		  // build phase polygon of the neighbor
                  Polygon2D< fvector > phasePolygon;

                  if( cellIsMixed[ neighborIndex ] )
                  {
                    phasePolygon.addVertex( reconstruction[ neighborIndex ][ 0 ] );
                    phasePolygon.addVertex( reconstruction[ neighborIndex ][ 1 ] );

                    Line2D< fvector > reconstLine( reconstruction[ neighborIndex ][ 2 ], reconstruction[ neighborIndex ][ 0 ] );
                    for( auto v : getInnerVertices( neighbor.geometry(), reconstLine ) )
                      phasePolygon.addVertex( v );
                  }
                  else if( c[ neighborIndex ] > 1 - eps )
                    for( int i = 0; i < neighborGeo.corners(); ++i )
                      phasePolygon.addVertex( neighborGeo.corner( i ) );

                  edgeVolumeFlux = polygonIntersectionVolume( timeIntegrationPolygon, phasePolygon );

                  update[ entityIndex ] += edgeVolumeFlux / entityGeo.volume();

                  intersectionFluxes.insert( std::pair<std::pair<int,int>,double>( fluxIndex , edgeVolumeFlux ) );
                }
            }
            else // intersection was already calculated 
            {
              update[ entityIndex ] -= intersectionFluxes.find( fluxIndex )->second / entityGeo.volume(); 
	    }
            
          }
        }
      }



      // Advance volume fractions f in time
      for( std::size_t i = 0; i < c.size(); ++i )
      {
        // discrete velocity divergence correction and advantage
        if ( cellIsMixed[ i ] || cellIsActive[ i ] )
        {
          //update[ i ] -= c[ i ] * divergence[ i ] * 0.5;

	  c[ i ] += update[ i ];
	  //c[ i ] /= ( 1 - divergence[ i ] * 0.5 );
	}
      }




      // Look if any volume fraction undershoots or overshoots
      for( unsigned int i = 0; i < c.size(); ++i )
      {
        c[ i ] = std::max( 0.0, c[ i ] );
        c[ i ] = std::min( 1.0, c[ i ] );
      }

    }


  }
}

#endif










