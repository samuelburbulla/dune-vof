#ifndef __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__
#define __DUNE_GRID_REC_VOL_TRACK_EVOLVE_HH__

#include <dune/common/fvector.hh>


#include "geometricaltoolbox.hh"
#include "reconstruct.hh"


namespace Dune
{
  namespace VoF
  {
	
    template < class G, class V, class R, class D >
    void evolve ( const G& grid, V& c, R& reconstruction, D& domain, const int numberOfCells, double t, double dt, 
		  const double eps, std::vector<bool>& cellIsMixed,  std::vector<bool>& cellIsActive )
    {
	    const int dimworld = G::dimensionworld;
	    typedef typename G::ctype ct;
	    typedef typename Dune::FieldVector<ct, dimworld> fvector;
	    typedef typename G::LeafGridView GridView;
	    typedef typename GridView::template Codim<0>::Iterator LeafIterator;
	    typedef typename LeafIterator::Entity Entity;
	    typedef typename Entity::Geometry Geometry;
	    typedef typename GridView::IntersectionIterator IntersectionIterator;
	    typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;
	    typedef typename GridView::Intersection Intersection;
	      
	    GridView gridView = grid.leafGridView();
	      
	      
	    // Calculate discrete velocity divergence
	    std::vector<double> divergence ( gridView.size(0) );
	      
	    for (LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf)
	    {
		
			const Entity &entity = *leaf;
			int entityIndex = gridView.indexSet().index( entity );
			
			divergence[ entityIndex ] = 0; 
			for ( IntersectionIterator is = gridView.ibegin( entity ); is != gridView.iend( entity ); ++is )
			{
				const IntersectionGeometry isGeo = is->geometry();
			  
				divergence[ entityIndex ] += is->centerUnitOuterNormal() * psi( isGeo.center(), t ) * isGeo.volume();
			}
		
	    }
	      
	    
	      
	    // Calculate volume fluxes for each cell
	    V update( c.size() );
	    for( std::size_t i = 0; i < c.size(); ++i ) update[i] = 0;
      
      
      	for (LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf )
      	{
      
			const Entity &entity = *leaf;
			int entityIndex = gridView.indexSet().index( entity );
	     
			if ( cellIsMixed[ entityIndex ] || cellIsActive[ entityIndex ] )
			{
		
		  
		  
			  	const Geometry geo = entity.geometry();
				double entityVolume = geo.volume();
				    
				  
				for ( IntersectionIterator is = gridView.ibegin( entity ); is != gridView.iend( entity ); ++is )
				{
				  	const Intersection &intersection = *is;
				  	const Entity &neighbor = intersection.outside();
				    
				  	int neighborIndex = gridView.indexSet().index( neighbor );  
				  	Geometry neighborGeo = neighbor.geometry();
				    
				    
				    
				  	const IntersectionGeometry isGeo = intersection.geometry();
				  
				  	fvector velocity = psi( isGeo.center(), t );
				  	velocity *= dt;
			    
				  	assert( std::abs(velocity[0]) < 1.0/numberOfCells && fabs(velocity[1]) < 1.0/numberOfCells );
				    
				  
				    //Compute edge volume fluxes       
				    double edgeVolumeFlux = 0.0;
				  
				    fvector outerNormal = intersection.centerUnitOuterNormal();
				    
				    
				    // build time integration polygon
				    Polygon2D<fvector> timeIntegrationPolygon;
				    
				    timeIntegrationPolygon.addVertex( isGeo.corner(0) );
				    timeIntegrationPolygon.addVertex( isGeo.corner(1) );
				    
				    timeIntegrationPolygon.addVertex( isGeo.corner(0) - velocity );
				    timeIntegrationPolygon.addVertex( isGeo.corner(1) - velocity );
				    
			      
				      //outflow
				    if ( velocity * outerNormal > 0 ) 
				    {

						// build phase polygon of entity  
						Polygon2D<fvector> phasePolygon;
				      
						if ( cellIsMixed[ entityIndex ] )
						{
						  	phasePolygon.addVertex( reconstruction[ entityIndex ][0] );
						  	phasePolygon.addVertex( reconstruction[ entityIndex ][1] );
					      
						  	Line2D<fvector> reconstLine ( reconstruction[ entityIndex ][2], reconstruction[ entityIndex ][0] );
						  	for ( auto v : getInnerVertices( geo, reconstLine ) )
						    	phasePolygon.addVertex( v );
						} 
					      
						else if ( c[ entityIndex ] > 1 - eps )
						{
						  	for ( int i = 0; i < geo.corners(); ++i )
						    	phasePolygon.addVertex( geo.corner(i) );
						}

						
						edgeVolumeFlux = polygonIntersectionVolume( timeIntegrationPolygon, phasePolygon );  
						
						update[entityIndex] -= edgeVolumeFlux / entityVolume;
				      
				    }
				  
				    //inflow
				    else if ( velocity * outerNormal < 0 ) 
				    {  
				      	if ( is -> neighbor() )
				      	{ 
							// build phase polygon of the neighbor
							Polygon2D<fvector> phasePolygon;
						      
							if ( cellIsMixed[ neighborIndex ] )
							{
							    phasePolygon.addVertex( reconstruction[ neighborIndex ][0] );
							    phasePolygon.addVertex( reconstruction[ neighborIndex ][1] );
							      
							    Line2D<fvector> reconstLine ( reconstruction[ neighborIndex ][2], reconstruction[ neighborIndex ][0] );
							    for ( auto v : getInnerVertices( neighbor.geometry(), reconstLine ) )
							        phasePolygon.addVertex( v );
							}
							else if ( c[ neighborIndex ] > 1 - eps )
							{
							    for ( int i = 0; i < neighborGeo.corners(); ++i )
								    phasePolygon.addVertex( neighborGeo.corner(i) );
							}
						    
							edgeVolumeFlux = polygonIntersectionVolume( timeIntegrationPolygon, phasePolygon );  
							      
							update[entityIndex] += edgeVolumeFlux / entityVolume;    

				      	}
				    }
				  
				}
			  

			}
      	}


      
      	// Advance volume fractions f in time
      	for ( std::size_t i = 0; i < c.size(); ++i)
      	{
		  	// discrete velocity divergence correction and advantage
		  	update[i] += c[i] * dt * divergence[i] * 0.5;
		  	c[i] += update[i] / ( 1 - dt * divergence[i] * 0.5 );
		  
      	}
      
      	
	
	
      	// Look if any volume fraction undershoots or overshoots
      	for (unsigned int i = 0; i < c.size(); ++i)
      	{
			c[i] = std::max( 0.0, c[i] );
			c[i] = std::min( 1.0, c[i] );
      	}
      
      	return;
    }


  }
}

#endif










