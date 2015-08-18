#ifndef DUNE_VOF_FLAGCELLS_HH
#define DUNE_VOF_FLAGCELLS_HH
 
#include <vector>

namespace Dune
{
  namespace VoF
  {
    
    template < class GV, class C, class R, class D >
    void flagCells ( const GV& gridView, const C& c, R& reconstruction, const D& domain, const int numberOfCells,
		     std::vector<bool> &cellIsMixed,  std::vector<bool> &cellIsActive, const double eps )
    {  
      typedef typename Dune::FieldVector<double, 2> fvector;
      typedef typename GV::template Codim<0>::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;
      typedef typename GV::Intersection Intersection;
    
      // clear cellIsMixed and cellIsActive vectors
      for ( std::size_t i = 0; i < cellIsMixed.size(); ++i )
	cellIsMixed[i] = false;
      for ( std::size_t i = 0; i < cellIsActive.size(); ++i )
	cellIsActive[i] = false;
     
      
      // mixed cells
      for (LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf)
      {
	
	const Entity &entity = *leaf;
	int entityIndex = gridView.indexSet().index( entity );
	
	if ( c[ entityIndex ] >= eps && c[ entityIndex ] <= 1 - eps )

	  cellIsMixed[ entityIndex ] = true;     
		
	else if ( c[ entityIndex ] > 1 - eps )
	{
	  for ( IntersectionIterator is = gridView.ibegin( entity ); is != gridView.iend( entity ); ++is )
	  {
	    const Intersection &intersection = *is;
	    const Entity &neighbor = intersection.outside();
	    int neighborIndex = gridView.indexSet().index( neighbor );
	    
	    if ( c[ neighborIndex ] < eps )
	    {
	      cellIsMixed[ entityIndex ] = true;
		      
	      const IntersectionGeometry isGeo = intersection.geometry();
	      auto n = intersection.centerUnitOuterNormal();
	      n *= -1.0;
	      
	      reconstruction[ entityIndex ] = std::array<fvector,3> ( { isGeo.corner(0), isGeo.corner(1), n } ); 
	      continue;
	    }
	  }
	}
      }
      
      // activate cells
      for ( std::size_t entityIndex = 0; entityIndex < cellIsMixed.size(); ++entityIndex )
      {
	if ( cellIsMixed[ entityIndex ] )
	  for ( int neighborIndex : domain.cellsInDomain[ entityIndex ] )
	    if ( !cellIsMixed[ neighborIndex ] )
	      cellIsActive[ neighborIndex ] = true;
      }
      
    }


  }
}

#endif