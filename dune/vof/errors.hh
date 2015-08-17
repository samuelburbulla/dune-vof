#ifndef __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__
#define __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__ 
 
#include <vector>

namespace Dune
{
  namespace VoF
  {
 
    template <class G>
    double l2error( const G& grid, const std::vector<double>& cStart, const std::vector<double>& cEnd )  
    {
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim<0>::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;
      
      double error = 0;
      
      GridView gridView = grid.leafGridView();
      
      for (LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf)
      {
	const Entity &entity = *leaf;
	const Geometry geo = entity.geometry();
	
	int i = gridView.indexSet().index( entity );
    
	error += (cStart[i] - cEnd[i]) * (cStart[i] - cEnd[i]) * geo.volume();   
      }
      return error;
    }


    template <class G>
    double l1error( const G& grid, const std::vector<double>& cStart, const std::vector<double>& cEnd )  
    {
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim<0>::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;
      
      double error = 0;
      
      GridView gridView = grid.leafGridView();
      
      for (LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf)
      {
	const Entity &entity = *leaf;
	const Geometry geo = entity.geometry();
	
	int i = gridView.indexSet().index( entity );
    
	error += fabs(cStart[i] - cEnd[i]) * geo.volume();   
      }
      return error;
    }


  }
}

#endif