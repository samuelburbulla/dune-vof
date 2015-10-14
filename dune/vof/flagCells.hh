#ifndef DUNE_VOF_FLAGCELLS_HH
#define DUNE_VOF_FLAGCELLS_HH

#include <vector>

namespace Dune
{
  namespace VoF
  {

    template< class GV, class C, class R, class D >
    void flagCells ( const GV &gridView, const C &c, R &reconstruction, const D &domain,
                     std::vector< bool > &cellIsMixed,  std::vector< bool > &cellIsActive, const double eps )
    {
      typedef typename Dune::FieldVector< double, 2 > fvector;
      typedef typename GV::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection::Geometry IntersectionGeometry;

      // clear cellIsMixed and cellIsActive vectors
      for( std::size_t i = 0; i < cellIsMixed.size(); ++i )
        cellIsMixed[ i ] = false;

      for( std::size_t i = 0; i < cellIsActive.size(); ++i )
        cellIsActive[ i ] = false;



      // mixed cells
      for( auto&& entity : elements( gridView ) )
      {

        int entityIndex = gridView.indexSet().index( entity );


        if( c[ entityIndex ] >= eps && c[ entityIndex ] <= 1 - eps )
          cellIsMixed[ entityIndex ] = true;

        else if( c[ entityIndex ] > 1 - eps )
          for( auto&& intersection : intersections( gridView, entity ) )
          {
            if ( intersection.neighbor() )
            {
              const Entity &neighbor = intersection.outside();
              int neighborIndex = gridView.indexSet().index( neighbor );

              if( c[ neighborIndex ] < eps )
              {
                cellIsMixed[ entityIndex ] = true;

                const IntersectionGeometry isGeo = intersection.geometry();
                auto n = intersection.centerUnitOuterNormal();
                n *= -1.0;

                reconstruction[ entityIndex ] = std::array< fvector, 3 >( {{ isGeo.corner( 0 ), isGeo.corner( 1 ), n }} );
              }
            }
          }
      }

      // activate cells
      for( auto&& entity : elements( gridView ) )
      {
        int entityIndex = gridView.indexSet().index( entity );

        if( cellIsMixed[ entityIndex ] )
          for( auto&& is : intersections( gridView,  entity ) )
          {
            if( is.neighbor() )
            {
              int neighborIndex = gridView.indexSet().index( is.outside() );

              if( !cellIsMixed[ neighborIndex ] )
                  cellIsActive[ neighborIndex ] = true;
            }
          }
      }
    }

  } // end of namespace VoF
} // end of namespace Dune

#endif
