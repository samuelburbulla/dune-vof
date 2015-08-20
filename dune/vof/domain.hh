#ifndef __DUNE_GRID_REC_VOL_TRACK_DOMAIN_HH__
#define __DUNE_GRID_REC_VOL_TRACK_DOMAIN_HH__

namespace Dune
{
  namespace VoF
  {


    template< class Grid >
    class DomainOfCells
    {

    public:

      const int dimworld = Grid::dimensionworld;
      typedef typename Grid::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Grid::template Codim< 0 >::EntitySeed EntitySeed;


      DomainOfCells ( const Grid &grid )
      {



        int n = grid.leafGridView().size( 0 );

        cellsInDomain = std::vector < std::vector < int>> ( n );
        seeds = std::vector< EntitySeed >( n );

        getCellsInDomain( grid );

      }

      std::vector < std::vector < int>> cellsInDomain;
      std::vector< EntitySeed > seeds;

    private:


      void getCellsInDomain ( const Grid &grid )
      {

        GridView gridView = grid.leafGridView();

        std::vector < std::vector < int>> cellsNextToThisVertex( gridView.size( dimworld ) );

        for( LeafIterator leaf = gridView.template begin< 0 >(); leaf != gridView.template end< 0 >(); ++leaf )
        {
          const Entity &entity = *leaf;
          int entityID = gridView.indexSet().index( entity );

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexID = gridView.indexSet().subIndex( entity, k, dimworld );
            cellsNextToThisVertex[ vertexID ].push_back( entityID );
          }
        }


        for( LeafIterator leaf = gridView.template begin< 0 >(); leaf != gridView.template end< 0 >(); ++leaf )
        {
          const Entity &entity = *leaf;
          int entityID = gridView.indexSet().index( entity );

          seeds[ entityID ] = entity.seed();

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexIndex = gridView.indexSet().subIndex( entity, k, dimworld );

            for( auto index : cellsNextToThisVertex[ vertexIndex ] )
              if( index != entityID )
                cellsInDomain[ entityID ].push_back( index );
          }


          // erase duplicate elements
          std::sort( cellsInDomain[ entityID ].begin(), cellsInDomain[ entityID ].end());
          auto last = std::unique( cellsInDomain[ entityID ].begin(), cellsInDomain[ entityID ].end() );
          cellsInDomain[ entityID ].erase( last, cellsInDomain[ entityID ].end());
        }

      }


    };

  }
}

#endif