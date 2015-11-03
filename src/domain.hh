#ifndef DUNE_GRID_REC_VOL_TRACK_DOMAIN_HH
#define DUNE_GRID_REC_VOL_TRACK_DOMAIN_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {


    template< class GridView >
    class DomainOfPointNeighbors
    {

    public:

      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::template Codim< 0 >::Entity Entity;
      typedef typename Entity::EntitySeed EntitySeed;


      explicit DomainOfPointNeighbors ( const GridView &gridView )  
        : _gridView ( gridView ), _mapper ( gridView ), _cellsInDomain ( _mapper.size() ), _seeds ( _mapper.size() )
      {       

        getCellsInDomain();

      }

      const std::vector< Entity > operator[] ( const Entity &entity ) const    // range-based!
      {
        std::vector< Entity > stencil;

        for ( auto&& index : _cellsInDomain[ _mapper.index( entity ) ] )
        {
          stencil.push_back( _gridView.grid().entity( _seeds[ index ] ) );
        }

        return stencil;
      }





    private:

      void getCellsInDomain()
      {

        std::vector < std::vector < int > > cellsNextToThisVertex( _gridView.size( dimworld ) );

        for( auto&& entity : elements( _gridView ) )
        {
          int entityID = _mapper.index( entity );
          _seeds[ entityID ] = entity.seed();

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexID = _gridView.indexSet().subIndex( entity, k, dimworld );  //mapper?
            cellsNextToThisVertex[ vertexID ].push_back( entityID );
          }
        }


        for( auto&& entity : elements( _gridView ) )
        {
          int entityID = _mapper.index( entity );
         
          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexIndex = _gridView.indexSet().subIndex( entity, k, dimworld );  //mapper?

            for( auto index : cellsNextToThisVertex[ vertexIndex ] )
              if( index != entityID )
                _cellsInDomain[ entityID ].push_back( index );
          }


          // erase duplicate elements
          std::sort( _cellsInDomain[ entityID ].begin(), _cellsInDomain[ entityID ].end());
          auto last = std::unique( _cellsInDomain[ entityID ].begin(), _cellsInDomain[ entityID ].end() );
          _cellsInDomain[ entityID ].erase( last, _cellsInDomain[ entityID ].end());
        }

      }

      GridView _gridView;
      Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > _mapper;

      std::vector < std::vector < int > > _cellsInDomain;
      std::vector < EntitySeed > _seeds;


    };

  }
}

#endif