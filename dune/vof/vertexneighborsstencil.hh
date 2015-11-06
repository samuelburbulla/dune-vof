#ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH
#define DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH

//- dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>


namespace Dune
{
  namespace VoF
  {


    // VertexNeighborsStencil
    // ----------------------

    template< class GV >
    struct VertexNeighborsStencil
    {
      using GridView = GV;
      using Entity = typename GridView::template Codim< 0 >::Entity;
      using EntitySeed = typename Entity::EntitySeed;

      static constexpr int dim = GridView::dimensionworld;

    private:
      using Mapper =
        Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout >;

    public:
      explicit VertexNeighborsStencil ( const GridView &gridView )
        : gridView_( gridView ), mapper_( gridView_ ), cellsInDomain_( mapper_.size() ), seeds_( mapper_.size() )
      {
        getCellsInDomain();
      }

      const std::vector< Entity > operator[] ( const Entity &entity ) const    // range-based!
      {
        std::vector< Entity > stencil;

        for ( const auto& index : cellsInDomain_[ mapper().index( entity ) ] )
          stencil.push_back( gridView().grid().entity( seeds_[ index ] ) );

        return stencil;
      }

    private:
      const GridView &gridView () const { return gridView_; }
      const Mapper &mapper () const { return mapper_; }

      void getCellsInDomain()
      {
        std::vector < std::vector < int > > cellsNextToThisVertex( gridView().size( dim ) );

        for( const auto &entity : elements( gridView() ) )
        {
          int entityID = mapper().index( entity );
          seeds_[ entityID ] = entity.seed();

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexID = gridView().indexSet().subIndex( entity, k, dim );  //mapper?
            cellsNextToThisVertex[ vertexID ].push_back( entityID );
          }
        }
        for( const auto &entity : elements( gridView() ) )
        {
          int entityID = mapper().index( entity );

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            int vertexIndex = gridView().indexSet().subIndex( entity, k, dim );  //mapper?

            for( auto index : cellsNextToThisVertex[ vertexIndex ] )
              if( index != entityID )
                cellsInDomain_[ entityID ].push_back( index );
          }


          // erase duplicate elements
          std::sort( cellsInDomain_[ entityID ].begin(), cellsInDomain_[ entityID ].end() );
          auto last = std::unique( cellsInDomain_[ entityID ].begin(), cellsInDomain_[ entityID ].end() );
          cellsInDomain_[ entityID ].erase( last, cellsInDomain_[ entityID ].end() );
        }
      }

      GridView gridView_;
      Mapper mapper_;

      std::vector< std::vector< int > > cellsInDomain_;
      std::vector< EntitySeed > seeds_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH