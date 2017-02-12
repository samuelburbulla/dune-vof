#ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH
#define DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH

namespace Dune
{
  namespace VoF
  {

    // VertexNeighborsStencil
    // ----------------------

     /**
     * \ingroup Method
     * \brief  set of vertex neighbors stencils
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct VertexNeighborsStencil
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Stencil = std::vector< Entity >;
      static constexpr int dim = GridView::dimension;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );

    public:
      explicit VertexNeighborsStencil ( const GridView& gridView )
       : gridView_( gridView ), stencils_( indexSet().size( 0 ) )
      {
        initialize();
      }

      const Stencil& operator[] ( const Entity& entity ) const
      {
        return stencils_[ indexSet().index( entity ) ];
      }

      const GridView& gridView() const { return gridView_; }
    private:
      const IndexSet& indexSet() const { return gridView().indexSet(); }

      void initialize()
      {
        std::vector< std::vector< Index > > cellsNextToThisVertex( indexSet().size( dim ) );
        std::vector< std::vector< Index > > cellsInDomain_( indexSet().size( 0 ) );
        std::vector< Entity > entities_( indexSet().size( 0 ) );

        for( const Entity& entity : elements( gridView(), Partitions::all ) )
        {
          Index id = indexSet().index( entity );
          entities_[ id ] = entity;

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            Index vId = indexSet().subIndex( entity, k, dim );
            cellsNextToThisVertex[ vId ].push_back( id );
          }
        }

        for( const Entity& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          Index id = indexSet().index( entity );

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            Index vId = indexSet().subIndex( entity, k, dim );

            for( Index index : cellsNextToThisVertex[ vId ] )
              if( index != id )
                cellsInDomain_[ id ].push_back( index );
          }


          // Erase duplicate elements
          std::sort( cellsInDomain_[ id ].begin(), cellsInDomain_[ id ].end() );
          auto last = std::unique( cellsInDomain_[ id ].begin(), cellsInDomain_[ id ].end() );
          cellsInDomain_[ id ].erase( last, cellsInDomain_[ id ].end() );

          // Create stencil object
          Stencil& stencil = stencils_[ id ];
          for ( Index index : cellsInDomain_[ id ] )
            stencil.push_back( entities_[ index ] );
        }
      }

      GridView gridView_;
      std::vector< Stencil > stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH
