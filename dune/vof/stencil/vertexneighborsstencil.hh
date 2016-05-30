#ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH
#define DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

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
      using Mapper = Dune::VoF::MCMGMapper< GridView, Dune::MCMGElementLayout >;
      using VertexMapper = Dune::VoF::MCMGMapper< GridView, Dune::MCMGVertexLayout >;

    public:
      explicit VertexNeighborsStencil ( const GridView& gridView )
       : gridView_( gridView ), mapper_( gridView_ ), vmapper_( gridView_ ), stencils_( mapper_.size() )
      {
        initialize();
      }

      const Stencil& operator[] ( const Entity& entity ) const
      {
        return stencils_[ mapper().index( entity ) ];
      }

    private:
      const GridView& gridView() const { return gridView_; }
      const Mapper& mapper() const { return mapper_; }
      const VertexMapper& vmapper() const { return vmapper_; }

      void initialize()
      {
        std::vector< std::vector< std::size_t > > cellsNextToThisVertex( gridView().indexSet().size( dim ) );
        std::vector< std::vector< std::size_t > > cellsInDomain_( mapper_.size() );
        std::vector< Entity > entities_( mapper_.size() );

        for( const Entity& entity : elements( gridView(), Partitions::all ) )
        {
          std::size_t id = mapper().index( entity );
          entities_[ id ] = entity;

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            std::size_t vId = vmapper().subIndex( entity, k, dim );
            cellsNextToThisVertex[ vId ].push_back( id );
          }
        }

        for( const Entity& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          std::size_t id = mapper().index( entity );

          for( int k = 0; k < entity.geometry().corners(); k++ )
          {
            std::size_t vId = vmapper().subIndex( entity, k, dim );

            for( std::size_t index : cellsNextToThisVertex[ vId ] )
              if( index != id )
                cellsInDomain_[ id ].push_back( index );
          }


          // Erase duplicate elements
          std::sort( cellsInDomain_[ id ].begin(), cellsInDomain_[ id ].end() );
          auto last = std::unique( cellsInDomain_[ id ].begin(), cellsInDomain_[ id ].end() );
          cellsInDomain_[ id ].erase( last, cellsInDomain_[ id ].end() );

          // Create stencil object
          Stencil& stencil = stencils_[ id ];
          for ( const std::size_t index : cellsInDomain_[ id ] )
            stencil.push_back( entities_[ index ] );
        }
      }

      GridView gridView_;
      Mapper mapper_;
      VertexMapper vmapper_;
      std::vector< Stencil > stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_VERTEXNEIGHBORSSTENCIL_HH
