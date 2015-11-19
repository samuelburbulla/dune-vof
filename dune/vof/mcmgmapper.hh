#ifndef DUNE_VOF_MCMGMAPPER_HH
#define DUNE_VOF_MCMGMAPPER_HH

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{

  namespace VoF
  {

    // MCMGMapper
    // ----------

    template< typename GV, template< int > class Layout >
    class MCMGMapper
    {
      typedef MCMGMapper< GV, Layout > This;

    public:
      typedef GV GridView;

      static const int dimension = GridView::dimension;

      typedef typename std::decay< decltype( std::declval< GridView >().indexSet() ) >::type IndexSet;

      typedef typename IndexSet::IndexType Index;

      MCMGMapper ( const GV &gridView, const Layout< GV::dimension > layout )
        : gridView_( gridView ),
          indexSet_( gridView_.indexSet() ),
          offset_( GlobalGeometryTypeIndex::size( GV::dimension ) ),
          layout_( layout )
      {
        update();
      }

      explicit MCMGMapper ( const GV &gridView )
        : gridView_( gridView ),
          indexSet_( gridView_.indexSet() ),
          offset_( GlobalGeometryTypeIndex::size( GV::dimension ) )
      {
        update();
      }

      template< class Entity >
      Index index ( const Entity &entity ) const
      {
        const GeometryType gt = entity.type();
        assert( layout().contains( gt ) );
        return indexSet().index( entity ) + offset_[ GlobalGeometryTypeIndex::index( gt ) ];
      }

      template< class Entity >
      Index subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        const GeometryType gt = ReferenceElements< double, Entity::mydimension >::general( entity.type() ).type( i, codim - Entity::codimension );
        assert( layout().contains( gt ) );
        return indexSet().subIndex( entity, i, codim ) + offset_[ GlobalGeometryTypeIndex::index( gt ) ];
      }

      Index size () const { return size_; }

      template< class Entity >
      bool contains ( const Entity &entity, Index &index ) const
      {
        const GeometryType gt = entity.type();
        const bool contained = (indexSet().contains( entity ) && layout().contains( gt ));
        index = (contained ? indexSet().index( entity ) + offset_[ gt ] : 0);
        return contained;
      }

      template< class Entity >
      bool contains ( const Entity &entity, int i, unsigned int codim, Index &index ) const
      {
        const GeometryType gt = ReferenceElements< double, Entity::mydimension >::general( entity.type() ).type( i, codim - Entity::codimension );
        const bool contained = (indexSet().contains( entity ) && layout().contains( gt ));
        index = (contained ? indexSet().subIndex( entity, i, codim ) + offset_[ gt ] : 0);
        return contained;
      }

      void update ()
      {
        std::fill( offset_.begin(), offset_.end(), 0 );
        size_ = 0;
        for( int codim = 0; codim <= GV::dimension; ++codim )
        {
          for( GeometryType type : indexSet().types( codim ) )
          {
            if( !layout().contains( type ) )
              continue;

            offset_[ GlobalGeometryTypeIndex::index( type ) ] = size_;
            size_ += indexSet().size( type );
          }
        }
      }

      const GridView &gridView () const { return gridView_; }
      const IndexSet &indexSet () const { return indexSet_; }

      const Layout< dimension > &layout () const { return layout_; }

    private:
      Index size_;
      GV gridView_;
      const IndexSet &indexSet_;
      std::vector< Index > offset_;
      Layout< GV::dimension > layout_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_MCMGMAPPER_HH
