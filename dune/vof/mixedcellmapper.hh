#ifndef DUNE_VOF_MIXEDCELLMAPPER_HH
#define DUNE_VOF_MIXEDCELLMAPPER_HH

#include <cassert>
#include <cstddef>

#include <limits>
#include <vector>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{

  namespace VoF
  {

    // MixedCellMapper
    // ---------------

    template< class GridView >
    class MixedCellMapper
    {
      typedef MultipleCodimMultipleGeomTypeMapper< GridView, MCMGElementLayout > ElementMapper;

      typedef typename ElementMapper::Index ElementIndex;

    public:
      static const int dimension = GridView::dimension;

      typedef std::size_t Index;

      explicit MixedCellMapper ( const Flags< GridView > &flags )
        : elementMapper_( flags.gridView() ), indices_( elementMapper_.size() )
      {
        update( flags );
      }

      template< class Entity >
      Index index ( const Entity &entity ) const
      {
        const ElementIndex elementIndex = elementMapper_.index( entity );
        assert( indices_[ elementMapper_.index( entity ) ] != invalidIndex() );
        return indices_[ elementMapper_.index( entity ) ];
      }

      template< class Entity >
      Index subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        const ElementIndex elementIndex = elementMapper_.subIndex( entity, i, codim );
        assert( indices_[ elementIndex ] != invalidIndex() );
        return indices_[ elementIndex ];
      }

      Index size () const { return size_; }

      template< class Entity >
      bool contains ( const Entity &entity, Index &index ) const
      {
        ElementIndex elementIndex;
        const bool contains = elementMapper_.contains( entity, elementIndex );
        index = (contains ? indices_[ elementIndex ] : invalidIndex());
        return contains;
      }

      template< class Entity >
      bool contains ( const Entity &entity, int i, unsigned int codim, Index &index ) const
      {
        ElementIndex elementIndex;
        const bool contains = elementMapper_.contains( entity, i, codim, elementIndex );
        index = (contains ? indices_[ elementIndex ] : invalidIndex());
        return contains;
      }

      void update ( const Flags< GridView > &flags )
      {
        elementMapper_.update();

        size_ = 0u;
        for( const auto &element : elements( flags.gridView(), Partitions::all ) )
          indices_[ elementMapper_.index( element ) ] = (flags.isMixed( element ) ? size_++ : invalidIndex());
      }

    private:
      static constexpr Index invalidIndex () noexcept { return std::numeric_limits< Index >::max(); }

      ElementMapper elementMapper_;
      std::vector< Index > indices_;
      Index size_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_MIXEDCELLMAPPER_HH
