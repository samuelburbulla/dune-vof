#ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
#define DUNE_VOF_INTERFACEGRID_INDEXSET_HH

#include <cstddef>

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/vof/interfacegrid/declaration.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIndexSet
    // --------------

    template< class Grid >
    class InterfaceGridIndexSet
      : public IndexSet< Grid, InterfaceGridIndexSet< Grid >, std::size_t >
    {
      typedef InterfaceGridIndexSet< Grid > This;
      typedef IndexSet< Grid, InterfaceGridIndexSet< Grid >, std::size_t > Base;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

      typedef typename Traits::Reconstruction Reconstruction;

    public:
      static const int dimension = Traits::dimension;

      typedef typename Base::IndexType IndexType;

      typedef std::array< GeometryType, 1 > Types;

      explicit InterfaceGridIndexSet ( const Flags< typename Reconstruction::GridView > &flags )
        : hostIndexSet_( flags.gridView().indexSet() ),
          indices_( hostIndexSet_.size( 0 ) )
      {
        update( flags );
      }

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      template< int cd >
      IndexType index ( const typename Traits::template Codim< cd >::Entity &entity ) const
      {
        // ...
      }

      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        return subIndex< Entity::codimension >( entity, i, codim );
      }

      template< int cd >
      IndexType subIndex ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        // ...
      }

      IndexType size ( GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        return (type == types( codim ) ? size( codim ) : static_cast< IndexType >( 0u ));
      }

      IndexType size ( int codim ) const
      {
        // ....
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        // ...
      }

      Types types( int codim ) const
      {
        assert( (codim >= 0) && (codim <= dimension) );
        const int mydim = dimension - codim;
        return {{ GeometryType( (mydim < 2 ? GeometryType::cube, GeometryType::none), mydim ) }};
      }

      void update ( const Flags< typename Reconstruction::GridView > &flags )
      {
        std::size_t size_;
        for( const auto &element : elements( flags.gridView(), Partition::all ) )
        {
          if( !flags.isMixed( element ) )
            continue;

          const std::size_t gtIndex = LocalGeometryTypeIndex::index( element.type() );
          indices_[ gtIndex ][ hostIndexSet_.index( element ) ] = size_++;
        }
      }

    private:
      const HostIndexSet &hostIndexSet_;
      IndexType size_;
      std::array< std::vector< IndexType >, LocalGeometryTypeIndex::size( dimension+1 ) > indices_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
