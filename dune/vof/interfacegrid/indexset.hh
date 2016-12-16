#ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
#define DUNE_VOF_INTERFACEGRID_INDEXSET_HH

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

    template< class Grid, class HostIndexSet >
    class InterfaceGridIndexSet
    : public IndexSet< Grid, InterfaceGridIndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType >
    {
    protected:
      typedef InterfaceGridIndexSet< Grid, HostIndexSet > This;
      typedef IndexSet< Grid, InterfaceGridIndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType > Base;

      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::HostGrid HostGrid;

    public:
      static const int dimension = Traits::dimension;

      typedef typename Base::IndexType IndexType;
      typedef typename HostIndexSet::Types Types;

      InterfaceGridIndexSet ()
      : hostIndexSet_( nullptr )
      {}

      explicit InterfaceGridIndexSet ( const HostIndexSet &hostIndexSet )
      : hostIndexSet_( &hostIndexSet )
      {}

      InterfaceGridIndexSet ( const This &other )
      : hostIndexSet_( other.hostIndexSet_ )
      {}

      const This &operator= ( const This &other )
      {
        hostIndexSet_ = other.hostIndexSet_;
        return *this;
      }

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      template< int cd >
      IndexType index ( const typename Traits::template Codim< cd >::Entity &entity ) const
      {
        return index( Grid::template getHostEntity< cd >( entity ) );
      }

      template< int cd >
      IndexType index ( const typename Traits::HostGrid::template Codim< cd >::Entity &entity ) const
      {
        return hostIndexSet().index( entity );
      }

      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        return subIndex< Entity::codimension >( entity, i, codim );
      }

      template< int cd >
      IndexType subIndex ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        return subIndex( Grid::template getHostEntity< cd >( entity ), i, codim );
      }

      template< int cd >
      IndexType subIndex ( const typename Traits::HostGrid::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        return hostIndexSet().subIndex( entity, i, codim );
      }

      IndexType size ( GeometryType type ) const
      {
        return hostIndexSet().size( type );
      }

      int size ( int codim ) const
      {
        return hostIndexSet().size( codim );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        static const int cc = Entity::codimension;
        return hostIndexSet().contains( Grid::template getHostEntity< cc >( entity ) );
      }

      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

      Types types( int codim ) const
      {
        return hostIndexSet().types( codim );
      }

      operator bool () const { return bool( hostIndexSet_ ); }

    protected:
      const HostIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ );
        return *hostIndexSet_;
      }

      const HostIndexSet *hostIndexSet_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
