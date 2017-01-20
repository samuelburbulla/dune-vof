#ifndef DUNE_VOF_INTERFACEGRID_IDSET_HH
#define DUNE_VOF_INTERFACEGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIdSet
    // -----------

    template< class Grid, class HostIdSet >
    class InterfaceGridIdSet
    : public IdSet< Grid, InterfaceGridIdSet< Grid, HostIdSet >, typename HostIdSet::IdType >
    {
    protected:
      typedef InterfaceGridIdSet< Grid, HostIdSet > This;
      typedef IdSet< Grid, InterfaceGridIdSet< Grid, HostIdSet >, typename HostIdSet::IdType > Base;

      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      typedef typename HostIdSet::IdType IdType;

      InterfaceGridIdSet ()
      : hostIdSet_( nullptr )
      {}

      explicit InterfaceGridIdSet ( const HostIdSet &hostIdSet )
      : hostIdSet_( &hostIdSet )
      {}

      InterfaceGridIdSet ( const This &other )
      : hostIdSet_( other.hostIdSet_ )
      {}

      const This &operator= ( const This &other )
      {
        hostIdSet_ = other.hostIdSet_;
        return *this;
      }

      //! id meethod for entity and specific codim
      template< int codim >
      IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
      {
        return id( Grid::template getHostEntity< codim >( entity ) );
      }

      //! id method for host entity (e.g. in ParallelGrid)
      template< int codim >
      IdType id ( const typename Traits::HostGrid::template Codim< codim >::Entity &entity ) const
      {
        return hostIdSet().id( entity );
      }

      //! id method of all entities
      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id< Entity::codimension >( entity );
      }

      template< int cd >
      IdType subId ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        return subId( Grid::template getHostEntity< cd >( entity ), i, codim );
      }

      //! subId method for host entity (e.g. in ParallelGrid)
      template< int cd >
      IdType subId ( const typename Traits::HostGrid::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        return hostIdSet().subId( entity, i, codim );
      }

      //! subId method for all entities
      template< class Entity >
      IdType subId ( const Entity &entity, int i, unsigned int codim ) const
      {
        return subId< Entity::codimension >( entity, i, codim );
      }

      operator bool () const { return bool( hostIdSet_ ); }

    protected:
      const HostIdSet &hostIdSet () const
      {
        assert( *this );
        return *hostIdSet_;
      }

      const HostIdSet *hostIdSet_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_IDSET_HH
