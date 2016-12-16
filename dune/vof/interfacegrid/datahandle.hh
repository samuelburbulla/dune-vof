#ifndef DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH
#define DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/datahandleif.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridDataHandle
    // ----------------

    template< class Grid, class WrappedHandle >
    class InterfaceGridDataHandle
      : public CommDataHandleIF< InterfaceGridDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef InterfaceGridDataHandle< Grid, WrappedHandle > This;

    protected:
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits :: ExtraData ExtraData ;

      template< int codim >
      struct Codim
      {
        typedef typename Traits::template Codim< codim >::Entity Entity;
      };

    public:
      // type of data to be communicated
      typedef typename WrappedHandle::DataType DataType;

      typedef CommDataHandleIF< This, DataType > DataHandleIF;

      InterfaceGridDataHandle ( ExtraData data, WrappedHandle &wrappedHandle )
        : wrappedHandle_( wrappedHandle ), data_( data )
      {}

      InterfaceGridDataHandle ( const This & ) = delete;

      bool contains ( int dim, int codim ) const
      {
        return wrappedHandle_.contains( dim, codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedsize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        typedef typename Codim< HostEntity::codimension >::Entity Entity;
        const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
        return wrappedHandle_.size( entity );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        typedef typename Codim< HostEntity::codimension >::Entity Entity;
        const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
        wrappedHandle_.gather( buffer, entity );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        typedef typename Codim< HostEntity::codimension >::Entity Entity;
        const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
        wrappedHandle_.scatter( buffer, entity, size );
      }

      ExtraData data() const { return data_; }

    protected:
      WrappedHandle &wrappedHandle_;
      ExtraData data_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH
