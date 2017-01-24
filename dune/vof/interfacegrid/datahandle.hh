#ifndef DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH
#define DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH

#include <type_traits>

#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/vof/interfacegrid/dataset.hh>
#include <dune/vof/interfacegrid/entity.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridDataHandle
    // -----------------------

    template< class Grid, class WrappedHandle >
    class InterfaceGridDataHandle
      : public CommDataHandleIF< InterfaceGridDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef InterfaceGridDataHandle< Grid, WrappedHandle > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      typedef typename WrappedHandle::DataType DataType;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;

      static const int dimension = Traits::dimension;

      typedef Dune::Entity< 0, dimension, const Grid, InterfaceGridEntity > Element;

      InterfaceGridDataHandle ( const DataSet &dataSet, WrappedHandle &wrappedHandle )
        : dataSet_( dataSet ), wrappedHandle_( wrappedHandle )
      {}

      InterfaceGridDataHandle ( const This & ) = delete;
      InterfaceGridDataHandle ( This && ) = delete;

      InterfaceGridDataHandle &operator= ( const This & ) = delete;
      InterfaceGridDataHandle &operator= ( This && ) = delete;

      bool contains ( int dim, int codim ) const { return (codim == 0); }

      bool fixedsize ( int dim, int codim ) const { return false; }

      template< class HostEntity, std::enable_if_t< (HostEntity::codimension == 0), int > = 0 >
      std::size_t size ( const HostEntity &hostEntity ) const
      {
        Element element( InterfaceGridEntity< 0, dimension, const Grid >( dataSet(), hostEntity ) );
        std::size_t size = 0;
        Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &element, &size ] ( auto codim ) {
            if( wrappedHandle_.contains( element.dimension, codim ) )
            {
              for( unsigned int i = 0; i < element.subEntities( codim ); ++i )
                size += wrappedHandle_.size( element.template subEntity< codim >( i ) );
            }
          } );
        return size;
      }

      template< class HostEntity, std::enable_if_t< (HostEntity::codimension > 0), int > = 0 >
      std::size_t size ( const HostEntity &hostEntity ) const
      {
        return 0;
      }

      template< class MessageBuffer, class HostEntity, std::enable_if_t< (HostEntity::codimension == 0), int > = 0 >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        Element element( InterfaceGridEntity< 0, dimension, const Grid >( dataSet(), hostEntity ) );
        Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &buffer, &element ] ( auto codim ) {
            if( wrappedHandle_.contains( element.dimension, codim ) )
            {
              for( unsigned int i = 0; i < element.subEntities( codim ); ++i )
                wrappedHandle_.gather( buffer, element.template subEntity< codim >( i ) );
            }
          } );
      }

      template< class MessageBuffer, class HostEntity, std::enable_if_t< (HostEntity::codimension > 0), int > = 0 >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {}

      template< class MessageBuffer, class HostEntity, std::enable_if_t< (HostEntity::codimension == 0), int > = 0 >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, std::size_t size )
      {
        Element element( InterfaceGridEntity< 0, dimension, const Grid >( dataSet(), hostEntity ) );
        Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &buffer, &element, size ] ( auto codim ) {
            if( wrappedHandle_.contains( element.dimension, codim ) )
            {
              // Bug:  We pass the same size to all elements, which is obviously wrong.
              //       Unfortunately, this cannot be helped, because we can only transfer
              //       objects of type "DataType", i.e., we cannot additionally transfer
              //       individual sizes.
              // Note: Usually, we will only transfer element data on InterfaceGrid, so this
              //       bug is not an issue.
              for( unsigned int i = 0; i < element.subEntities( codim ); ++i )
                wrappedHandle_.scatter( buffer, element.template subEntity< codim >( i ), size );
            }
          } );
      }

      template< class MessageBuffer, class HostEntity, std::enable_if_t< (HostEntity::codimension > 0), int > = 0 >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, std::size_t size )
      {}

      const DataSet &dataSet () const { return dataSet_; }

    protected:
      const DataSet &dataSet_;
      WrappedHandle &wrappedHandle_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_DATAHANDLE_HH
