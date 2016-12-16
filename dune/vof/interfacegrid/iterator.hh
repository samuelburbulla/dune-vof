#ifndef DUNE_VOF_INTERFACEGRID_ITERATOR_HH
#define DUNE_VOF_INTERFACEGRID_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entityiterator.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIterator
    // --------------

    template< class Grid, class HostIterator >
    class InterfaceGridIterator
    {
      typedef InterfaceGridIterator< Grid, HostIterator > This;

    protected:
      typedef typename remove_const< Grid >::type::Traits Traits;

    public:
      /** \brief grid dimension */
      static const int dimension = remove_const< Grid >::type::dimension;
      /** \brief world dimension */
      static const int codimension = HostIterator::codimension;

      /** \brief type of entity */
      typedef typename Traits::template Codim< codimension >::Entity Entity;

    protected:
      typedef typename Traits::ExtraData ExtraData;

      typedef typename Traits::template Codim< codimension > :: EntityImpl EntityImpl;

    public:
      InterfaceGridIterator ( ExtraData data, const HostIterator &hostIterator )
      : data_( data ),
        hostIterator_( hostIterator )
      {}

      explicit InterfaceGridIterator ( const EntityImpl &entity )
      : data_( entity.data() ),
        hostIterator_( entity.hostEntity() )
      {}

      InterfaceGridIterator ( const This &other )
      : data_( other.data_ ),
        hostIterator_( other.hostIterator_ )
      {}

      template< class HI >
      explicit InterfaceGridIterator ( const InterfaceGridIterator< Grid, HI > &other )
      : data_( other.data_ ),
        hostIterator_( other.hostIterator_ )
      {}

      const This &operator= ( const This &other )
      {
        data_ = other.data_;
        hostIterator_ = other.hostIterator_;
        return *this;
      }

      template< class HI >
      const This &operator= ( const InterfaceGridIterator< Grid, HI > &other )
      {
        data_ = other.data_;
        hostIterator_ = other.hostIterator_;
        return *this;
      }

      /** \brief check for equality */
      template< class HI >
      bool equals ( const InterfaceGridIterator< Grid, HI > &other ) const
      {
        return (hostIterator() == other.hostIterator());
      }

      /** \brief dereference entity */
      Entity dereference () const
      {
        return Entity( EntityImpl( data(), *hostIterator() ) );
      }

      /** \brief increment */
      void increment ()
      {
        ++hostIterator_;
      }

      /** \brief obtain level */
      int level () const { return hostIterator().level(); }

      /** \brief obtain host iterator */
      const HostIterator &hostIterator() const { return hostIterator_; }

    protected:
      ExtraData data () const { return data_; }

    protected:
      ExtraData data_;
      HostIterator hostIterator_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ITERATOR_HH
