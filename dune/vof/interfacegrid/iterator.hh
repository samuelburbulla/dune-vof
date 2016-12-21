#ifndef DUNE_VOF_INTERFACEGRID_ITERATOR_HH
#define DUNE_VOF_INTERFACEGRID_ITERATOR_HH

#include <cassert>

#include <type_traits>

#include <dune/grid/common/entityiterator.hh>

#include <dune/vof/flags.hh>
#include <dune/vof/interfacegrid/entity.hh>

namespace Dune
{

  namespace VoF
  {

    // Internal Forward Declarations
    // -----------------------------

    template< int codim, class Grid, class HostElementIterator >
    class InterfaceGridIterator;



    // InterfaceGridIterator
    // ---------------------

    template< class Grid, class HostElementIterator >
    class InterfaceGridIterator< 0, Grid, HostElementIterator >
    {
      typedef InterfaceGridIterator< 0, Grid, HostElementIterator > This;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = std::remove_const< Grid >::type::dimension;

      typedef Dune::Entity< codimension, dimension, Grid, InterfaceGridEntity > Entity;

      typedef VoF::Flags< typename Traits::Reconstruction::GridView > Flags;

      InterfaceGridIterator () = default;

      InterfaceGridIterator ( const Flags &flags, const HostElementIterator &hostElementBegin, const HostElementIterator &hostElementEnd )
        : flags_( &flags ), hostElementIterator_( hostElementBegin ), hostElementEnd_( hostElementEnd )
      {
        for( ; (hostElementIterator_ != hostElementEnd_) && !flags.isMixed( *hostElementIterator_ ); ++hostElementIterator_ )
          continue;
      }

      InterfaceGridIterator ( const Flags &flags, const HostElementIterator &hostElementEnd )
        : flags_( &flags ), hostElementIterator_( hostElementEnd ), hostElementEnd_( hostElementEnd )
      {}

      operator bool () const { return flags_ && (hostElementIterator_ != hostElementEnd_); }

      bool equals ( const This &other ) const { return (hostElementIterator_ == other.hostElementIterator_); }

      Entity dereference () const { return InterfaceGridEntity< codimension, dimension, Grid >( *hostElementIterator() ); }

      void increment ()
      {
        for( ++hostElementIterator_; (hostElementIterator_ != hostElementEnd_) && !flags().isMixed( *hostElementIterator_ ); ++hostElementIterator_ )
          continue;
      }

      const Flags &flags () const { assert( flags_ ); return *flags_; }
      const HostIterator &hostElementIterator() const { return hostElementIterator_; }

    protected:
      const Flags *flags_;
      HostIterator hostElementIterator_, hostElementEnd_;
    };



    // InterfaceGridIterator
    // ---------------------

    template< int codim, class Grid, class HostElementIterator >
    class InterfaceGridIterator
    {
      typedef InterfaceGridIterator< codim, Grid, HostIterator > This;

      typedef typename std::remove_const< Grid >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = std::remove_const< Grid >::type::dimension;

      typedef Dune::Entity< codimension, dimension, Grid, InterfaceGridEntity > Entity;

      typedef VoF::Flags< typename Traits::Reconstruction::GridView > Flags;

      InterfaceGridIterator () = default;

      InterfaceGridIterator ( const Flags &flags, const HostElementIterator &hostElementBegin, const HostElementIterator &hostElementEnd )
        : elementIterator_( flags, hostElementBegin, hostElementEnd )
      {}

      InterfaceGridIterator ( const Flags &flags, const HostElementIterator &hostElementEnd )
        : elementIterator_( flags, hostElementEnd )
      {}

      operator bool () const { return static_cast< bool >( elementIterator_ ); }

      bool equals ( const This &other ) const { return elementIterator_.equals( other.elementIterator_) && (subEntity_ == other.subEntity_); }

      Entity dereference () const { return InterfaceGridEntity< codimension, dimension, Grid >( *hostElementIterator(), subEntity_ ); }

      void increment ()
      {
        // TODO: extend to (dimension == 2), i.e., the 3d case
        if( ++subEntity < 2 )
          return;

        subEntity = 0;
        elementIterator_.increment();
      }

      const Flags &flags () const { return elementIterator_.flags(); }
      const HostElementIterator &hostElementIterator() const { return elementIterator_.hostElementIterator(); }

    private:
      ElementIterator elementIterator_;
      int subEntity_ = 0;
    };



    // InterfaceGridHierarchicIterator
    // -------------------------------

    template< class Grid >
    class InterfaceGridHierarchicIterator
    {
      typedef InterfaceGridHierarchicIterator< Grid > This;

    public:
      static const int codimension = 0;
      static const int dimension = std::remove_const< Grid >::type::dimension;

      typedef Dune::Entity< codimension, dimension, Grid, InterfaceGridEntity > Entity;

      InterfaceGridHierarchicIterator () = default;

      operator bool () const { return false; }

      bool equals ( const This &other ) const { return true; }

      Entity dereference () const { DUNE_THROW( GridError, "InterfaceGrid consists of only one level" ); }

      void increment () { DUNE_THROW( GridError, "Cannot increment beyond the end iterator" ); }
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ITERATOR_HH
