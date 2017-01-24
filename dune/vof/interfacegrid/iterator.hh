#ifndef DUNE_VOF_INTERFACEGRID_ITERATOR_HH
#define DUNE_VOF_INTERFACEGRID_ITERATOR_HH

#include <cassert>

#include <type_traits>

#include <dune/grid/common/entityiterator.hh>

#include <dune/vof/interfacegrid/dataset.hh>
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

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = Traits::dimension;

      typedef Dune::Entity< codimension, dimension, Grid, InterfaceGridEntity > Entity;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;

      InterfaceGridIterator () = default;

      InterfaceGridIterator ( const DataSet &dataSet, const HostElementIterator &hostElementBegin, const HostElementIterator &hostElementEnd )
        : dataSet_( &dataSet ), hostElementIterator_( hostElementBegin ), hostElementEnd_( hostElementEnd )
      {
        for( ; (hostElementIterator_ != hostElementEnd_) && !dataSet.flags().isMixed( *hostElementIterator_ ); ++hostElementIterator_ )
          continue;
      }

      InterfaceGridIterator ( const DataSet &dataSet, const HostElementIterator &hostElementEnd )
        : dataSet_( &dataSet ), hostElementIterator_( hostElementEnd ), hostElementEnd_( hostElementEnd )
      {}

      operator bool () const { return dataSet_ && (hostElementIterator_ != hostElementEnd_); }

      bool equals ( const This &other ) const { return (hostElementIterator_ == other.hostElementIterator_); }

      Entity dereference () const { return InterfaceGridEntity< codimension, dimension, Grid >( dataSet(), *hostElementIterator() ); }

      void increment ()
      {
        for( ++hostElementIterator_; (hostElementIterator_ != hostElementEnd_) && !dataSet().flags().isMixed( *hostElementIterator_ ); ++hostElementIterator_ )
          continue;
      }

      const DataSet &dataSet () const { assert( dataSet_ ); return *dataSet_; }
      const HostElementIterator &hostElementIterator() const { return hostElementIterator_; }

    protected:
      const DataSet *dataSet_;
      HostElementIterator hostElementIterator_, hostElementEnd_;
    };



    // InterfaceGridIterator
    // ---------------------

    template< int codim, class Grid, class HostElementIterator >
    class InterfaceGridIterator
    {
      typedef InterfaceGridIterator< codim, Grid, HostElementIterator > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

      typedef InterfaceGridIterator< 0, Grid, HostElementIterator > ElementIterator;

    public:
      static const int codimension = codim;
      static const int dimension = Traits::dimension;

      typedef Dune::Entity< codimension, dimension, Grid, InterfaceGridEntity > Entity;

      typedef typename ElementIterator::DataSet DataSet;

      InterfaceGridIterator () = default;

      InterfaceGridIterator ( const DataSet &dataSet, const HostElementIterator &hostElementBegin, const HostElementIterator &hostElementEnd )
        : elementIterator_( dataSet, hostElementBegin, hostElementEnd ),
          subEntities_(*this ? dataSet.numVertices( *hostElementIterator() ) : 0)
      {}

      InterfaceGridIterator ( const DataSet &dataSet, const HostElementIterator &hostElementEnd )
        : elementIterator_( dataSet, hostElementEnd )
      {}

      operator bool () const { return static_cast< bool >( elementIterator_ ); }

      bool equals ( const This &other ) const { return elementIterator_.equals( other.elementIterator_) && (subEntity_ == other.subEntity_); }

      Entity dereference () const { return InterfaceGridEntity< codimension, dimension, Grid >( dataSet(), *hostElementIterator(), subEntity_ ); }

      void increment ()
      {
        if( ++subEntity_ < subEntities_ )
          return;

        elementIterator_.increment();
        subEntity_ = 0;
        subEntities_ = (*this ? dataSet().numVertices( *hostElementIterator() ) : 0);
      }

      const DataSet &dataSet () const { return elementIterator_.dataSet(); }
      const HostElementIterator &hostElementIterator() const { return elementIterator_.hostElementIterator(); }

    private:
      ElementIterator elementIterator_;
      int subEntity_ = 0, subEntities_ = 0;
    };



    // InterfaceGridHierarchicIterator
    // -------------------------------

    template< class Grid >
    class InterfaceGridHierarchicIterator
    {
      typedef InterfaceGridHierarchicIterator< Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = Traits::dimension;

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
