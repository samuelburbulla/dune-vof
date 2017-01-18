#ifndef DUNE_VOF_INTERFACEGRID_GRIDVIEW_HH
#define DUNE_VOF_INTERFACEGRID_GRIDVIEW_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/vof/interfacegrid/datahandle.hh>
#include <dune/vof/interfacegrid/entity.hh>
#include <dune/vof/interfacegrid/indexset.hh>
#include <dune/vof/interfacegrid/intersection.hh>
#include <dune/vof/interfacegrid/intersectioniterator.hh>
#include <dune/vof/interfacegrid/iterator.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridView
    // -----------------

    template< class Reconstruction >
    class InterfaceGridView;

    

    // InterfaceGridViewTraits
    // -----------------------

    template< class Reconstruction >
    struct InterfaceGridViewTraits
    {
      typedef InterfaceGridView< Reconstruction > GridViewImp;

      typedef InterfaceGrid< Reconstruction > Grid;

      static const int dimension = Reconstruction::GridView::dimension-1;
      static const int dimensionworld = Reconstruction::GridView::dimensionworld;

      typedef InterfaceGridIndexSet< const Grid > IndexSet;

      typedef Dune::Intersection< const Grid, InterfaceGridIntersection< const Grid > > Intersection;

      typedef Dune::IntersectionIterator< const Grid, InterfaceGridIntersectionIterator< const Grid >, InterfaceGridIntersection< const Grid > > IntersectionIterator;

      typedef typename Reconstruction::GridView::CollectiveCommunication CollectiveCommunication;

      template< PartitionIteratorType pitype >
      using HostIterator = typename Reconstruction::GridView::template Codim< 0 >::template Partition< pitype >::Iterator;

      template< int codim >
      struct Codim
      {
        typedef Dune::Entity< codim, dimension, const Grid, InterfaceGridEntity > Entity;

        typedef Dune::Geometry< dimension - codim, dimensionworld, const Grid, InterfaceGridGeometry > Geometry;
        typedef Dune::Geometry< dimension - codim, dimension, const Grid, InterfaceGridGeometry > LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::EntityIterator< codim, const Grid, InterfaceGridIterator< codim, const Grid, HostIterator< pitype > > > Iterator;
        };

        typedef typename Partition< All_Partition >::Iterator Iterator;
      };

      static const bool conforming = true;
    };



    // InterfaceGridView
    // -----------------

    template< class Reconstruction >
    class InterfaceGridView
    {
      typedef InterfaceGridView< Reconstruction > This;

      typedef InterfaceGridViewTraits< Reconstruction > Traits;

    public:
      typedef typename Traits::Grid Grid;

      typedef typename Traits::IndexSet IndexSet;

      typedef typename Traits::Intersection Intersection;
      typedef typename Traits::IntersectionIterator IntersectionIterator;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< int codim >
      struct Codim
        : public Traits::template Codim< codim >
      {};

      static const bool conforming = Traits::conforming;

      explicit InterfaceGridView ( const Grid &grid ) : grid_( &grid ) {}

      const Grid &grid () const { assert( grid_ ); return *grid_; }

      const IndexSet &indexSet () const { return grid().leafIndexSet(); }
 
      int size ( int codim ) const { return indexSet().size( codim ); }

      int size ( const GeometryType &type ) const { return indexSet().size( type ); }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator begin () const
      {
        typedef InterfaceGridIterator< codim, const Grid, typename Traits::template HostIterator< pitype > > Impl;
        return Impl( grid().dataSet(), hostGridView().template begin< 0, pitype >(), hostGridView().template end< 0, pitype >() );
      }

      template< int codim >
      typename Codim< codim >::Iterator begin () const
      {
        return begin< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::Iterator end () const
      {
        typedef InterfaceGridIterator< codim, const Grid, typename Traits::template HostIterator< pitype > > Impl;
        return Impl( grid().dataSet(), hostGridView().template end< 0, pitype >() );
      }

      template< int codim >
      typename Codim< codim >::Iterator end () const
      {
        return end< codim, All_Partition >();
      }

      IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
      {
        return InterfaceGridIntersectionIterator< const Grid >( entity, 0 );
      }

      IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
      {
        return InterfaceGridIntersectionIterator< const Grid >( entity, entity.subEntities( 1 ) );
      }

      const CollectiveCommunication &comm () const { return grid().comm(); }

      int overlapSize ( int codim ) const { return hostGridView().overlapSize( codim ); }
      int ghostSize ( int codim ) const { return hostGridView().ghostSize( codim ); }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle, InterfaceType interface, CommunicationDirection direction ) const
      {
        // TODO: Please implement me
      }

    protected:
      decltype( auto ) hostGridView () const { return grid().dataSet().gridView(); }

      const Grid *grid_ = nullptr;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GRIDVIEW_HH
