#ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
#define DUNE_VOF_INTERFACEGRID_GRID_HH

#include <cstddef>
#include <utility>

#include <dune/grid/common/grid.hh>

#include <dune/vof/interfacegrid/capabilities.hh>
#include <dune/vof/interfacegrid/datahandle.hh>
#include <dune/vof/interfacegrid/dataset.hh>
#include <dune/vof/interfacegrid/declaration.hh>
#include <dune/vof/interfacegrid/entity.hh>
#include <dune/vof/interfacegrid/entityseed.hh>
#include <dune/vof/interfacegrid/geometry.hh>
#include <dune/vof/interfacegrid/gridview.hh>
#include <dune/vof/interfacegrid/idset.hh>

namespace Dune
{

  namespace VoF
  {

    // External Forward Declarations
    // -----------------------------

    template< class Grid >
    struct HostGridAccess;



    // InterfaceGridFamily
    // -------------------

    template< class R >
    struct InterfaceGridFamily
    {
      struct Traits
      {
        typedef InterfaceGrid< R > Grid;

        typedef R Reconstruction;

        static const int dimension = Reconstruction::GridView::dimension-1;
        static const int dimensionworld = Reconstruction::GridView::dimensionworld;

        typedef typename Reconstruction::GridView::ctype ctype;

        typedef Dune::Intersection< const Grid, InterfaceGridIntersection< const Grid > > LeafIntersection;
        typedef LeafIntersection LevelIntersection;

        typedef Dune::IntersectionIterator< const Grid, InterfaceGridIntersectionIterator< const Grid >, InterfaceGridIntersection< const Grid > > LeafIntersectionIterator;
        typedef LeafIntersectionIterator LevelIntersectionIterator;

        typedef Dune::EntityIterator< 0, const Grid, InterfaceGridHierarchicIterator< const Grid > > HierarchicIterator;

        typedef Dune::GridView< InterfaceGridViewTraits< Reconstruction > > LeafGridView;
        typedef LeafGridView LevelGridView;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, InterfaceGridGeometry > Geometry;
          typedef Dune::Geometry< dimension-codim, dimension, const Grid, InterfaceGridGeometry > LocalGeometry;

          typedef Dune::Entity< codim, dimension, const Grid, InterfaceGridEntity > Entity;
          typedef Dune::EntitySeed< const Grid, InterfaceGridEntitySeed< codim, const Grid > > EntitySeed;

          template< PartitionIteratorType pitype >
          struct Partition
          {
            typedef typename InterfaceGridViewTraits< Reconstruction >::template Codim< codim >::template Partition< pitype >::Iterator LeafIterator;
            typedef typename InterfaceGridViewTraits< Reconstruction >::template Codim< codim >::template Partition< pitype >::Iterator LevelIterator;
          };

          typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
          typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        };

        typedef InterfaceGridIndexSet< const Grid > LeafIndexSet;
        typedef LeafIndexSet LevelIndexSet;

        typedef InterfaceGridIdSet< const Grid, typename Reconstruction::GridView::Grid::GlobalIdSet > GlobalIdSet;
        typedef InterfaceGridIdSet< const Grid, typename Reconstruction::GridView::Grid::LocalIdSet > LocalIdSet;

        typedef typename Reconstruction::GridView::CollectiveCommunication CollectiveCommunication;
      };
    };



    // InterfaceGrid
    // -------------

    template< class R >
    class InterfaceGrid
      : public GridDefaultImplementation< R::GridView::dimension-1, R::GridView::dimensionworld, typename R::GridView::ctype, InterfaceGridFamily< R > >
    {
      typedef InterfaceGrid< R > This;
      typedef GridDefaultImplementation< R::GridView::dimension-1, R::GridView::dimensionworld, typename R::GridView::ctype, InterfaceGridFamily< R > > Base;

      template< int, int, class > friend class InterfaceGridEntity;
      template< class > friend class InterfaceGridIntersection;
      template< class, class > friend class InterfaceGridIdSet;
      template< class > friend class InterfaceGridIndexSet;
      friend struct HostGridAccess< This >;

    public:
      typedef InterfaceGridFamily< R > GridFamily;

      typedef typename GridFamily::Traits Traits;

      typedef R Reconstruction;

      typedef InterfaceGridDataSet< Reconstruction > DataSet;

      typedef typename DataSet::Flags Flags;

      static const int dimension = Traits::dimension;

      typedef typename Traits::LevelGridView LevelGridView;
      typedef typename Traits::LeafGridView LeafGridView;

      typedef typename Traits::LeafIndexSet LeafIndexSet;
      typedef typename Traits::LevelIndexSet LevelIndexSet;

      typedef typename Traits::GlobalIdSet GlobalIdSet;
      typedef typename Traits::LocalIdSet LocalIdSet;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< class ColorFunction, class... Args >
      explicit InterfaceGrid ( const ColorFunction &colorFunction, Args &&... args )
        : leafIndexSet_( colorFunction, std::forward< Args >( args )... )
      {}

      int maxLevel () const { return 0; }

      int size ( int level, int codim ) const { return levelGridView( level ).size( codim ); }
      int size ( int codim ) const { return leafGridView().size( codim ); }

      int size ( int level, GeometryType type ) const { return levelGridView( level ).size( type ); }
      int size ( GeometryType type ) const { return leafGridView().size( type ); }

      std::size_t numBoundarySegments () const { return dataSet().offsets()[ dataSet().indices().size() ]; }

      const GlobalIdSet &globalIdSet () const
      {
        if( !globalIdSet_ )
          globalIdSet_ = GlobalIdSet( dataSet().gridView().grid().globalIdSet() );
        return globalIdSet_;
      }

      const LocalIdSet &localIdSet () const
      {
        if( !localIdSet_ )
          localIdSet_ = LocalIdSet( dataSet().gridView().grid().localIdSet() );
        return localIdSet_;
      }

      const LevelIndexSet &levelIndexSet ( int level ) const { assert( level == 0 ); return leafIndexSet(); }
      const LeafIndexSet &leafIndexSet () const { return leafIndexSet_; }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle, InterfaceType interface, CommunicationDirection direction, int level ) const
      {
        levelGridView( level ).communicate( dataHandle, interface, direction );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle, InterfaceType interface, CommunicationDirection direction ) const
      {
        leafGridView().communicate( dataHandle, interface, direction );
      }

      const CollectiveCommunication &comm () const { return dataSet().gridView().comm(); }

      LevelGridView levelGridView ( int level ) const
      {
        return LevelGridView( InterfaceGridView< Reconstruction >( *this ) );
      }

      LeafGridView leafGridView () const
      {
        return LeafGridView( InterfaceGridView< Reconstruction >( *this ) );
      }

      template< class EntitySeed, std::enable_if_t< (EntitySeed::codimension == 0), int > = 0 >
      typename Traits::template Codim< EntitySeed::codimension >::Entity entity ( const EntitySeed &seed ) const
      {
        typedef InterfaceGridEntity< EntitySeed::codimension, dimension, const This > Impl;
        return Impl( dataSet(), dataSet().gridView().grid().entity( getRealImplementation( seed ).hostElementSeed() ) );
      }

      template< class EntitySeed, std::enable_if_t< (EntitySeed::codimension > 0), int > = 0 >
      typename Traits::template Codim< EntitySeed::codimension >::Entity entity ( const EntitySeed &seed ) const
      {
        typedef InterfaceGridEntity< EntitySeed::codimension, dimension, const This > Impl;
        return Impl( dataSet(), dataSet().gridView().grid().entity( getRealImplementation( seed ).hostElementSeed() ), getRealImplementation( seed ).subEntity() );
      }

      const DataSet &dataSet () const { return leafIndexSet().dataSet(); }

      const Reconstruction &reconstruction () const { return dataSet().reconstruction(); }
      const Flags &flags () const { return dataSet().flags(); }

      template< class ColorFunction >
      void update ( const ColorFunction &colorFunction ) { leafIndexSet_.update( colorFunction ); }

    protected:
      using Base::getRealImplementation;

      LeafIndexSet leafIndexSet_;
      mutable LocalIdSet localIdSet_;
      mutable GlobalIdSet globalIdSet_;
    };



    // interfaceGrid
    // -------------
    
    template< class ColorFunction, class Reconstruction >
    inline static auto interfaceGrid ( const ColorFunction &color, Reconstruction reconstruction ) 
      -> InterfaceGrid< Reconstruction >
    {
      return InterfaceGrid< Reconstruction >( color, std::move( reconstruction ) );
    }

  } // namespace VoF

} // namespace Dune

#include <dune/vof/interfacegrid/persistentcontainer.hh>

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
