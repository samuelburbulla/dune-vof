#ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
#define DUNE_VOF_INTERFACEGRID_GRID_HH

#include <cstddef>

#include <dune/grid/common/grid.hh>

#include <dune/vof/interfacegrid/capabilities.hh>
#include <dune/vof/interfacegrid/datahandle.hh>
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

        typedef InterfaceGridIterator< const Grid, typename HostGrid::HierarchicIterator > HierarchicIteratorImpl;
        typedef Dune::EntityIterator< 0, const Grid, HierarchicIteratorImpl > HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, InterfaceGridGeometry > Geometry;
          typedef Dune::Geometry< dimension-codim, dimension, const Grid, InterfaceGridLocalGeometry > LocalGeometry;

          typedef Dune::Entity< codim, dimension, const Grid, InterfaceGridEntity > Entity;
          typedef Dune::EntitySeed< const Grid, InterfaceGridEntitySeed< codim, const Grid > > EntitySeed;

          template< PartitionIteratorType pitype >
          struct Partition
          {
            typedef typename HostGrid::template Codim< codim >::template Partition< pitype > HostPartition;

            typedef InterfaceGridIterator< const Grid, typename HostPartition::LeafIterator > LeafIteratorImpl;
            typedef Dune::EntityIterator< codim, const Grid, LeafIteratorImpl > LeafIterator;

            typedef InterfaceGridIterator< const Grid, typename HostPartition::LevelIterator > LevelIteratorImpl;
            typedef Dune::EntityIterator< codim, const Grid, LevelIteratorImpl > LevelIterator;
          };

          typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
          typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        };

        typedef InterfaceGridIndexSet< const Grid > LeafIndexSet;
        typedef LeafIndexSet LevelIndexSet;

        typedef InterfaceGridGlobalIdSet< const Grid > GlobalIdSet;
        typedef InterfaceGridLocalIdSet< const Grid > LocalIdSet;

        typedef typename Reconstruction::GridView::CollectiveCommunication CollectiveCommunication;

        typedef Dune::GridView< InterfaceGridViewTraits< Reconstruction > > LeafGridView;
        typedef LeafGridView LevelGridView;
      };
    };



    // InterfaceGrid
    // -------------

    template< class R >
    class InterfaceGrid
      : public GridDefaultImplementation< R::GridView::dimension-1, R::GridView::dimensionworld, typename R::GridView::ctype, InterfaceGridFamily< R > >
    {
      typedef InterfaceGrid< R > Grid;
      typedef GridDefaultImplementation< R::GridView::dimension-1, R::GridView::dimensionworld, typename R::GridView::ctype, InterfaceGridFamily< R > > Base;

      typedef InterfaceGridFamily< R >::Traits Traits;

      template< int, int, class > friend class InterfaceGridEntity;
      template< class, class > friend class InterfaceGridIntersection;
      template< class, class > friend class InterfaceGridIntersectionIterator;
      template< class, class > friend class InterfaceGridIdSet;
      template< class, class > friend class InterfaceGridIndexSet;
      template< class > friend class HostGridAccess;

    public:
      typedef R Reconstruction;

      typedef VoF::Flags< typename Reconstruction::GridView > Flags;

      typedef typename Reconstruction::ColorFunction ColorFunction;

      typedef typename Traits::LevelGridView LevelGridView;
      typedef typename Traits::LeafGridView LeafGridView;

      typedef typename Traits::LeafIndexSet LeafIndexSet;
      typedef typename Traits::LevelIndexSet LevelIndexSet;

      typedef typename Traits::GlobalIdSet GlobalIdSet;
      typedef typename Traits::LocalIdSet LocalIdSet;

      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      template< class... Args >
      explicit InterfaceGrid ( const ColorFunction &colorFunction, Args &&... args )
        : reconstruction_( std::forward< Args >( args )... ),
          flags_( reconstruction_.gridView() ),
          reconstructionSet_( reconstruction_.gridView() )
      {
        update( colorFunction );
      }

      int maxLevel () const { return 0; }

      int size ( int level, int codim ) const { return levelGridView( level ).size( codim ); }
      int size ( int codim ) const { return leafGridView().size( codim ); }

      int size ( int level, GeometryType type ) const { return levelGridView( level ).size( type ); }
      int size ( GeometryType type ) const { return leafGridView().size( type ); }

      std::size_t numBoundarySegments () const { return leafIndexSet().size( 1 ); }

      const GlobalIdSet &globalIdSet () const
      {
        if( !globalIdSet_ )
          globalIdSet_ = GlobalIdSet( hostGrid().globalIdSet() );
        return globalIdSet_;
      }

      const LocalIdSet &localIdSet () const
      {
        if( !localIdSet_ )
          localIdSet_ = LocalIdSet( hostGrid().localIdSet() );
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

      const CollectiveCommunication &comm () const { return reconstruction().gridView().comm(); }

      LevelGridView levelGridView ( int level ) const
      {
        return LevelGridView( InterfaceGridView< Reconstruction >( *this ) );
      }

      LeafGridView leafGridView () const
      {
        return LeafGridView( InterfaceGridView< Reconstruction >( *this ) );
      }

      template< class EntitySeed >
      typename Traits::template Codim< EntitySeed::codimension >::Entity entity ( const EntitySeed &seed ) const
      {
        typedef typename Traits::template Codim< EntitySeed::codimension >::EntityImpl EntityImpl;
        return EntityImpl( extraData(), hostGrid().entity( seed.impl().hostEntitySeed() ) );
      }

      const Reconstruction &reconstruction () { return reconstruction_; }
      const Flags &flags () const { return flags_; }

      void update ( const ColorFunction &colorFunction )
      {
        flags_.reflag( colorFunction, 1e-6 );
        reconstruction_( colorFunction, reconstructions_, flags_ );

        // further updates
      }

    protected:
      using Base::getRealImplementation;

      Reconstruction reconstruction_;
      typename Reconstruction::ReconstructionSet reconstructions_;
      Flags flags_;

      LeafIndexSet leafIndexSet_;
      LocalIdSet localIdSet_;
      GlobalIdSet globalIdSet_;
    };

  } // namespace VoF

} // namespace Dune

#include <dune/vof/interfacegrid/persistentcontainer.hh>

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
