#ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
#define DUNE_VOF_INTERFACEGRID_GRID_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/common/hostgridinoutstreams.hh>
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

    // InterfaceGridExportParams
    // ------------------

    template< class HG >
    struct InterfaceGridExportParams
    : public HostGridHasInOutStreams< HG, Conversion< HG, HasObjectStream >::exists >
    {
      typedef HG HostGrid;
    };



    // InterfaceGridFamily
    // ------------

    template< class HostGrid >
    struct InterfaceGridFamily
    {
      struct Traits
      : public InterfaceGridExportParams< HostGrid >
      {
        typedef InterfaceGrid< HostGrid > Grid;

        typedef typename HostGrid::ctype ctype;

        struct EmptyData{};

        // type of data passed to entities, intersections, and iterators
        // for InterfaceGrid this is just an empty place holder
        typedef EmptyData ExtraData;

        static const int dimension = HostGrid::dimension;
        static const int dimensionworld = HostGrid::dimensionworld;

        typedef InterfaceGridIntersection< const Grid, typename HostGrid::LeafIntersection > LeafIntersectionImpl;
        typedef Dune::Intersection< const Grid, LeafIntersectionImpl > LeafIntersection;

        typedef InterfaceGridIntersection< const Grid, typename HostGrid::LevelIntersection > LevelIntersectionImpl;
        typedef Dune::Intersection< const Grid, LevelIntersectionImpl > LevelIntersection;

        typedef InterfaceGridIntersectionIterator< const Grid, typename HostGrid::LeafIntersectionIterator > LeafIntersectionIteratorImpl;
        typedef Dune::IntersectionIterator< const Grid, LeafIntersectionIteratorImpl, LeafIntersectionImpl > LeafIntersectionIterator;

        typedef InterfaceGridIntersectionIterator< const Grid, typename HostGrid::LevelIntersectionIterator > LevelIntersectionIteratorImpl;
        typedef Dune::IntersectionIterator< const Grid, LevelIntersectionIteratorImpl, LevelIntersectionImpl > LevelIntersectionIterator;

        typedef InterfaceGridIterator< const Grid, typename HostGrid::HierarchicIterator > HierarchicIteratorImpl;
        typedef Dune::EntityIterator< 0, const Grid, HierarchicIteratorImpl > HierarchicIterator;

        template< int codim >
        struct Codim
        {
          typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, InterfaceGridGeometry > Geometry;
          typedef Dune::Geometry< dimension-codim, dimension, const Grid, InterfaceGridLocalGeometry > LocalGeometry;

          typedef InterfaceGridEntity< codim, dimension, const Grid > EntityImpl;
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

        typedef InterfaceGridIndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet > LeafIndexSet;
        typedef InterfaceGridIndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet > LevelIndexSet;

        typedef InterfaceGridIdSet< const Grid, typename HostGrid::Traits::GlobalIdSet > GlobalIdSet;
        typedef InterfaceGridIdSet< const Grid, typename HostGrid::Traits::LocalIdSet > LocalIdSet;

        typedef typename HostGrid::Traits::CollectiveCommunication CollectiveCommunication;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::GridView< InterfaceGridViewTraits< typename HostGrid::LeafGridView,  pitype > >  LeafGridView;
          typedef Dune::GridView< InterfaceGridViewTraits< typename HostGrid::LevelGridView, pitype > > LevelGridView;
        };

        typedef typename Partition< All_Partition > :: LeafGridView   LeafGridView;
        typedef typename Partition< All_Partition > :: LevelGridView  LevelGridView;
      };
    };



    // InterfaceGrid
    // ------

    /** \class InterfaceGrid
     *  \brief identical grid wrapper
     *  \ingroup InterfaceGrid
     *
     *  \tparam  HostGrid   DUNE grid to be wrapped (called host grid)
     *
     *  \nosubgrouping
     */
    template< class HostGrid >
    class InterfaceGrid
    /** \cond */
    : public GridDefaultImplementation
        < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, InterfaceGridFamily< HostGrid > >,
      public InterfaceGridExportParams< HostGrid >
    /** \endcond */
    {
      typedef InterfaceGrid< HostGrid > Grid;

      typedef GridDefaultImplementation
        < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, InterfaceGridFamily< HostGrid > >
        Base;

      template< int, int, class > friend class InterfaceGridEntity;
      template< class, class > friend class InterfaceGridIntersection;
      template< class, class > friend class InterfaceGridIntersectionIterator;
      template< class, class > friend class InterfaceGridIdSet;
      template< class, class > friend class InterfaceGridIndexSet;
      template< class > friend class HostGridAccess;

    public:
      /** \cond */
      typedef InterfaceGridFamily< HostGrid > GridFamily;
      /** \endcond */

      /** \name Traits
       *  \{ */

      //! type of the grid traits
      typedef typename GridFamily::Traits Traits;

      /** \brief traits structure containing types for a codimension
       *
       *  \tparam codim  codimension
       *
       *  \nosubgrouping
       */
      template< int codim >
      struct Codim;

      /** \} */

      /** \name Iterator Types
       *  \{ */

      //! iterator over the grid hierarchy
      typedef typename Traits::HierarchicIterator HierarchicIterator;
      //! iterator over intersections with other entities on the leaf level
      typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
      //! iterator over intersections with other entities on the same level
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      /** \} */

      /** \name Grid View Types
       *  \{ */

      /** \brief Types for GridView */
      typedef typename Traits::LevelGridView LevelGridView;
      typedef typename Traits::LeafGridView LeafGridView;

      /** \} */

      /** \name Index and Id Set Types
       *  \{ */

      /** \brief type of leaf index set
       *
       *  The index set assigns consecutive indices to the entities of the
       *  leaf grid. The indices are of integral type and can be used to access
       *  arrays.
       *
       *  The leaf index set is a model of Dune::IndexSet.
       */
      typedef typename Traits::LeafIndexSet LeafIndexSet;

      /** \brief type of level index set
       *
       *  The index set assigns consecutive indices to the entities of a grid
       *  level. The indices are of integral type and can be used to access
       *  arrays.
       *
       *  The level index set is a model of Dune::IndexSet.
       */
      typedef typename Traits::LevelIndexSet LevelIndexSet;

      /** \brief type of global id set
       *
       *  The id set assigns a unique identifier to each entity within the
       *  grid. This identifier is unique over all processes sharing this grid.
       *
       *  \note Id's are neither consecutive nor necessarily of an integral
       *        type.
       *
       *  The global id set is a model of Dune::IdSet.
       */
      typedef typename Traits::GlobalIdSet GlobalIdSet;

      /** \brief type of local id set
       *
       *  The id set assigns a unique identifier to each entity within the
       *  grid. This identifier needs only to be unique over this process.
       *
       *  Though the local id set may be identical to the global id set, it is
       *  often implemented more efficiently.
       *
       *  \note Ids are neither consecutive nor necessarily of an integral
       *        type.
       *  \note Local ids need not be compatible with global ids. Also, no
       *        mapping from local ids to global ones needs to exist.
       *
       *  The global id set is a model of Dune::IdSet.
       */
      typedef typename Traits::LocalIdSet LocalIdSet;

      /** \} */

      /** \name Miscellaneous Types
       * \{ */

      //! type of vector coordinates (e.g., double)
      typedef typename Traits::ctype ctype;

      //! communicator with all other processes having some part of the grid
      typedef typename Traits::CollectiveCommunication CollectiveCommunication;

      /** \} */

      /** \name Construction and Destruction
       *  \{ */

      /** \brief constructor
       *
       *  The references to host grid and coordinate function are stored in the
       *  grid. Therefore, they must remain valid until the grid is destroyed.
       *
       *  \param[in]  hostGrid       reference to the grid to wrap
       */
      explicit InterfaceGrid ( HostGrid &hostGrid )
      : hostGrid_( &hostGrid ),
        levelIndexSets_( hostGrid.maxLevel()+1, nullptr )
      {}

      /** \brief destructor
       */
      ~InterfaceGrid ()
      {
        for( unsigned int i = 0; i < levelIndexSets_.size(); ++i )
        {
          if( levelIndexSets_[ i ] )
            delete( levelIndexSets_[ i ] );
        }
      }

      /** \} */

      /** \name Size Methods
       *  \{ */

      /** \brief obtain maximal grid level
       *
       *  Grid levels are numbered 0, ..., L, where L is the value returned by
       *  this method.
       *
       *  \returns maximal grid level
       */
      int maxLevel () const
      {
        return hostGrid().maxLevel();
      }

      /** \brief obtain number of entites on a level
       *
       *  \param[in]  level  level to consider
       *  \param[in]  codim  codimension to consider
       *
       *  \returns number of entities of codimension \em codim on grid level
       *           \em level.
       */
      int size ( int level, int codim ) const
      {
        return hostGrid().size( level, codim );
      }

      /** \brief obtain number of leaf entities
       *
       *  \param[in]  codim  codimension to consider
       *
       *  \returns number of leaf entities of codimension \em codim
       */
      int size ( int codim ) const
      {
        return hostGrid().size( codim );
      }

      /** \brief obtain number of entites on a level
       *
       *  \param[in]  level  level to consider
       *  \param[in]  type   geometry type to consider
       *
       *  \returns number of entities with a geometry of type \em type on grid
       *           level \em level.
       */
      int size ( int level, GeometryType type ) const
      {
        return hostGrid().size( level, type );
      }

      /** \brief returns the number of boundary segments within the macro grid
       *
       *  \returns number of boundary segments within the macro grid
       */
      int size ( GeometryType type ) const
      {
        return hostGrid().size( type );
      }

      /** \brief obtain number of leaf entities
       *
       *  \param[in]  type   geometry type to consider
       *
       *  \returns number of leaf entities with a geometry of type \em type
       */
      size_t numBoundarySegments () const
      {
        return hostGrid().numBoundarySegments( );
      }
      /** \} */

      template< int codim >
      typename Codim< codim >::LevelIterator lbegin ( int level ) const
      {
        return lbegin< codim, All_Partition >( level );
      }

      template< int codim >
      typename Codim< codim >::LevelIterator lend ( int level ) const
      {
        return lend< codim, All_Partition >( level );
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::LevelIterator
      lbegin ( int level ) const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::LevelIteratorImpl Impl;
        return Impl( extraData(), hostGrid().template lbegin< codim, pitype >( level ) );
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::LevelIterator
      lend ( int level ) const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::LevelIteratorImpl Impl;
        return Impl( extraData(), hostGrid().template lend< codim, pitype >( level ) );
      }

      template< int codim >
      typename Codim< codim >::LeafIterator leafbegin () const
      {
        return leafbegin< codim, All_Partition >();
      }

      template< int codim >
      typename Codim< codim >::LeafIterator leafend () const
      {
        return leafend< codim, All_Partition >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::LeafIterator
      leafbegin () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::LeafIteratorImpl Impl;
        return Impl( extraData(), hostGrid().template leafbegin< codim, pitype >() );
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::LeafIterator
      leafend () const
      {
        typedef typename Traits::template Codim< codim >::template Partition< pitype >::LeafIteratorImpl Impl;
        return Impl( extraData(), hostGrid().template leafend< codim, pitype >() );
      }

      const GlobalIdSet &globalIdSet () const
      {
        if( !globalIdSet_ )
          globalIdSet_ = GlobalIdSet( hostGrid().globalIdSet() );
        assert( globalIdSet_ );
        return globalIdSet_;
      }

      const LocalIdSet &localIdSet () const
      {
        if( !localIdSet_ )
          localIdSet_ = LocalIdSet( hostGrid().localIdSet() );
        assert( localIdSet_ );
        return localIdSet_;
      }

      const LevelIndexSet &levelIndexSet ( int level ) const
      {
        assert( levelIndexSets_.size() == (size_t)(maxLevel()+1) );
        if( (level < 0) || (level > maxLevel()) )
        {
          DUNE_THROW( GridError, "LevelIndexSet for nonexisting level " << level
                                 << " requested." );
        }

        LevelIndexSet *&levelIndexSet = levelIndexSets_[ level ];
        if( !levelIndexSet )
          levelIndexSet = new LevelIndexSet( hostGrid().levelIndexSet( level ) );
        assert( levelIndexSet );
        return *levelIndexSet;
      }

      const LeafIndexSet &leafIndexSet () const
      {
        if( !leafIndexSet_ )
          leafIndexSet_ = LeafIndexSet( hostGrid().leafIndexSet() );
        assert( leafIndexSet_ );
        return leafIndexSet_;
      }

      void globalRefine ( int refCount )
      {
        hostGrid().globalRefine( refCount );
        // update overall status
        update();
      }

      bool mark ( int refCount, const typename Codim< 0 >::Entity &entity )
      {
        return hostGrid().mark( refCount, getHostEntity< 0 >( entity ) );
      }

      int getMark ( const typename Codim< 0 >::Entity &entity ) const
      {
        return hostGrid().getMark( getHostEntity< 0 >( entity ) );
      }

      /** \brief  @copydoc Dune::Grid::preAdapt() */
      bool preAdapt ()
      {
        return hostGrid().preAdapt();
      }

      /** \brief  @copydoc Dune::Grid::adapt() */
      bool adapt ()
      {
        bool ret = hostGrid().adapt();
        update();
        return ret;
      }

      /** \brief  @copydoc Dune::Grid::adapt()
          \param handle handler for restriction and prolongation operations
          which is a Model of the AdaptDataHandleInterface class.
      */
      template< class GridImp, class DataHandle >
      bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &datahandle )
      {
        typedef InterfaceGridAdaptDataHandle< Grid,
                AdaptDataHandleInterface< GridImp, DataHandle > > WrappedDataHandle;

        WrappedDataHandle wrappedDataHandle( extraData(), datahandle );
        const bool ret = hostGrid().adapt( wrappedDataHandle );
        update();
        return ret;
      }

      /** \brief  @copydoc Dune::Grid::postAdapt() */
      void postAdapt ()
      {
        hostGrid().postAdapt();
      }

      /** \name Parallel Data Distribution and Communication Methods
       *  \{ */

      /** \brief obtain size of overlap region for the leaf grid
       *
       *  \param[in]  codim  codimension for with the information is desired
       */
      int overlapSize ( int codim ) const
      {
        return hostGrid().overlapSize( codim );
      }

      /** \brief obtain size of ghost region for the leaf grid
       *
       *  \param[in]  codim  codimension for with the information is desired
       */
      int ghostSize( int codim ) const
      {
        return hostGrid().ghostSize( codim );
      }

      /** \brief obtain size of overlap region for a grid level
       *
       *  \param[in]  level  grid level (0, ..., maxLevel())
       *  \param[in]  codim  codimension (0, ..., dimension)
       */
      int overlapSize ( int level, int codim ) const
      {
        return hostGrid().overlapSize( level, codim );
      }

      /** \brief obtain size of ghost region for a grid level
       *
       *  \param[in]  level  grid level (0, ..., maxLevel())
       *  \param[in]  codim  codimension (0, ..., dimension)
       */
      int ghostSize ( int level, int codim ) const
      {
        return hostGrid().ghostSize( level, codim );
      }

      /** \brief communicate information on a grid level
       *
       *  \param      dataHandle  communication data handle (user defined)
       *  \param[in]  interface   communication interface (one of
       *                          InteriorBorder_InteriorBorder_Interface,
       *                          InteriorBorder_All_Interface,
       *                          Overlap_OverlapFront_Interface,
       *                          Overlap_All_Interface,
       *                          All_All_Interface)
       *  \param[in]  direction   communication direction (one of
       *                          ForwardCommunication or BackwardCommunication)
       *  \param[in]  level       grid level to communicate
       */
      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction,
                         int level ) const
      {
        levelGridView( level ).communicate( dataHandle, interface, direction );
      }

      /** \brief communicate information on leaf entities
       *
       *  \param      dataHandle  communication data handle (user defined)
       *  \param[in]  interface   communication interface (one of
       *                          InteriorBorder_InteriorBorder_Interface,
       *                          InteriorBorder_All_Interface,
       *                          Overlap_OverlapFront_Interface,
       *                          Overlap_All_Interface,
       *                          All_All_Interface)
       *  \param[in]  direction   communication direction (one of
       *                          ForwardCommunication, BackwardCommunication)
       */
      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                         InterfaceType interface,
                         CommunicationDirection direction ) const
      {
        leafGridView().communicate( dataHandle, interface, direction );
      }

      /** \brief obtain CollectiveCommunication object
       *
       *  The CollectiveCommunication object should be used to globally
       *  communicate information between all processes sharing this grid.
       *
       *  \note The CollectiveCommunication object returned is identical to the
       *        one returned by the host grid.
       */
      const CollectiveCommunication &comm () const
      {
        return hostGrid().comm();
      }

      // data handle interface different between geo and interface

      /** \brief rebalance the load each process has to handle
       *
       *  A parallel grid is redistributed such that each process has about
       *  the same load (e.g., the same number of leaf entites).
       *
       *  \note DUNE does not specify, how the load is measured.
       *
       *  \returns \b true, if the grid has changed.
       */
      bool loadBalance ()
      {
        const bool gridChanged = hostGrid().loadBalance();
        if( gridChanged )
          update();
        return gridChanged;
      }

      /** \brief rebalance the load each process has to handle
       *
       *  A parallel grid is redistributed such that each process has about
       *  the same load (e.g., the same number of leaf entites).
       *
       *  The data handle is used to communicate the data associated with
       *  entities that move from one process to another.
       *
       *  \note DUNE does not specify, how the load is measured.
       *
       *  \param  datahandle  communication data handle (user defined)
       *
       *  \returns \b true, if the grid has changed.
       */

      template< class DataHandle, class Data >
      bool loadBalance ( CommDataHandleIF< DataHandle, Data > &datahandle )
      {
        typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
        typedef InterfaceGridDataHandle< DataHandleIF, Grid > WrappedDataHandle;

        WrappedDataHandle wrappedDataHandle( datahandle );
        typename WrappedDataHandle::DataHandleIF &wrappedDataHandleIF = wrappedDataHandle;
        const bool gridChanged = hostGrid().loadBalance( wrappedDataHandleIF );
        if ( gridChanged )
          update();
        return gridChanged;
      }

      /** \brief rebalance the load each process has to handle
       *
       *  A parallel grid is redistributed such that each process has about
       *  the same load (e.g., the same number of leaf entites).
       *
       *  The data handle is used to communicate the data associated with
       *  entities that move from one process to another.
       *
       *  \note DUNE does not specify, how the load is measured.
       *
       *  \param  dataHandle  data handle following the ALUGrid interface
       *
       *  \returns \b true, if the grid has changed.
       */
      template< class DofManager >
      bool loadBalance ( DofManager &dofManager )
      {
        InterfaceGridWrappedDofManager< DofManager, Grid > wrappedDofManager( dofManager );

        const bool gridChanged = hostGrid().loadBalance( wrappedDofManager );
        if( gridChanged )
          update();
        return gridChanged;
      }

      /** \brief View for a grid level for All_Partition */
      LevelGridView levelGridView ( int level ) const
      {
        typedef typename LevelGridView::GridViewImp ViewImp;
        return LevelGridView( ViewImp( *this, hostGrid().levelGridView( level ) ) );
      }

      /** \brief View for the leaf grid for All_Partition */
      LeafGridView leafGridView () const
      {
        typedef typename LeafGridView::GridViewImp ViewImp;
        return LeafGridView( ViewImp( *this, hostGrid().leafGridView() ) );
      }

      /** \brief obtain Entity from EntitySeed. */
      template< class EntitySeed >
      typename Traits::template Codim< EntitySeed::codimension >::Entity
      entity ( const EntitySeed &seed ) const
      {
        typedef typename Traits::template Codim< EntitySeed::codimension >::EntityImpl EntityImpl;
        return EntityImpl( extraData(), hostGrid().entity( seed.impl().hostEntitySeed() ) );
      }

      /** \} */

      /** \name Miscellaneous Methods
       *  \{ */

      const HostGrid &hostGrid () const { return *hostGrid_; }
      HostGrid &hostGrid () { return *hostGrid_; }

      /** \brief update grid caches
       *
       *  This method has to be called whenever the underlying host grid changes.
       *
       *  \note If you adapt the host grid through this geometry grid's
       *        adaptation or load balancing methods, update is automatically
       *        called.
       */
      void update ()
      {
        const int newNumLevels = maxLevel()+1;
        const int oldNumLevels = levelIndexSets_.size();

        for( int i = newNumLevels; i < oldNumLevels; ++i )
        {
          if( !levelIndexSets_[ i ] )
            delete levelIndexSets_[ i ];
        }
        levelIndexSets_.resize( newNumLevels, nullptr );
      }

      /** \} */

    public:
      using Base::getRealImplementation;

      template< int codim >
      static const typename HostGrid::template Codim< codim >::Entity &
      getHostEntity( const typename Codim< codim >::Entity &entity )
      {
        return getRealImplementation( entity ).hostEntity();
      }

      typedef typename Traits :: ExtraData ExtraData;
      ExtraData extraData () const  { return ExtraData(); }

    protected:
      HostGrid *const hostGrid_;
      mutable std::vector< LevelIndexSet * > levelIndexSets_;
      mutable LeafIndexSet leafIndexSet_;
      mutable GlobalIdSet globalIdSet_;
      mutable LocalIdSet localIdSet_;
    };



    // InterfaceGrid::Codim
    // -------------

    template< class HostGrid >
    template< int codim >
    struct InterfaceGrid< HostGrid >::Codim
    : public Base::template Codim< codim >
    {
      /** \name Entity and Entity Pointer Types
       *  \{ */

      /** \brief type of entity
       *
       *  The entity is a model of Dune::Entity.
       */
      typedef typename Traits::template Codim< codim >::Entity Entity;

      /** \} */

      /** \name Geometry Types
       *  \{ */

      /** \brief type of world geometry
       *
       *  Models the geomtry mapping of the entity, i.e., the mapping from the
       *  reference element into world coordinates.
       *
       *  The geometry is a model of Dune::Geometry, implemented through the
       *  generic geometries provided by dune-grid.
       */
      typedef typename Traits::template Codim< codim >::Geometry Geometry;

      /** \brief type of local geometry
       *
       *  Models the geomtry mapping into the reference element of dimension
       *  \em dimension.
       *
       *  The local geometry is a model of Dune::Geometry, implemented through
       *  the generic geometries provided by dune-grid.
       */
      typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

      /** \} */

      /** \name Iterator Types
       *  \{ */

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename Traits::template Codim< codim >
          ::template Partition< pitype >::LeafIterator
          LeafIterator;
        typedef typename Traits::template Codim< codim >
          ::template Partition< pitype >::LevelIterator
          LevelIterator;
      };

      /** \brief type of level iterator
       *
       *  This iterator enumerates the entites of codimension \em codim of a
       *  grid level.
       *
       *  The level iterator is a model of Dune::LevelIterator.
       */
      typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

      /** \brief type of leaf iterator
       *
       *  This iterator enumerates the entites of codimension \em codim of the
       *  leaf grid.
       *
       *  The leaf iterator is a model of Dune::LeafIterator.
       */
      typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

      /** \} */
    };

  } // namespace VoF

} // namespace Dune

#include <dune/vof/interfacegrid/persistentcontainer.hh>

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GRID_HH
