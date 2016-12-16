#ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
#define DUNE_VOF_INTERFACEGRID_ENTITY_HH

#include <type_traits>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/entity.hh>

namespace Dune
{

  namespace VoF
  {

    // External Forward Declarations
    // -----------------------------

    template< class, class > class InterfaceGridIntersection;
    template< class, class > class InterfaceGridIntersectionIterator;

    template< class, class > class InterfaceGridIterator;


    // InterfaceGridEntityBasic
    // -----------------

    /** \copydoc InterfaceGridEntity
     *
     *  \nosubgrouping
     */
    template< int codim, int dim, class Grid >
    class InterfaceGridEntityBasic
    {
    protected:
      typedef typename std::remove_const< Grid >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = codim;
      //! dimension of the grid
      static const int dimension = Traits::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = Traits::dimensionworld;

      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename Traits::ctype ctype;

      //! type of corresponding entity seed
      typedef typename Grid::template Codim< codimension >::EntitySeed EntitySeed;
      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      /** \} */

    protected:
      // type of the host grid
      typedef typename Traits::HostGrid  HostGrid;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

      // type of entity implementation
      typedef typename Traits::template Codim< codimension > :: EntityImpl  EntityImpl;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity */
      explicit InterfaceGridEntityBasic ( ExtraData data )
      : hostEntity_(),
        data_ ( data )
      {}

      /** \brief construct an initialized entity
       *
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note The reference to the host entity must remain valid  as long as
       *        this entity is in use.
       */
      InterfaceGridEntityBasic ( ExtraData data, const HostEntity &hostEntity )
      : hostEntity_( hostEntity ),
        data_( data )
      {}

      /** \} */

      /** \brief return true if entity hold a vaild host entity */
      // operator bool () const { return bool( hostEntity_ ); }

      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        return hostEntity().type();
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostEntity().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      /** obtain the geometry of this entity */
      Geometry geometry () const
      {
        return Geometry( hostEntity().geometry() );
      }

      /** \brief return EntitySeed of host grid entity */
      EntitySeed seed () const { return typename EntitySeed::Implementation( hostEntity().seed() ); }

      bool equals( const EntityImpl& other ) const { return hostEntity() ==  other.hostEntity(); }

      /** \} */


      /** \name Methods Supporting the Grid Implementation
       *  \{ */

      const HostEntity &hostEntity () const
      {
        return hostEntity_;
      }

      ExtraData data() const { return data_; }

      /** \} */

    protected:
      HostEntity hostEntity_;
      ExtraData        data_;
    };



    // InterfaceGridEntity
    // ------------

    template< int codim, int dim, class Grid >
    class InterfaceGridEntity : public InterfaceGridEntityBasic< codim, dim, Grid >
    {
      typedef InterfaceGridEntityBasic< codim, dim, Grid > Base ;
    protected:
      typedef typename std::remove_const< Grid >::type::Traits Traits;

    protected:
      // type of the host grid
      typedef typename Traits::HostGrid  HostGrid;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

    public:
      using Base :: codimension ;

      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
      /** \} */

      explicit InterfaceGridEntity ( ExtraData data )
      : Base( data )
      {}

      InterfaceGridEntity ( ExtraData data, const HostEntity &hostEntity )
      : Base( data, hostEntity )
      {}
    };


    // InterfaceGridEntity for codimension 0
    // ----------------------------------

    /** \copydoc InterfaceGridEntity
     *
     *  \nosubgrouping
     */
    template< int dim, class Grid >
    class InterfaceGridEntity< 0, dim, Grid > : public InterfaceGridEntityBasic< 0, dim, Grid >
    {
      typedef InterfaceGridEntityBasic< 0, dim, Grid > Base ;
    protected:
      typedef typename Base::Traits Traits;
      typedef typename Base::HostGrid HostGrid;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Base::ExtraData ExtraData;

    public:
      using Base::codimension ;
      using Base::data ;
      using Base::hostEntity ;
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
      /** \} */

    protected:
      typedef typename Traits :: LeafIntersectionIteratorImpl  LeafIntersectionIteratorImpl;
      typedef typename Traits :: LevelIntersectionIteratorImpl LevelIntersectionIteratorImpl;

      typedef typename Traits :: LeafIntersectionImpl          LeafIntersectionImpl;
      typedef typename Traits :: LevelIntersectionImpl         LevelIntersectionImpl;

    public:
      /** \name Types Required by DUNE
       *  \{ */

      //! type of corresponding local geometry
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

      //! type of Entity interface
      typedef typename Traits ::template Codim< codimension >::Entity Entity;

      //! type of hierarchic iterator
      typedef typename Traits::HierarchicIterator        HierarchicIterator;
      //! type of leaf intersection iterator
      typedef typename Traits::LeafIntersectionIterator  LeafIntersectionIterator;
      //! type of level intersection iterator
      typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity */
      explicit InterfaceGridEntity ( ExtraData data )
      : Base( data )
      {}

      /** \brief construct an initialized entity
       *
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note The reference to the host entity must remain valid as long as
       *        this entity is in use.
       */
      InterfaceGridEntity ( ExtraData data, const HostEntity &hostEntity )
      : Base( data, hostEntity )
      {}

      /** \} */

      unsigned int subEntities( const unsigned int codim ) const
      {
        return hostEntity().subEntities( codim );
      }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      template< int codim >
      typename Grid::template Codim< codim >::Entity
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::EntityImpl EntityImpl;
        return EntityImpl( data(), hostEntity().template subEntity< codim >( i ) );
      }

      LevelIntersectionIterator ilevelbegin () const
      {
        return LevelIntersectionIteratorImpl( data(), hostEntity().ilevelbegin() );
      }

      LevelIntersectionIterator ilevelend () const
      {
        return LevelIntersectionIteratorImpl( data(), hostEntity().ilevelend() );
      }

      LeafIntersectionIterator ileafbegin () const
      {
        return LeafIntersectionIteratorImpl( data(), hostEntity().ileafbegin() );
      }

      LeafIntersectionIterator ileafend () const
      {
        return LeafIntersectionIteratorImpl( data(), hostEntity().ileafend() );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      bool isLeaf () const
      {
        return hostEntity().isLeaf();
      }

      Entity father () const { DUNE_THROW( GridError, "InterfaceGrid consists of only one level" ); }

      bool hasFather () const { return false; }

      LocalGeometry geometryInFather () const
      {
        return LocalGeometry( hostEntity().geometryInFather() );
      }

      HierarchicIterator hbegin ( int maxLevel ) const
      {
        typedef typename Traits :: HierarchicIteratorImpl HierarchicIteratorImpl ;
        return HierarchicIteratorImpl( data(), hostEntity().hbegin( maxLevel ) );
      }

      HierarchicIterator hend ( int maxLevel ) const
      {
        typedef typename Traits :: HierarchicIteratorImpl HierarchicIteratorImpl ;
        return HierarchicIteratorImpl( data(), hostEntity().hend( maxLevel ) );
      }

      bool isRegular () const { return true; }
      bool isNew () const { return false; }
      bool mightVanish () const { return false; }
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
