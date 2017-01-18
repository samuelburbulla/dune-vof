#ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
#define DUNE_VOF_INTERFACEGRID_ENTITY_HH

#include <cassert>

#include <type_traits>
#include <utility>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/entity.hh>

#include <dune/vof/interfacegrid/dataset.hh>
#include <dune/vof/interfacegrid/entityseed.hh>

namespace Dune
{

  namespace VoF
  {

    // External Forward Declarations
    // -----------------------------

    template< class Grid >
    class InterfaceGridHierarchicIterator;



    // InterfaceGridEntity
    // -------------------

    template< int cd, int dim, class Grid >
    class InterfaceGridEntity< cd, dim, Grid >
    {
      typedef InterfaceGridEntity< cd, dim, Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

      static_assert( dim == Traits::dimension, "Internal Dune Error" );

    public:
      static const int codimension = cd;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;

      typedef Dune::EntitySeed< Grid, InterfaceGridEntitySeed< codimension, Grid > > EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;
      typedef typename Traits::Reconstruction::GridView::template Codim< 0 >::Entity HostElement;

      InterfaceGridEntity () = default;

      InterfaceGridEntity ( const DataSet &dataSet, const HostElement &hostElement, int subEntity )
        : dataSet_( &dataSet ), hostElement_( hostElement ), subEntity_( subEntity )
      {}

      InterfaceGridEntity ( const DataSet &dataSet, HostElement &&hostElement, int subEntity )
        : dataSet_( &dataSet ), hostElement_( std::move( hostElement ) ), subEntity_( subEntity )
      {}

      bool equals ( const This &other ) const { return (hostElement() == other.hostElement()) && (subEntity() == other.subEntity()); }

      Geometry geometry () const
      {
        // TODO: Please implement me
      }

      int level () const { return 0; }

      PartitionType partitionType () const { return hostElement().partitionType(); }

      EntitySeed seed () const { return InterfaceGridEntitySeed< codimension, Grid >( hostElement_.seed(), subEntity_ ); }

      unsigned int subEntities ( unsigned int codim )
      {
        assert( (codim >= static_cast< unsigned int >( codimension )) && (codim <= static_cast< unsigned int >( dimension )) );
        return 1;
      }

      GeometryType type () const { return GeometryType( (mydimension < 2 ? GeometryType::cube, GeometryType::none), mydimension ); }

      const DataSet &dataSet () const { assert( dataSet_ ); return *dataSet_; }
      const HostElement &hostElement () const { return hostElement_; }

      int subEntity () const { return subEntity_; }

    private:
      const DataSet *dataSet_ = nullptr;
      HostElement hostElement_;
      int subEntity_;
    };



    // InterfaceGridEntity for codimension 0
    // -------------------------------------

    template< int dim, class Grid >
    class InterfaceGridEntity< 0, dim, Grid >
    {
      typedef InterfaceGridEntity< cd, dim, Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

      static_assert( dim == Traits::dimension, "Internal Dune Error" );

    public:
      static const int codimension = 0;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;

      typedef Dune::EntitySeed< Grid, InterfaceGridEntitySeed< codimension, Grid > > EntitySeed;
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      template< int codim >
      struct Codim
      {
        typedef Dune::Entity< codim, dimension, Grid, InterfaceGridEntity > Entity;
      }

      typedef Dune::EntityIterator< 0, Grid, InterfaceGridHierarchicIterator< Grid > > HierarchicIterator;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;
      typedef typename Traits::Reconstruction::GridView::template Codim< 0 >::Entity HostElement;

      InterfaceGridEntity () = default;

      explicit InterfaceGridEntity ( const HostElement &hostElement )
        : hostElement_( hostElement )
      {}

      explicit InterfaceGridEntity ( HostElement &&hostElement )
        : hostElement_( std::move( hostElement ) )
      {}

      bool equals ( const This &other ) const { return (hostElement() == other.hostElement()) && (subEntity() == other.subEntity()); }

      typename Codim< 0 >::Entity father () const { DUNE_THROW( GridError, "InterfaceGrid consists of only one level" ); }

      Geometry geometry () const
      {
        // TODO: Please implement me
      }

      LocalGeometry geometryInFather () const { DUNE_THROW( GridError, "InterfaceGrid consists of only one level" ); }

      bool hasBoundaryIntersections () const { return true; }
      bool hasFather () const { return false; }

      HierarchicIterator hbegin ( int maxLevel ) const { return InterfaceGridHierarchicIterator< Grid >(); }
      HierarchicIterator hend ( int maxLevel ) const { return InterfaceGridHierarchicIterator< Grid >(); }

      bool isLeaf () const { return true; }
      bool isNew () const { return false; }
      bool isRegular () const { return true; }
      bool mightVanish () const { return false; }

      int level () const { return 0; }

      PartitionType partitionType () const { return hostElement().partitionType(); }

      EntitySeed seed () const { return InterfaceGridEntitySeed< codimension, Grid >( hostElement_.seed() ); }

      template< int codim >
      typename Codim< codim >::Entity subEntity ( int i ) const
      {
        assert( (i >= 0) && (i < static_cast< int >( subEntities( codim ) ) );
        return subEntity( i, Dune::Codim< codim >() );
      }

      unsigned int subEntities ( unsigned int codim )
      {
        // TODO: extend to (dimension == 2), i.e., the 3d case
        assert( (codim >= static_cast< unsigned int >( codimension )) && (codim <= static_cast< unsigned int >( dimension )) );
        return (codim+1);
      }

      GeometryType type () const { return GeometryType( (mydimension < 2 ? GeometryType::cube, GeometryType::none), mydimension ); }

      const DataSet &dataSet () const { assert( dataSet_ ); return *dataSet_; }
      const HostElement &hostElement () const { return hostElement_; }

    private:
      typename Codim< 0 >::Entity subEntity ( int i, Dune::Codim< 0 > ) const { return *this; }

      template< int codim >
      typename Codim< codim >::Entity subEntity ( int i, Dune::Codim< codim > ) const
      {
        return InterfaceGridEntity< codim, dimension, Grid >( dataSet(), hostElement(), i );
      }

      const DataSet *dataSet_ = nullptr;
      HostElement hostElement_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
