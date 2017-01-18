#ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
#define DUNE_VOF_INTERFACEGRID_ENTITY_HH

#include <cassert>
#include <cstddef>

#include <limits>
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



    // BasicInterfaceGridEntity
    // ------------------------

    template< int cd, int dim, class Grid >
    class BasicInterfaceGridEntity
    {
      typedef InterfaceGridEntity< cd, dim, Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

      static_assert( dim == Traits::dimension, "Internal Dune Error" );

    protected:
      static const int dimensionworld = Traits::dimensionworld;

    public:
      static const int codimension = cd;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;
      typedef typename Traits::Reconstruction::GridView::template Codim< 0 >::Entity HostElement;

      BasicInterfaceGridEntity () = default;

      BasicInterfaceGridEntity ( const DataSet &dataSet, const HostElement &hostElement )
        : dataSet_( &dataSet ), hostElement_( hostElement )
      {}

      BasicInterfaceGridEntity ( const DataSet &dataSet, HostElement &&hostElement )
        : dataSet_( &dataSet ), hostElement_( std::move( hostElement ) )
      {}

      int level () const { return 0; }

      PartitionType partitionType () const { return hostElement().partitionType(); }

      GeometryType type () const { return GeometryType( (mydimension < 2 ? GeometryType::cube, GeometryType::none), mydimension ); }

      const DataSet &dataSet () const { assert( dataSet_ ); return *dataSet_; }
      const HostElement &hostElement () const { return hostElement_; }

    private:
      const DataSet *dataSet_ = nullptr;
      HostElement hostElement_;
    };



    // InterfaceGridEntity
    // -------------------

    template< int cd, int dim, class Grid >
    class InterfaceGridEntity
      : public BasicInterfaceGridEntity< cd, dim, Grid >
    {
      typedef InterfaceGridEntity< cd, dim, Grid > This;
      typedef BasicInterfaceGridEntity< cd, dim, Grid > Base;

    public:
      using Base::codimension;
      using Base::dimension;

      typedef Dune::EntitySeed< Grid, InterfaceGridEntitySeed< codimension, Grid > > EntitySeed;
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      typedef typename Base::DataSet DataSet;
      typedef typename Base::HostElement HostElement;

      using Base::hostElement;

      InterfaceGridEntity () = default;

      InterfaceGridEntity ( const DataSet &dataSet, const HostElement &hostElement, int subEntity )
        : Base( dataSet, hostElement ), subEntity_( subEntity )
      {}

      InterfaceGridEntity ( const DataSet &dataSet, HostElement &&hostElement, int subEntity )
        : Base( dataSet, std::move( hostElement ) ), subEntity_( subEntity )
      {}

      bool equals ( const This &other ) const { return (hostElement() == other.hostElement()) && (subEntity() == other.subEntity()); }

      Geometry geometry () const
      {
        // TODO: Please implement me
      }

      EntitySeed seed () const { return InterfaceGridEntitySeed< codimension, Grid >( hostElement().seed(), subEntity() ); }

      unsigned int subEntities ( unsigned int codim ) const
      {
        assert( (codim >= static_cast< unsigned int >( codimension )) && (codim <= static_cast< unsigned int >( dimension )) );
        return (codim - static_cast< unsigned int >( codimension ) + 1u);
      }

      int subEntity () const { return subEntity_; }

    private:
      int subEntity_;
    };



    // InterfaceGridEntity for codimension 0
    // -------------------------------------

    template< int dim, class Grid >
    class InterfaceGridEntity< 0, dim, Grid >
      : public BasicInterfaceGridEntity< 0, dim, Grid >
    {
      typedef InterfaceGridEntity< cd, dim, Grid > This;
      typedef BasicInterfaceGridEntity< cd, dim, Grid > Base;

    protected:
      using Base::dimensionworld;

    public:
      using Base::codimension;
      using Base::mydimension;
      using Base::dimension;

      typedef Dune::EntitySeed< Grid, InterfaceGridEntitySeed< codimension, Grid > > EntitySeed;
      typedef Dune::Geometry< mydimension, dimensionworld, Grid, InterfaceGridGeometry > Geometry;
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

      template< int codim >
      struct Codim
      {
        typedef Dune::Entity< codim, dimension, Grid, InterfaceGridEntity > Entity;
      }

      typedef Dune::EntityIterator< 0, Grid, InterfaceGridHierarchicIterator< Grid > > HierarchicIterator;

      typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

      typedef Base::DataSet DataSet;
      typedef Base::HostElement HostElement;

      using Base::dataSet;
      using Base::hostElement;

      InterfaceGridEntity () = default;

      InterfaceGridEntity ( const DataSet &dataSet, const HostElement &hostElement )
        : Base( dataSet, hostElement )
      {}

      InterfaceGridEntity ( const DataSet &dataSet, HostElement &&hostElement )
        : Base( dataSet, hostElement )
      {}

      bool equals ( const This &other ) const { return (hostElement() == other.hostElement()); }

      typename Codim< 0 >::Entity father () const { DUNE_THROW( GridError, "InterfaceGrid consists of only one level" ); }

      Geometry geometry () const
      {
        typedef InterfaceGridGeometry< mydimension, dimensionworld, Grid > Impl;
        const auto elementIndex = dataSet().indices().index( hostElement() );
        const std::size_t index = dataSet().offsets()[ elementIndex ];
        const std::size_t size = dataSet().offsets()[ elementIndex + 1 ] - index;
        return Impl( normal(), dataSet().vertices().data() + index, size );
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

      EntitySeed seed () const { return InterfaceGridEntitySeed< codimension, Grid >( hostElement_.seed() ); }

      template< int codim >
      typename Codim< codim >::Entity subEntity ( int i ) const
      {
        assert( (i >= 0) && (i < static_cast< int >( subEntities( codim ) ) );
        return subEntity( i, Dune::Codim< codim >() );
      }

      unsigned int subEntities ( unsigned int codim ) const
      {
        if( codim > 0u )
        {
          assert( codim <= static_cast< unsigned int >( dimension ) );
          const auto elementIndex = dataSet().indices().index( hostElement() );
          return (dataSet().offsets()[ elementIndex + 1 ] - dataSet().offsets()[ elementIndex ]);
        }
        else
          return 1;
      }

      GlobalCoordinate normal () const { return dataSet().reconstructionSet()[ hostElement() ].innerNormal(); }

    private:
      typename Codim< 0 >::Entity subEntity ( int i, Dune::Codim< 0 > ) const { return *this; }

      template< int codim >
      typename Codim< codim >::Entity subEntity ( int i, Dune::Codim< codim > ) const
      {
        return InterfaceGridEntity< codim, dimension, Grid >( dataSet(), hostElement(), i );
      }
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ENTITY_HH
