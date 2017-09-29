#ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
#define DUNE_VOF_INTERFACEGRID_INDEXSET_HH

#include <cstddef>

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/vof/interfacegrid/dataset.hh>
#include <dune/vof/interfacegrid/declaration.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIndexSetTraits
    // ---------------------------

    template< class Grid >
    struct InterfaceGridIndexSetTraits
    {
      typedef typename std::remove_const_t< Grid >::Traits GridTraits;

      static const int dimension = GridTraits::dimension;

      typedef typename GridTraits::Reconstruction Reconstruction;

      typedef InterfaceGridDataSet< Reconstruction > DataSet;

      typedef typename DataSet::Indices::Index IndexType;

      typedef std::array< GeometryType, 1 > Types;

      template< int codim >
      struct Codim
      {
        typedef typename GridTraits::template Codim< codim >::Entity Entity;
      };
    };



    // InterfaceGridIndexSet
    // ---------------------

    template< class Grid >
    class InterfaceGridIndexSet
      : public IndexSet< Grid, InterfaceGridIndexSet< Grid >, typename InterfaceGridIndexSetTraits< Grid >::IndexType >
    {
      typedef InterfaceGridIndexSet< Grid > This;
      typedef IndexSet< Grid, InterfaceGridIndexSet< Grid >, typename InterfaceGridIndexSetTraits< Grid >::IndexType > Base;

    public:
      static const int dimension = InterfaceGridIndexSetTraits< Grid >::dimension;

      typedef typename InterfaceGridIndexSetTraits< Grid >::DataSet DataSet;

      typedef typename InterfaceGridIndexSetTraits< Grid >::IndexType IndexType;

      typedef typename InterfaceGridIndexSetTraits< Grid >::Types Types;

      template< int codim >
      using Codim = typename InterfaceGridIndexSetTraits< Grid >::template Codim< codim >;

      template< class ColorFunction, class... Args >
      explicit InterfaceGridIndexSet ( const ColorFunction &colorFunction, Args &&... args )
        : dataSet_( colorFunction, std::forward< Args >( args )... )
      {}

      InterfaceGridIndexSet ( This &&other )
        : dataSet_( std::move( other.dataSet_ ) )
      {}

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      template< int cd >
      IndexType index ( const typename Codim< cd >::Entity &entity ) const
      {
        return index ( entity, Dune::Codim< cd >() );
      }

      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const
      {
        return subIndex< Entity::codimension >( entity, i, codim );
      }

      template< int cd >
      IndexType subIndex ( const typename Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        assert( (static_cast< unsigned int >( cd ) <= codim) && (codim <= static_cast< unsigned int >( dimension )) );
        assert( (i >= 0) && (static_cast< unsigned int >( i ) < entity.subEntities( codim )) );
        return subIndex( entity, i, codim, Dune::Codim< cd >() );
      }

      IndexType size ( GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        return (type == types( codim )[ 0 ] ? size( codim ) : static_cast< IndexType >( 0u ));
      }

      IndexType size ( int codim ) const
      {
        assert( (codim >= 0) && (codim <= dimension) );
        const IndexType size = dataSet().indices().size();
        return (codim == 0 ? size : dataSet().offsets()[ size ]);
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return true;
      }

      Types types ( int codim ) const
      {
        assert( (codim >= 0) && (codim <= dimension) );
        const int mydim = dimension - codim;
        return {{ mydim < 2 ? GeometryTypes::cube( mydim ) : GeometryTypes::none( mydim ) }};
      }

      template< class ColorFunction >
      void update ( const ColorFunction &colorFunction )
      {
        dataSet_.update( colorFunction );
      }

      const DataSet &dataSet () const { return dataSet_; }

    private:
      IndexType index ( const typename Codim< 0 >::Entity &entity, Dune::Codim< 0 > ) const
      {
        return dataSet().indices().index( Grid::getRealImplementation( entity ).hostElement() );
      }

      template< int cd >
      IndexType index ( const typename Codim< cd >::Entity &entity, Dune::Codim< cd > ) const
      {
        const IndexType elementIndex = dataSet().indices().index( Grid::getRealImplementation( entity ).hostElement() );
        return dataSet().offsets()[ elementIndex ] + static_cast< IndexType >( Grid::getRealImplementation( entity ).subEntity() );
      }

      IndexType subIndex ( const typename Codim< 0 >::Entity &entity, int i, unsigned int codim, Dune::Codim< 0 > ) const
      {
        const IndexType elementIndex = dataSet().indices().index( Grid::getRealImplementation( entity ).hostElement() );
        return (codim == 0 ? elementIndex : dataSet().offsets()[ elementIndex ] + static_cast< IndexType >( i ));
      }

      template< int cd >
      IndexType subIndex ( const typename Codim< cd >::Entity &entity, int i, unsigned int codim, Dune::Codim< cd > ) const
      {
        const IndexType elementIndex = dataSet().indices().index( Grid::getRealImplementation( entity ).hostElement() );
        const IndexType offset = dataSet().offsets()[ elementIndex ];
        const IndexType size = dataSet().offsets()[ elementIndex + 1 ];
        return offset + (static_cast< IndexType >( i + Grid::getRealImplementation( entity ).subEntity() ) % size);
      }

      DataSet dataSet_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
