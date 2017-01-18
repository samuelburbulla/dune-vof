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
#include <dune/vof/mixedcellmapper.hh>

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

      typedef typename DataSet::ColorFunction ColorFunction;

      typedef MixedCellMapper< typename Reconstruction::GridView > Indices;

      typedef typename Indices::Index IndexType;

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

      typedef typename InterfaceGridIndexSetTraits< Grid >::Indices Indices;

    public:
      static const int dimension = InterfaceGridIndexSetTraits< Grid >::dimension;

      typedef typename InterfaceGridIndexSetTraits< Grid >::DataSet DataSet;

      typedef typename InterfaceGridIndexSetTraits< Grid >::ColorFunction ColorFunction;

      typedef typename InterfaceGridIndexSetTraits< Grid >::IndexType IndexType;

      typedef typename InterfaceGridIndexSetTraits< Grid >::Types Types;

      template< int codim >
      using Codim = InterfaceGridIndexSetTraits< Grid >::template Codim< codim >;

      explicit InterfaceGridIndexSet ( const ColorFunction &colorFunction, Args &&... args )
        : dataSet_( colorFunction, std::forward< Args >( args )... ), indices_( dataSet_.flags() )
      {}

      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      template< int cd >
      IndexType index ( const typename Codim< cd >::Entity &entity ) const
      {
        // TODO: extend to (dimension == 2), i.e., the 3d case
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
        // TODO: extend to (dimension == 2), i.e., the 3d case
        return subIndex( entity, i, codim, Dune::Codim< cd >() );
      }

      IndexType size ( GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        return (type == types( codim ) ? size( codim ) : static_cast< IndexType >( 0u ));
      }

      IndexType size ( int codim ) const
      {
        // TODO: extend to (dimension == 2), i.e., the 3d case
        assert( dimension == 1 );

        assert( (codim >= 0) && (codim <= dimension) )
        return (codim+1) * size_;
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
        return {{ GeometryType( (mydim < 2 ? GeometryType::cube, GeometryType::none), mydim ) }};
      }

      void update ( const ColorFunction &colorFunction )
      {
        dataSet_.update( colorFunction );
        indices_.update( dataSet_.flags() );
      }

      const DataSet &dataSet () const { return dataSet_; }

    private:
      IndexType index ( const typename Codim< 0 >::Entity &entity, Dune::Codim< 0 > ) const
      {
        return indices_.index( Grid::getRealImplementation( entity ).hostElement() );
      }

      IndexType index ( const typename Codim< 1 >::Entity &entity, Dune::Codim< 1 > ) const
      {
        return 2*indices_.index( Grid::getRealImplementation( entity ).hostElement() ) + Grid::getRealImplementation( entity ).subEntity();
      }

      IndexType subIndex ( const typename Codim< 0 >::Entity &entity, int i, unsigned int codim, Dune::Codim< 0 > ) const
      {
        assert( codim <= static_cast< unsigned int >( dimension ) );
        assert( (i >= 0) && (i < static_cast< int >( codim+1 )) );
        return (codim+1)*indices_.index( Grid::getRealImplementation( entity ).hostElement() ) + i;
      }

      IndexType subIndex ( const typename Codim< 1 >::Entity &entity, int i, unsigned int codim, Dune::Codim< 1 > ) const
      {
        assert( (codim == 1u) && (i == 0) );
        return index( entity, Dune::Codim< 1 >() );
      }

      DataSet dataSet_;
      Indices indices_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INDEXSET_HH
