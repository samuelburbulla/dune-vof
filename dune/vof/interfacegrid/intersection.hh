#ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
#define DUNE_VOF_INTERFACEGRID_INTERSECTION_HH

#include <cstddef>

#include <type_traits>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/vof/interfacegrid/dataset.hh>
#include <dune/vof/interfacegrid/entity.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIntersection
    // -------------------------

    template< class Grid >
    class InterfaceGridIntersection
    {
      typedef InterfaceGridIntersection< Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension-1;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::ctype ctype;

      typedef Dune::Entity< 0, dimension, Grid, InterfaceGridEntity > Entity;
      typedef Dune::Geometry< mydimension, dimensionworld, Grid, InterfaceGridGeometry > Geometry;
      typedef Dune::Geometry< mydimension, dimension, Grid, InterfaceGridGeometry > LocalGeometry;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, dimensionworld > GlobalCoordinate;

      typedef InterfaceGridDataSet< typename Traits::Reconstruction > DataSet;

      InterfaceGridIntersection () = default;

      InterfaceGridIntersection ( const Entity &inside, int indexInInside ) : inside_( inside ), indexInInside_( indexInInside ) {}

      bool equals ( const This &other ) const { return (inside_ == other.inside_) && (indexInInside_ == other.indexInInside_); }

      Entity inside () const { return inside_; }
      Entity outside () const { return Entity(); }

      bool boundary () const { return true; }
      bool conforming () const { return true; }
      bool neighbor () const { return false; }

      int boundaryId () const { return 1; }

      std::size_t boundarySegmentIndex () const
      {
        const auto elementIndex = dataSet().indices().index( Grid::getRealImplementation( inside_ ).hostElement() );
        return dataSet().offsets()[ elementIndex ] + static_cast< std::size_t >( indexInInside() );
      }

      LocalGeometry geometryInInside () const
      {
        return geometryInInside( Dune::Dim< mydimension >() );
      }

      LocalGeometry geometryInOutside () const
      {
        DUNE_THROW( InvalidStateException, "InterfaceGrid does not have neighbors." );
      }

      Geometry geometry () const
      {
        return Geometry( dataSet().geometry( Grid::getRealImplementation( inside_ ).hostElement(), indexInInside(), Dune::Dim< mydimension >() ) );
      }

      GeometryType type () const { return mydimension < 2 ? GeometryTypes::cube( mydimension ) : GeometryTypes::none( mydimension ); }

      int indexInInside () const { return indexInInside_; }
      int indexInOutside () const { return -1; }

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate &local ) const
      {
        return (mydimension > 0 ? outerNormal( local ) : unitOuterNormal( local ));
      }

      GlobalCoordinate outerNormal ( const LocalCoordinate & ) const
      {
        GlobalCoordinate normal;
        dataSet().covariantOuterNormal( Grid::getRealImplementation( inside_ ).hostElement(), indexInInside_, normal );
        return normal;
      }

      GlobalCoordinate unitOuterNormal ( const LocalCoordinate & ) const { return centerUnitOuterNormal(); }

      GlobalCoordinate centerUnitOuterNormal () const
      {
        GlobalCoordinate normal;
        dataSet().covariantOuterNormal( Grid::getRealImplementation( inside_ ).hostElement(), indexInInside_, normal );
        return normal /= normal.two_norm();
      }

      const DataSet &dataSet () const { return Grid::getRealImplementation( inside_ ).dataSet(); }

    protected:
      LocalGeometry geometryInInside ( Dune::Dim< 0 > ) const
      {
        return LocalGeometry( InterfaceGridGeometry< 0, dimension, Grid >( FieldVector< ctype, dimension >{ ctype( indexInInside() ) } ) );
      }

      LocalGeometry geometryInInside ( Dune::Dim< 1 > ) const
      {
        DUNE_THROW( InvalidStateException, "Intersection::geometryInInside does not make for arbitrary polytopes." );
      }

      Entity inside_;
      int indexInInside_ = -1;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
