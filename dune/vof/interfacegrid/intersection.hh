#ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
#define DUNE_VOF_INTERFACEGRID_INTERSECTION_HH

#include <cstddef>

#include <type_traits>

#include <dune/grid/common/intersection.hh>

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
      typedef typename Traits::template Codim< 1 >::Geometry       Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry  LocalGeometry;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, dimensionworld > GlobalCoordinate;

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
        const IndexType elementIndex = dataSet().indices().index( Grid::getRealImplementation( entity ).hostElement() );
        return dataSet().offsets()[ elementIndex ] + static_cast< std::size_t >( indexInInside() );
      }

      LocalGeometry geometryInInside () const
      {
        // TODO: Please implement me
      }

      LocalGeometry geometryInOutside () const
      {
        // TODO: Please implement me
      }

      Geometry geometry () const
      {
        // TODO: Please implement me
      }

      GeometryType type () const
      {
        GeometryType( (mydimension < 2 ? GeometryType::cube, GeometryType::none), mydimension );
      }

      int indexInInside () const { return indexInInside_; }
      int indexInOutside () const { return -1; }

      GlobalCoordinate integrationOuterNormal ( const LocalCoordinate &local ) const
      {
        // TODO: Please implement me
      }

      GlobalCoordinate outerNormal ( const LocalCoordinate & ) const { return centerUnitOuterNormal(); }
      GlobalCoordinate unitOuterNormal ( const LocalCoordinate & ) const { return centerUnitOuterNormal(); }

      GlobalCoordinate centerUnitOuterNormal () const
      {
        // TODO: Please implement me
      }

    protected:
      Entity inside_;
      int indexInInside_ = -1;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INTERSECTION_HH
