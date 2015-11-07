#ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
#define DUNE_VOF_RECONSTRUCTIONSET_HH

//- dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>

//- dune-vof includes
#include <dune/vof/hyperplane.hh>

namespace Dune
{
  namespace VoF
  {

    // ReconstructionSet
    // -----------------

    template< class GV >
    struct ReconstructionSet
    {
      using GridView = GV;
      using Entity = typename GridView::template Codim< 0 >::Entity;

    private:
      using Mapper = MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout >;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      using Reconstruction = Hyperplane< Coordinate >;
      using Intersections = std::vector< Coordinate >;

      using iterator = typename std::vector< Reconstruction >::iterator;
      using const_iterator = typename std::vector< Reconstruction >::const_iterator;

      explicit ReconstructionSet ( const GridView &gridView )
       : mapper_( gridView ), reconstructionSet_( mapper().size() ), intersectionsSet_( mapper().size() )
       {}

      const Reconstruction& operator[] ( const Entity &entity ) const { return reconstructionSet_[ mapper().index( entity ) ]; }
      Reconstruction& operator[] ( const Entity &entity ) { return reconstructionSet_[ mapper().index( entity ) ]; }

      iterator begin () { return reconstructionSet_.begin(); }
      const_iterator begin () const { return reconstructionSet_.begin(); }

      iterator end () { return reconstructionSet_.end(); }
      const_iterator end () const { return reconstructionSet_.end(); }

      const Intersections& intersections ( const Entity &entity ) const { return intersectionsSet_[ mapper().index( entity ) ]; }
      Intersections& intersections ( const Entity &entity ) { return intersectionsSet_[ mapper().index( entity ) ]; }

      const std::vector< Intersections >& intersectionsSet () const { return intersectionsSet_; }

      const void clear()
      {
        std::fill( reconstructionSet_.begin(), reconstructionSet_.end(), Reconstruction() );
        std::fill( intersectionsSet_.begin(), intersectionsSet_.end(), Intersections() );
      }

    private:
      const Mapper &mapper () const { return mapper_; }

      Mapper mapper_;
      std::vector< Reconstruction > reconstructionSet_;
      std::vector< Intersections > intersectionsSet_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
