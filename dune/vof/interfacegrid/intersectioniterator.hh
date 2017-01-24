#ifndef DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH
#define DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH

#include <type_traits>

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/vof/interfacegrid/entity.hh>
#include <dune/vof/interfacegrid/intersection.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIntersectionIterator
    // ---------------------------------

    template< class Grid >
    class InterfaceGridIntersectionIterator
    {
      typedef InterfaceGridIntersectionIterator< Grid > This;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int dimension = Traits::dimension;

      typedef Dune::Intersection< Grid, InterfaceGridIntersection< Grid > > Intersection;
      typedef Dune::Entity< 0, dimension, Grid, InterfaceGridEntity > Entity;

      InterfaceGridIntersectionIterator () = default;

      InterfaceGridIntersectionIterator ( const Entity &inside, int indexInInside ) : inside_( inside ), indexInInside_( indexInInside ) {}

      bool equals ( const This &other ) const { return (inside_ == other.inside_) && (indexInInside_ == other.indexInInside_); }

      void increment () { ++indexInInside_; }

      Intersection dereference () const { return InterfaceGridIntersection< Grid >( inside_, indexInInside_ ); }

    protected:
      Entity inside_;
      int indexInInside_ = -1;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_INTERSECTIONITERATOR_HH
