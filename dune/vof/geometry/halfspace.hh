#ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
#define DUNE_VOF_GEOMETRY_HALFSPACE_HH

#include <cassert>

#include <dune/common/deprecated.hh>

namespace Dune
{
  namespace VoF
  {

    // forward declarations
    // --------------------

    template< class > struct HalfSpace;
    template< class > struct HyperPlane;


    // HalfSpace
    // ---------

    template < class Coord >
    struct HalfSpace
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;
      using Boundary = HyperPlane< Coordinate >;

      static constexpr int dimension = Coordinate::dimension;
      static constexpr int dimensionworld = Coordinate::dimension;

      HalfSpace () = default;
      HalfSpace ( const HalfSpace &other ) = default;

      HalfSpace ( const Coordinate& normal, const ctype& dist )
       : innerNormal_( normal ), distanceToOrigin_( dist )
      {}

      HalfSpace ( const Boundary& plane )
      : HalfSpace( plane.normal_, plane.distanceToOrigin_ )
      {}

      HalfSpace ( const Coordinate& normal, const Coordinate& point ) DUNE_DEPRECATED_MSG( "HalfSpace( normal, point ) might be inconsistent." )
       : innerNormal_( normal ), distanceToOrigin_( -1.0*( normal * point ) )
      {}

      const Coordinate& normal () const DUNE_DEPRECATED_MSG( "Use innerNormal() instead." ) { return innerNormal_; }
      const ctype& distance () const DUNE_DEPRECATED_MSG( "distance() will be removed." ) { return distanceToOrigin_; }

      const Coordinate& innerNormal () const { return innerNormal_; }
      Coordinate outerNormal () const
      {
        Coordinate outer( innerNormal() );
        outer *= -1.0;
        return outer;
      }

      Boundary boundary () const { return Boundary( innerNormal(), distanceToOrigin_ ); }

      ctype levelSet ( const Coordinate& point ) const { return innerNormal() * point + distanceToOrigin_; }

      Coordinate& normal () DUNE_DEPRECATED_MSG( "Mutability will be removed." ) { return innerNormal_; }
      ctype& distance () DUNE_DEPRECATED_MSG( "Mutability will be removed." ) { return distanceToOrigin_; }

    private:
      Coordinate innerNormal_;
      ctype distanceToOrigin_;
    };


    // HyperPlane
    // ----------

    template< class Coord >
    struct HyperPlane
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;

      static constexpr int dimension = Coordinate::dimension - 1;
      static constexpr int dimensionworld = Coordinate::dimension;

      HyperPlane () = default;
      HyperPlane ( const HyperPlane &other ) = default;

      HyperPlane ( const Coordinate& normal, const ctype& dist )
       : normal_( normal ), distanceToOrigin_( dist )
      {}

      ctype levelSet ( const Coordinate& point ) const { return normal_ * point + distanceToOrigin_; }

      operator HalfSpace< Coordinate > () const { return HalfSpace< Coordinate >( *this ); }

    private:
      friend HalfSpace< Coordinate >;

      Coordinate normal_;
      ctype distanceToOrigin_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
