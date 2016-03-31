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

    /**
     * \ingroup Geometry
     * \brief half space represented as a plane in hesse normal form
     *
     * \tparam  Coord  global coordinate type
     */
    template < class Coord >
    struct HalfSpace
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;
      using Boundary = HyperPlane< Coordinate >;

      static constexpr int dimension = Coordinate::dimension;
      static constexpr int dimensionworld = Coordinate::dimension;

      HalfSpace ()
      : innerNormal_( 0 ), distanceToOrigin_( 0.0 )
      {}

      HalfSpace ( const Coordinate& normal, const ctype& dist )
       : innerNormal_( normal ), distanceToOrigin_( dist )
      {}

      HalfSpace ( const Boundary& plane )
      : HalfSpace( plane.normal_, plane.distanceToOrigin_ )
      {}

      HalfSpace ( const Coordinate& normal, const Coordinate& point )
       : innerNormal_( normal ), distanceToOrigin_( -1.0*( normal * point ) )
      {}

      /**
       * \brief inner normal (used in the representation)
       */
      const Coordinate& innerNormal () const { return innerNormal_; }

      /**
       * \brief bounding hyperplane
       */
      Boundary boundary () const { return Boundary( innerNormal(), distanceToOrigin_ ); }

      /**
       * \brief level set function of the bounding plane
       */
      ctype levelSet ( const Coordinate& point ) const { return innerNormal() * point + distanceToOrigin_; }

      /**
       * \brief outer normal (returns a copy)
       */
      Coordinate outerNormal () const
      {
        Coordinate outer( innerNormal() );
        outer *= -1.0;
        return outer;
      }

      explicit operator bool () const
      {
        return innerNormal_.two_norm2() > 0.5;
      }

    private:
      Coordinate innerNormal_;
      ctype distanceToOrigin_;
    };


    // HyperPlane
    // ----------

    /**
     * \ingroup Geometry
     * \brief plane in hesse normal form
     *
     * \tparam  Coord  global coordinate type
     */
    template< class Coord >
    struct HyperPlane
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;

      static constexpr int dimension = Coordinate::dimension - 1;
      static constexpr int dimensionworld = Coordinate::dimension;

      HyperPlane ()
      : normal_( 0 ), distanceToOrigin_( 0.0 )
      {}

      HyperPlane ( const Coordinate& normal, const ctype& dist )
       : normal_( normal ), distanceToOrigin_( dist )
      {}

      /**
       * \brief level set
       */
      ctype levelSet ( const Coordinate& point ) const { return normal_ * point + distanceToOrigin_; }

      /**
       * \brief cast into half space bounded by this plane
       */
      operator HalfSpace< Coordinate > () const { return HalfSpace< Coordinate >( *this ); }

      explicit operator bool () const
      {
        return normal_.two_norm2() > 0.5;
      }

    private:
      friend HalfSpace< Coordinate >;

      Coordinate normal_;
      ctype distanceToOrigin_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
