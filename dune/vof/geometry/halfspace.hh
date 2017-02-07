  #ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
#define DUNE_VOF_GEOMETRY_HALFSPACE_HH

#include <cassert>

#include <dune/common/deprecated.hh>

#include <dune/vof/geometry/utility.hh>

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

      HalfSpace ( const std::vector< Coordinate > &points )
      {
        assert( points.size() == dimension );

        std::array< Coordinate, dimension-1 > vectors;
        for ( std::size_t i = 1; i < dimension; ++i )
          vectors[ i-1 ] = points[ i ] - points[ 0 ];

        innerNormal_ = generalizedCrossProduct( vectors );
        innerNormal_ /= innerNormal_.two_norm();

        distanceToOrigin_ = -1.0*( innerNormal_ * points[ 0 ] );
      }

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

      /**
       * \brief flip halfspace (returns a copy)
       */
      HalfSpace< Coordinate > flip () const
      {
        return HalfSpace< Coordinate > ( outerNormal(), - distanceToOrigin_ );
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
       * \brief normal
       */
      const Coordinate& normal () const { return normal_; }

      /**
       * \brief distance to origin
       */
      ctype distance () const { return distanceToOrigin_; }

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
