#ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
#define DUNE_VOF_GEOMETRY_HALFSPACE_HH

namespace Dune
{
  namespace VoF
  {

    // HalfSpace
    // ---------

    template < class Coord >
    struct HalfSpace
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;

      HalfSpace () = default;
      HalfSpace ( const HalfSpace &other ) = default;

      HalfSpace ( const Coordinate& normal, const ctype dist )
       : normal_( normal ), dist_( dist )
      {}

      HalfSpace ( const Coordinate& normal, const Coordinate& point )
       : normal_( normal ), dist_( -1.0*( normal * point ) )
      {}

      const Coordinate& normal () const { return normal_; }
      const ctype& distance () const { return dist_; }

      Coordinate& normal () DUNE_DEPRECATED_MSG( "Mutability will be removed." ) { return normal_; }
      ctype& distance () DUNE_DEPRECATED_MSG( "Mutability will be removed." ) { return dist_; }

    private:
      Coordinate normal_;
      ctype dist_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_HALFSPACE_HH
