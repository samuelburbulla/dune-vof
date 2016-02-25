#ifndef DUNE_VOF_GEOMETRY_HYPERPLANE_HH
#define DUNE_VOF_GEOMETRY_HYPERPLANE_HH

namespace Dune
{
  namespace VoF
  {

    // TODO:
    // - move some stuff from geometricutility over here
    // - normalize normals nad distance?

    // Hyperplane
    // ----------

    template < class Coord >
    struct Hyperplane
    {
      using Coordinate = Coord;
      using ctype = typename Coordinate::value_type;  // should be real_type

      Hyperplane ()
       : normal_( 0.0 ), dist_( 0.0 )
      {};

      Hyperplane ( const Hyperplane &other )
       : normal_( other.normal() ), dist_( other.distance() )
      {};

      Hyperplane ( const Coordinate& normal, const ctype dist )
       : normal_( normal ), dist_( dist )
      {}

      Hyperplane ( const Coordinate& normal, const Coordinate& point )
       : normal_( normal ), dist_( -1.0*( normal * point ) )
      {}

      const Coordinate& normal () const { return normal_; }
      Coordinate& normal () { return normal_; }

      const ctype& distance () const { return dist_; }
      ctype& distance () { return dist_; }

    private:
      Coordinate normal_;
      ctype dist_;
    };

  } // end of namespace VoF

} // end of namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_HYPERPLANE_HH
