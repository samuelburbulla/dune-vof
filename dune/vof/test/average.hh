# ifndef DUNE_VOF_AVERAGE_HH
# define DUNE_VOF_AVERAGE_HH

#include <iostream>

#include <dune/common/exceptions.hh>

#include "../interpolation/interpolations.hh"
#include "problems/problems.hh"
#include "../geometry/intersect.hh"
#include "../geometry/halfspace.hh"

namespace Dune
{
  namespace VoF
  {

    // Average
    // -------

    /**
     * \ingroup   Method
     * \brief     volume of fluid intial data interpolation method
     *
     * \tparam  PR  problem type
     */

    template< class PR >
    struct Average
    {
      using Problem = PR;

      Average ( const Problem& problem ) : problem_( problem ) {}

      template< class ColorFunction >
      void operator() ( ColorFunction& uh, const double time = 0.0, bool verbose = false )
      {
        using GridView = typename std::remove_reference< decltype( uh.gridView() ) >::type;

        #if GridView == SPGrid
          using Interpolation = RecursiveInterpolationCube< GridView >;
        #else
          using Interpolation = RecursiveInterpolation< GridView >;
        #endif

        if( verbose )
          std::cout << " -- average using recursive algorithm" << std::endl;

        auto indicator = [ this, time ]( const auto &x ) { FieldVector< double, 1 > u; this->problem_.evaluate( x, time, u ); return u; };

        Interpolation interpolation ( uh.gridView(), 3 );
        interpolation( indicator, uh );

        uh.communicate();
      }

      const Problem& problem_;
    };


    template<>
    struct Average< Ellipse< double, 2 > >
    {
      using Problem = Ellipse< double, 2 >;

      Average ( const Problem& problem ) : problem_( problem ) {}

      template< class ColorFunction >
      void operator() ( ColorFunction& uh, const double time = 0.0, bool verbose = false )
      {
        if( verbose )
          std::cout << " -- average using intersection algorithm" << std::endl;

        circleInterpolation( problem_.referenceMap(), problem_.volumeElement(), uh );

        uh.communicate();
      }

      const Problem& problem_;
    };

    template<>
    struct Average< RotatingCircle< double, 2 > >
    {
      using Problem = RotatingCircle< double, 2 >;

      Average ( const Problem& problem ) : problem_( problem ) {}

      template< class ColorFunction >
      void operator() ( ColorFunction& uh, const double time = 0.0, bool verbose = false )
      {
        if( verbose )
          std::cout << " -- average using intersection algorithm" << std::endl;

        circleInterpolation( problem_.center( time ), problem_.radius( time ), uh );

        uh.communicate();
      }

      const Problem& problem_;
    };

    template<>
    struct Average< Slope< double, 2 > >
    {
      using Problem = Slope< double, 2 >;

      Average ( const Problem& problem ) : problem_( problem ) {}

      template< class ColorFunction >
      void operator() ( ColorFunction& uh, const double time = 0.0, bool verbose = false )
      {
        if( verbose )
          std::cout << " -- average using intersection algorithm" << std::endl;

        using Coordinate = FieldVector< double, 2 >;

        HalfSpace< Coordinate > halfspace ( problem_.normal(), Coordinate( { 0.5, 0.5 } ) );

        for ( const auto& entity : elements( uh.gridView(), Partitions::interior ) )
        {
          const auto& geo = entity.geometry();
          uh[ entity ] = intersect( makePolytope( geo ), halfspace, eager ).volume() / geo.volume();
        }

        uh.communicate();
      }

      const Problem& problem_;
    };


    template<>
    struct Average< LinearWall< double, 1 > >
    {
      using Problem = LinearWall< double, 1 >;

      Average ( const Problem& problem ) : problem_( problem ) {}

      template< class ColorFunction >
      void operator() ( ColorFunction& uh, const double time = 0.0, bool verbose = false )
      {
        if( verbose )
          std::cout << " -- average using intersection algorithm" << std::endl;

        for ( const auto& entity : elements( uh.gridView(), Partitions::interior ) )
        {
          typename Problem::RangeType jump = problem_.jump( time );

          const auto& geo = entity.geometry();
          if ( geo.corner(0) < jump && geo.corner(1) < jump )
            uh[ entity ] = 1;
          else if ( geo.corner(0) > jump && geo.corner(1) > jump )
            uh[ entity ] = 0;
          else
            uh[ entity ] = ( jump - geo.corner(0) ) / geo.volume();
        }

        uh.communicate();
      }

      const Problem& problem_;
    };

  } // namespace VoF

} // namespace Dune

# endif // #ifndef DUNE_VOF_AVERAGE_HH
