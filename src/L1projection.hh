# ifndef __DUNE__REC_VOL_TRACK_L1PROJECTION_HH__
# define __DUNE__REC_VOL_TRACK_L1PROJECTION_HH__

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune
{
  namespace VoF
  {


    template< class GridView, template <class> class ColorFunction, class ContinuousFunction >
    void L1projection ( ColorFunction<GridView> &colorFunction, const ContinuousFunction &f )
    {
      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::ctype ctype;

      for( auto&& entity : elements ( colorFunction.gridView() ) )
      {
        const auto geo = entity.geometry();
        const Dune::GeometryType gt = geo.type();

        // get quadrature rule of order p
        int p = 10;
        const Dune::QuadratureRule< ctype, dimworld > &rule 
         = Dune::QuadratureRules< ctype, dimworld >::rule( gt, p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, " order   not   available " );

        // compute approximate integral
        ctype result = 0;

        for( typename Dune::QuadratureRule< ctype, dimworld >::const_iterator i = rule.begin(); i != rule.end(); ++i )
        {
          ctype fval = f( geo.global( i->position() ) );
          ctype weight = i->weight();
          ctype detjac = geo.integrationElement( i->position() );
          result += fval * weight * detjac;
        }

        colorFunction[ entity ] = result / geo.volume();
      }

    }

  }
}

# endif