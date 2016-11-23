# ifndef DUNE_VOF_AVERAGE_HH
# define DUNE_VOF_AVERAGE_HH

#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

namespace Dune
{
  namespace VoF
  {

    // Average
    // -------

    template< class DF, class F >
    void average ( DF &u, const F &f, const double time = 0.0 )
    {
      typedef typename DF::GridView::ctype ctype;
      using RangeType = Dune::FieldVector< ctype, 1 >;


      for( const auto& entity : elements ( u.gridView() ) )
      {
        const auto geo = entity.geometry();

        // get quadrature rule of order p
        int p = 19;
        const auto &rule = Dune::QuadratureRules< ctype, DF::GridView::dimension >::rule( geo.type(), p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, "Requested quadrature order not available." );

        // compute approximate integral
        ctype result = 0;
        for ( const auto qp : rule )
        {
          RangeType u;
          f.evaluate( geo.global( qp.position() ), time, u );
          result += u * qp.weight() * geo.integrationElement( qp.position() );
        }

        u[ entity ] = result / geo.volume();
      }
    }


  } // namespace VoF

} // namespace Dune

# endif // #ifndef DUNE_VOF_AVERAGE_HH
