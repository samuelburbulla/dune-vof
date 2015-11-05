# ifndef AVERAGE_HH
# define AVERAGE_HH

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

template< class DF, class F >
void average ( DF &u, const F &f )
{
  typedef typename DF::GridView::ctype ctype;

  for( auto&& entity : elements ( u.gridView() ) )
  {
    const auto geo = entity.geometry();

    // get quadrature rule of order p
    int p = 10;
    const auto &rule = Dune::QuadratureRules< ctype, DF::GridView::dimensionworld >::rule( geo.type(), p );

    // ensure that rule has at least the requested order
    if( rule.order() < p )
      DUNE_THROW( Dune::Exception, "Requested quadrature order not available." );

    // compute approximate integral
    ctype result = 0;

    for ( auto i = rule.begin(); i != rule.end(); ++i )
    {
      ctype val = f( geo.global( i->position() ) );
      ctype weight = i->weight() * geo.integrationElement( i->position() );
      result += val * weight;
    }

    u[ entity ] = result / geo.volume();
  }
}

# endif // #ifndef AVERAGE_HH
