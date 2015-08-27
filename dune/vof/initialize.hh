# ifndef __DUNE__REC_VOL_TRACK_INITIALIZE_HH__
# define __DUNE__REC_VOL_TRACK_INITIALIZE_HH__

#include "testproblem.hh"
#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune
{
  namespace VoF
  {


    template< class G, class V, class F >
    void initialize ( const G &grid, V &c, F &f )
    {
      const int dimworld = G::dimensionworld;

      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity::Geometry Geometry;

      GridView gridView = grid.leafGridView();

      for( LeafIterator it = gridView.template begin< 0 >(); it != gridView.template end< 0 >(); ++it )
      {
        const Geometry geo = it->geometry();

        const Dune::GeometryType gt = geo.type();

        // get quadrature rule of order p
        int p = 20;
        const Dune::QuadratureRule< double, dimworld > &rule = Dune::QuadratureRules< double, dimworld >::rule( gt, p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, " order   not   available " );

        // compute approximate integral
        double result = 0;

        for( typename Dune::QuadratureRule< double, dimworld >::const_iterator i = rule.begin(); i != rule.end(); ++i )
        {
          double fval = f( geo.global( i->position() ) );
          double weight = i->weight();
          double detjac = geo.integrationElement( i->position() );
          result += fval * weight * detjac;
        }

        const int indexi = gridView.indexSet().index( *it );
        c[ indexi ] = result / geo.volume();
      }

    }



    template< class G, class V, class F >
    void L1projection ( const G &grid, V &v, F &f )
    {
      const int dimworld = G::dimensionworld;
      typedef typename G::ctype ctype;
      typedef typename Dune::FieldVector< ctype, dimworld > fvector;
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity::Geometry Geometry;

      GridView gridView = grid.leafGridView();

      for( LeafIterator it = gridView.template begin< 0 >(); it != gridView.template end< 0 >(); ++it )
      {
        const Geometry geo = it->geometry();

        const Dune::GeometryType gt = geo.type();

        // get quadrature rule of order p
        int p = 20;
        const Dune::QuadratureRule< double, dimworld > &rule = Dune::QuadratureRules< double, dimworld >::rule( gt, p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, " order   not   available " );

        // compute approximate integral
        fvector result (0);

        for( typename Dune::QuadratureRule< double, dimworld >::const_iterator i = rule.begin(); i != rule.end(); ++i )
        {
          fvector fval = f( geo.global( i->position() ) );
          
          fval *= i->weight();
          fval *= geo.integrationElement( i->position() );
          result += fval;
        }

        const int indexi = gridView.indexSet().index( *it );
        v[ indexi ] = result;
        v[ indexi ] /= geo.volume();
      }

    }

  }
}

# endif