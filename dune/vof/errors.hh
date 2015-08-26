#ifndef __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__
#define __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune
{
  namespace VoF
  {

    template< class G, class F >
    double l1error ( const G &grid, const std::vector< double > &concentration, const F &ft )
    {
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;

      double error = 0;

      GridView gridView = grid.leafGridView();

      for( auto&& entity : elements( gridView ) )
      {
        const Geometry geo = entity.geometry();

        int i = gridView.indexSet().index( entity );


        const Dune::GeometryType gt = geo.type();
        int p = 20;
        const Dune::QuadratureRule< double, 2 > &rule = Dune::QuadratureRules< double, 2 >::rule( gt, p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, " order   not   available " );

        // compute approximate integral
        double result = 0;
        for( typename Dune::QuadratureRule< double, 2 >::const_iterator qp = rule.begin(); qp != rule.end(); ++qp )
        {
          double fval = std::abs ( ft( geo.global( qp->position() ) ) - concentration[ i ] );
          double weight = qp->weight();
          double detjac = geo.integrationElement( qp->position() );
          result += fval * weight * detjac;
        }

        error += result;
      }

      return error;
    }


    template< class G, class F >
    double l2error ( const G &grid, const std::vector< double > &concentration, const F &ft )
    {
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim< 0 >::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;

      double error = 0;

      GridView gridView = grid.leafGridView();

      for( auto&& entity : elements( gridView ) )
      {
        const Geometry geo = entity.geometry();

        int i = gridView.indexSet().index( entity );

        const Dune::GeometryType gt = geo.type();
        int p = 20;
        const Dune::QuadratureRule< double, 2 > &rule = Dune::QuadratureRules< double, 2 >::rule( gt, p );

        // ensure that rule has at least the requested order
        if( rule.order() < p )
          DUNE_THROW( Dune::Exception, " order   not   available " );

        // compute approximate integral
        double result = 0;
        for( typename Dune::QuadratureRule< double, 2 >::const_iterator qp = rule.begin(); qp != rule.end(); ++qp )
        {
          double fval = std::abs ( ft( geo.global( qp->position() ) ) - concentration[ i ] );
          fval *= fval;
          double weight = qp->weight();
          double detjac = geo.integrationElement( qp->position() );
          result += fval * weight * detjac;
        }

        error += result;
      }

      return std::sqrt( error );
    }


  } // end of namespace VoF
} // end of namespace Dune

#endif