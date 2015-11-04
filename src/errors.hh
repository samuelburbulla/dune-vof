#ifndef __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__
#define __DUNE_GRID_REC_VOL_TRACK_ERRORS_HH__

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

#include "L1projection.hh"
#include "geometricaltoolbox.hh"

namespace Dune
{
  namespace VoF
  {

    template< class GridView, class ColorFunction, class F >
    double l1error ( const GridView &gridView, const ColorFunction &colorFunction, const F &ft )
    {
      ColorFunction colorFunctionEnd( gridView );
      Dune::VoF::L1projection( colorFunctionEnd, ft );

      double error = 0;

      for( auto&& entity : elements( gridView ) )
      {
        const auto geo = entity.geometry();

        error += geo.volume() * std::abs( colorFunction[ entity ] - colorFunctionEnd[ entity ] );
      }

      return error;
    }


    template< class GridView, class ColorFunction, class F >
    double l2error ( const GridView &gridView, const ColorFunction &colorFunction, const F &ft )
    {
      ColorFunction colorFunctionEnd( gridView );
      Dune::VoF::L1projection( colorFunctionEnd, ft );

      double error = 0;

      for( auto&& entity : elements( gridView ) )
      {
        const auto geo = entity.geometry();

        double diff = colorFunction[ entity ] - colorFunctionEnd[ entity ];
        diff *= diff;
        error += geo.volume() * diff;
      }

      return std::sqrt( error );
    }


/*
    template< class GridView, class R >
    double recError ( const GridView &grid, const R &rec )
    {
      typedef typename G::LeafGridView GridView;
      typedef typename Dune::FieldVector<double,2> V;

      V center { 0.47, 0.54 };

      double error = 0;

      GridView gridView = grid.leafGridView();

      for( auto&& entity : elements( gridView ) )
      {
        int i = gridView.indexSet().index( entity );
        double result = 0;

        if ( rec[i][2].two_norm() > 0.01 )
        {
          V discreteNormal = rec[i][2];
          discreteNormal /= discreteNormal.two_norm();
          
          
          V exactNormalc = rec[i][0];
          exactNormalc += rec[i][1]; 
          exactNormalc *= 0.5;
          exactNormalc -= center;
          exactNormalc /= exactNormalc.two_norm();

          V exactNormal0 = rec[i][0];
          exactNormal0 -= center;
          exactNormal0 /= exactNormal0.two_norm();

          V exactNormal1 = rec[i][1];
          exactNormal1 -= center;
          exactNormal1 /= exactNormal1.two_norm();
          


          result += 0.5 * (discreteNormal - exactNormalc).one_norm();
          result += 0.25 * (discreteNormal - exactNormal0).one_norm();
          result += 0.25 * (discreteNormal - exactNormal1).one_norm(); 

          result *= (rec[i][0] - rec[i][1]).two_norm();

          error += result;
        }
      }

      return error;
    }
*/

  } // end of namespace VoF
} // end of namespace Dune

#endif