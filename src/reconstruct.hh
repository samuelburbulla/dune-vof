#ifndef __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__
#define __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__

#include <stdexcept>
#include <stdio.h>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>

#include "geometricaltoolbox.hh"
#include "secondyoungsnormalguessing.hh"
#include "swartzmethod.hh"
#include "hypersurface.hh"

namespace Dune
{
  namespace VoF
  {

    template< class GV, class E, class V >
    double computeInterfaceLinePosition( const GV &gridView, const E &entity, const V &n, double concentration );

    template< class GV, class E, class V >
    V interfaceLineCentroid( const GV &gridView, const E &entity, const V &normal, const double pEntity );




    template< class GridView, class ColorFunction, class ReconstructionSet, class Flags, class Domain >
    void reconstruct ( const GridView &gridView, const ColorFunction &colorFunction, 
                       ReconstructionSet  &reconstructionSet, 
                       const Domain &domain, const Flags &flags, const double eps )
    {

      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::ctype ctype;
      typedef typename Dune::FieldVector< ctype, dimworld > fvector;
      typedef typename Dune::VoF::HyperSurface< fvector > ReconstructionType; 


      reconstructionSet.clear();

      ReconstructionSet guessedNormals ( gridView );
      Dune::VoF::SecondYoungsNormalGuessing ( gridView, colorFunction, flags, domain, guessedNormals );



      for( auto&& entity : elements( gridView ) )
      {

        if( flags.isMixed( entity ) )
        {
          ReconstructionType improvedRec;
          Dune::VoF::SwartzMethod ( gridView, entity, guessedNormals, colorFunction, flags, domain, improvedRec );
           
          reconstructionSet[ entity ] = improvedRec;
          Dune::VoF::computeInterfaceLinePosition( gridView, entity, entity.geometry(), colorFunction[ entity ], improvedRec, 
            reconstructionSet.intersections( entity ) ); 
        }

        // reconstructions, which are given by edges of elements
        else if( flags.isFullAndMixed( entity ) )
        {
          for( auto&& intersection : intersections( gridView, entity ) )
          {
            if ( intersection.neighbor() )
            {
              auto neighbor = intersection.outside();
              
              if ( colorFunction[ neighbor ] < eps )               
              {
                auto isGeo = intersection.geometry();
                auto n = intersection.centerUnitOuterNormal();
                n *= -1.0;

                auto p = isGeo.corner(0) * n;
                p *= -1.0;
          
                reconstructionSet[ entity ] = HyperSurface< fvector > ( n, p );
                reconstructionSet.intersections( entity ) = std::vector< fvector > ( { isGeo.corner(0), isGeo.corner(1) } );
              }
            }
          }
        } 
        

      }



    }






  } // end of namespace VoF
} // end of namespace Dune


#endif

