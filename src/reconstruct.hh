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
#include "brent.hh"
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
                       const Domain &domain, const Flags &flags )
    {

      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::ctype ctype;
      typedef typename Dune::FieldVector< ctype, dimworld > fvector;
      typedef typename Dune::VoF::HyperSurface< fvector > ReconstructionType; 


      reconstructionSet.clear();

      ReconstructionSet guessedNormals ( gridView );
      Dune::VoF::SecondYoungsNormalGuessing ( gridView, colorFunction, flags, domain, guessedNormals );


      int count = 0;


      fvector normalOld;
      ReconstructionType g;

      for( auto&& entity : elements( gridView ) )
      {
        count = 0;

        if( flags.isMixed( entity ) )
        {
          auto geoEntity = entity.geometry();

          g.normal() = guessedNormals[ entity ].normal();


          double sumCount;
          fvector centroidLineNormal, sumNormals;

          do
          {

            sumCount = 0;
            sumNormals = 0;

            
            std::vector< fvector > entityIntersections;
            computeInterfaceLinePosition( gridView, entity, geoEntity, colorFunction[ entity ], g, entityIntersections );


            fvector centroid1;
            
            if( entityIntersections.size() >= 2 )
            {
              int c = 0;
              for( std::size_t i = 0; i < entityIntersections.size(); ++i )
              {
                centroid1 += entityIntersections[ i ];
                c++;
              }
              centroid1 *= 1.0 / c;
            }
            else continue;


            for( auto&& neighbor : domain[ entity ] )

              if( flags.isMixed( neighbor ) )
              {

                ReconstructionType h;

                // nimm diese Zelle nur, wenn die geratene Normal schon in eine aehnliche Richtung zeigt
                h.normal() = guessedNormals[ neighbor ].normal();

                if( g.normal() * h.normal() < 0 )
                  continue;


                auto geoNeighbor = neighbor.geometry();

                h.normal() = g.normal();


                std::vector< fvector > neighborIntersections;
                computeInterfaceLinePosition( gridView, neighbor, geoNeighbor, colorFunction[ neighbor ], h, neighborIntersections );

                fvector centroid2;
                // reconstruction doesn't intersect
                if( neighborIntersections.size() == 0 )
                {
                  continue;
                }
                // get middle of intersections
                else
                {
                  int c = 0;
                  for( std::size_t i = 0; i < neighborIntersections.size(); ++i )
                  {
                    centroid2 += neighborIntersections[ i ];
                    c++;
                  }
                  centroid2 *= 1.0 / c;
                }

                centroidLineNormal = centroid2 - centroid1;
                rotate90degreesCounterClockwise( centroidLineNormal );


                assert( centroidLineNormal.two_norm() > 0 );



                // errechnete Normale muss in eine aehnliche Richtung wie geschaetze Normale auf jeder Zelle zeigen
                if( centroidLineNormal * g.normal() < 0 )
                  centroidLineNormal *= -1.0;


                sumNormals.axpy( 1.0 / (geoEntity.center() - centroid2).two_norm(), centroidLineNormal );
                sumCount++;
              }



            normalOld = g.normal();

            if( sumCount > 0 )
            {
              g.normal() = sumNormals;
              g.normal() *= 1.0 / sumCount;
            }

            assert( g.normal().two_norm() > 1e-14 );

            count++;

          } while( (g.normal() - normalOld).two_norm2() > 1e-8 && count < 30 );            // limit number of loops
          

          std::vector< fvector > entityIntersections;
          computeInterfaceLinePosition( gridView, entity, geoEntity, colorFunction[ entity ], g, entityIntersections );

          if( entityIntersections.size() != 0 )
          {
            reconstructionSet[ entity ] = g;
            reconstructionSet.intersections( entity ) = entityIntersections;
          }

        }
      }

    }




    template< class GridView, class Entity, class Geometry, class ReconstructionType, class PointList>
    void computeInterfaceLinePosition ( const GridView &gridView, const Entity &entity, const Geometry &geo, 
      const double concentration, ReconstructionType &g, PointList &intersections )
    {
      double pMin = 0, pMax = 0, volume;
      //use bigger range than [0,1] initially
      double volMin = -1;
      double volMax = 2;

      // Initial guess for p
      for( int i = 0; i < geo.corners(); ++i )
      {
        g.p() = geo.corner( i ) * g.normal();
        g.p() *= -1.0;

        volume = getVolumeFraction( gridView, entity, geo, g );


        if( volume <= volMax && volume >= concentration )
        {
          pMax = g.p();
          volMax = volume;
        }

        if( volume >= volMin && volume <= concentration )
        {
          pMin = g.p();
          volMin = volume;
        }
      }


      g.p() = Dune::VoF::brentsMethod( pMin, pMax,
                          [ &gridView, &entity, &geo, &concentration, &g ] ( double p ) -> double
                          { ReconstructionType h( g.normal(), p ); return getVolumeFraction( gridView, entity, geo, h ) - concentration; } );

      intersections = lineCellIntersections( gridView, entity, geo, g );
    }



  } // end of namespace VoF
} // end of namespace Dune


#endif

