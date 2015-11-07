#ifndef DUNE_VOF_SWARTZMETHOD_HH
#define DUNE_VOF_SWARTZMETHOD_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include "brents.hh"
#include "geometricutility.hh"
#include "hyperplane.hh"


namespace Dune
{
  namespace VoF
  {

    template< class Entity, class ReconstructionSet, class ColorFunction, class Flags, class Domain, class ReconstructionType >
    void SwartzMethod ( const Entity &entity, const ReconstructionSet &guessedNormals, const ColorFunction &colorFunction,
      const Flags &flags, const Domain &domain, ReconstructionType &g )
    {
      using fvector = typename Entity::Geometry::GlobalCoordinate;

      fvector normalOld, centroidLineNormal, sumNormals, centroid1, centroid2;
      double sumCount;
      int count = 0;

      ReconstructionType h;
      g.normal() = guessedNormals[ entity ].normal();

      std::vector< fvector > entityIntersections, neighborIntersections;


      auto geoEntity = entity.geometry();
      do
      {

        sumCount = 0;
        sumNormals = 0;


        computeInterfaceLinePosition( geoEntity, colorFunction[ entity ], g, entityIntersections );

        centroid1 = 0;

        int c = 0;
        for( std::size_t i = 0; i < entityIntersections.size(); ++i )
        {
          centroid1 += entityIntersections[ i ];
          c++;
        }

        if ( c == 0) continue;
        centroid1 *= 1.0 / c;


        for( auto&& neighbor : domain[ entity ] )

          if( flags.isMixed( neighbor ) )
          {

              // nimm diese Zelle nur, wenn die geratene Normal schon in eine aehnliche Richtung zeigt
              h.normal() = guessedNormals[ neighbor ].normal();

              if( g.normal() * h.normal() < 0 )
                continue;


              auto geoNeighbor = neighbor.geometry();

              h.normal() = g.normal();


              computeInterfaceLinePosition( geoNeighbor, colorFunction[ neighbor ], h, neighborIntersections );

              fvector centroid2;

              int c = 0;
              for( std::size_t i = 0; i < neighborIntersections.size(); ++i )
              {
                centroid2 += neighborIntersections[ i ];
                c++;
              }
              if ( c == 0 ) continue;
              centroid2 *= 1.0 / c;


              centroidLineNormal = centroid2;
              centroidLineNormal -= centroid1;
              rotate90degreesCounterClockwise( centroidLineNormal );

              assert( centroidLineNormal.two_norm() > 0 );


              // errechnete Normale muss in eine aehnliche Richtung wie geschaetze Normale auf jeder Zelle zeigen
              if( centroidLineNormal * g.normal() < 0 )
                centroidLineNormal *= -1.0;


              sumNormals.axpy( 1.0 / (geoEntity.center() - centroid2).two_norm(), centroidLineNormal );
              sumCount++;
          }

          normalOld = g.normal();

          if( sumCount == 0 ) continue;

          g.normal() = sumNormals;
          g.normal() *= 1.0 / sumCount;

          assert( g.normal().two_norm() > 1e-14 );

          count++;

      } while( (g.normal() - normalOld).two_norm2() > 1e-8 && count < 30 );            // limit number of loops

    }

    template< class Geometry, class ReconstructionType, class PointList>
    void computeInterfaceLinePosition ( const Geometry &geo, double concentration, ReconstructionType &g, PointList &intersections )
    {
      intersections.clear();

      double pMin = 0, pMax = 0, volume;
      //use bigger range than [0,1] initially
      double volMin = -1;
      double volMax = 2;

      // Initial guess for p
      for( int i = 0; i < geo.corners(); ++i )
      {
        g.distance() = geo.corner( i ) * g.normal();
        g.distance() *= -1.0;

        volume = getVolumeFraction( geo, g );


        if( volume <= volMax && volume >= concentration )
        {
          pMax = g.distance();
          volMax = volume;
        }

        if( volume >= volMin && volume <= concentration )
        {
          pMin = g.distance();
          volMin = volume;
        }
      }

      g.distance() = brentsMethod( [ &geo, &concentration, &g ] ( double p ) -> double {
                              ReconstructionType h( g.normal(), p );
                              return ( getVolumeFraction( geo, h ) - concentration );
                            }, pMin, pMax, 1e-12 );

      intersections = lineCellIntersections( geo, g );
    }


  } // namespace VoF

} // namespace Dune

#endif
