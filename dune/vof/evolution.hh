#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include "geometricutility.hh"


namespace Dune
{
  namespace VoF
  {


    template< class GridView, class ColorFunction, class ReconstructionSet, class Flags, class VelocityField >
    void evolve ( const GridView &gridView, const ColorFunction &colorFunction, ReconstructionSet &reconstructionSet,
                  const double t, const double dt, const Flags &flags, const VelocityField &velocityField, const double eps, ColorFunction &update )
    {

      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::ctype ctype;
      typedef typename Dune::FieldVector< ctype, dimworld > fvector;

      update.clear();


      for( const auto &entity : elements( gridView ) )
      {
        if( flags.isMixed( entity ) || flags.isActive( entity ) )
        {

          auto entityGeo = entity.geometry();

          for( const auto &intersection : intersections( gridView, entity ) )
          {
            if( intersection.neighbor() )
            {

              const auto &neighbor = intersection.outside();
              const auto isGeo = intersection.geometry();

              fvector outerNormal = intersection.centerUnitOuterNormal();
              fvector velocity = velocityField( isGeo.center() );

                // build time integration polygon
                velocity *= dt;

                Polygon2D< fvector > timeIntegrationPolygon;

                timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) );
                timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) );

                timeIntegrationPolygon.addVertex( isGeo.corner( 0 ) - velocity );
                timeIntegrationPolygon.addVertex( isGeo.corner( 1 ) - velocity );


                // outflows
                if( velocity * outerNormal > 0 )
                {

                  // build flux polygon with present material in cell
                  Polygon2D< fvector > fluxPolygon;
                  double fluxVolume = 0;

                  if( flags.isMixed( entity ) )
                  {
                    polygonLineIntersection( timeIntegrationPolygon, reconstructionSet[ entity ], fluxPolygon );
                    polyAddInnerVertices( timeIntegrationPolygon, reconstructionSet[ entity ], fluxPolygon );

                    fluxVolume = fluxPolygon.volume();
                  }

                  else if( colorFunction[ entity ] >= 1 - eps )
                    fluxVolume = timeIntegrationPolygon.volume();


                  update[ entity ] -= fluxVolume / entityGeo.volume();
                }

                //inflow
                else if( velocity * outerNormal < 0 )
                {
                  // build phase polygon of the neighbor
                  Polygon2D< fvector > fluxPolygon;
                  double fluxVolume = 0;

                  if( flags.isMixed( neighbor ) )
                  {
                    polygonLineIntersection( timeIntegrationPolygon, reconstructionSet[ neighbor ], fluxPolygon );
                    polyAddInnerVertices( timeIntegrationPolygon, reconstructionSet[ neighbor ], fluxPolygon );

                    fluxVolume = fluxPolygon.volume();
                  }
                  else if( colorFunction[ neighbor ] > 1 - eps )
                    fluxVolume = timeIntegrationPolygon.volume();

                  update[ entity ] += fluxVolume / entityGeo.volume();

                }
            } // end of intersection is neighbor
          }
        }
      }
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
