#ifndef DUNE_VOF_RECONSTRUCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HH

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

//- dune-grid includes
#include <dune/grid/common/geometry.hh>

//- local includes
#include "brents.hh"
#include "geometricutility.hh"
#include "hyperplane.hh"
#include "secondyoungsnormalguessing.hh"
#include "modifiedswartz.hh"


namespace Dune
{
  namespace VoF
  {

    template< class GridView, class ColorFunction, class ReconstructionSet, class Flags, class Domain >
    void reconstruct ( const GridView &gridView, const ColorFunction &colorFunction,
                       ReconstructionSet  &reconstructionSet,
                       const Domain &domain, const Flags &flags, const double eps )
    {

      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::ctype ctype;
      typedef typename Dune::FieldVector< ctype, dimworld > fvector;
      typedef typename Dune::VoF::Hyperplane< fvector > ReconstructionType;


      reconstructionSet.clear();

      ReconstructionSet guessedNormals ( gridView );
      Dune::VoF::SecondYoungsNormalGuessing ( gridView, colorFunction, flags, domain, guessedNormals );



      for( auto&& entity : elements( gridView ) )
      {

        if( flags.isMixed( entity ) )
        {
          ReconstructionType improvedRec;
          SwartzMethod ( entity, guessedNormals, colorFunction, flags, domain, improvedRec );

          reconstructionSet[ entity ] = improvedRec;
          computeInterfaceLinePosition( entity.geometry(), colorFunction[ entity ], improvedRec, reconstructionSet.intersections( entity ) );
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

                reconstructionSet[ entity ] = Hyperplane< fvector > ( n, p );
                reconstructionSet.intersections( entity ) = std::vector< fvector > ( { isGeo.corner(0), isGeo.corner(1) } );
              }
            }
          }
        }
      }
    }

  } // namespace VoF

} // namespace Dune


#endif // #ifndef DUNE_VOF_RECONSTRUCTION_HH
