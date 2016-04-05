#ifndef DUNE_VOF_GEOMETRY_ROTATION_HH
#define DUNE_VOF_GEOMETRY_ROTATION_HH

/* dune includes */
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/* local includes */
#include "polyhedron.hh"

namespace Dune
{
  namespace VoF
  {

    template < class Coord >
    const Polyhedron< Coord > rotateToReferenceFrame ( const Coord& n, const Polyhedron< Coord >& cell )
    {
      // If n == e_z, there is nothing to do.
      if ( n == Coord ( { 0.0, 0.0, 1.0 } ) )
        return cell;


      Dune::FieldMatrix< double, 3, 3 > rotationMatrix;

      // If n == - e_z, rotate of angle pi around x-axis.
      if ( n == Coord ( { 0.0, 0.0, -1.0 } ) )
      {
        rotationMatrix[0][0] = 1.0;
        rotationMatrix[1][1] = -1.0;
        rotationMatrix[2][2] = -1.0;
      }
      else
      {
        double gamma = ( 1 - n[2] ) / ( n[0] * n[0] + n[1] * n[1] );
        rotationMatrix[0][0] = n[1] * n[1] * gamma + n[2];
        rotationMatrix[0][1] = - n[0] * n[1] * gamma;
        rotationMatrix[1][0] = rotationMatrix[0][1];
        rotationMatrix[0][2] = - n[0];
        rotationMatrix[1][1] = n[0] * n[0] * gamma + n[2];
        rotationMatrix[1][2] = - n[1];
        rotationMatrix[2][0] = n[0];
        rotationMatrix[2][1] = n[1];
        rotationMatrix[2][2] = n[2];
      }

      std::vector< Coord > newNodes ( cell.nodes().size() );

      for ( std::size_t i = 0; i < cell.nodes().size(); ++i )
        rotationMatrix.mv( cell.node( i ), newNodes[ i ] );

      return Polyhedron< Coord > ( cell, newNodes );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_ROTATION_HH
