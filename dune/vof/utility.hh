#ifndef DUNE_VOF_UTILITY_HH
#define DUNE_VOF_UTILITY_HH

// dune-common includes
#include <dune/common/fmatrix.hh>

// dune-vof includes
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/polytope.hh>
#include <dune/vof/geometry/intersect.hh>

namespace Dune
{
  namespace VoF
  {

    // Generate input for interface grid
    // ---------------------------------
    template< class ReconstructionSet, class Flags >
    void getInterfaceVertices(
      const ReconstructionSet &reconstructions,
      const Flags &flags,
      std::vector< typename ReconstructionSet::DataType::Coordinate > &vertices,
      std::vector< std::size_t > &offsets
    )
    {
      vertices.clear();
      offsets.clear();

      std::size_t offset = 0;
      offsets.push_back( 0 );

      for( const auto &entity : elements( reconstructions.gridView() ) )
      {
        if ( !flags.isMixed( entity ) )
          continue;

        auto segment = interface( entity, reconstructions );

        assert( segment.size() > 0 );
        assert( segment.vertex( 0 ) != segment.vertex( 1 ) );

        offset += segment.size();
        offsets.push_back( offset );

        for ( std::size_t i = 0; i < segment.size(); ++i )
          vertices.push_back( segment.vertex( i ) );
      }
    }


    template< class ctype >
    ctype clamp( ctype v, ctype lo, ctype hi )
    {
      return std::max( std::min( v, hi ), lo );
    }


  } // namespace VoF

}  // namespace Dune

#endif // #ifndef DUNE_VOF_UTILITY_HH
