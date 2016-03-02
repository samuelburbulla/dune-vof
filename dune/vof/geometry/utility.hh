#ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
#define DUNE_VOF_GEOMETRY_UTILITY_HH

#include <cassert>

namespace Dune
{
  namespace VoF
  {


    template< class DomainVector >
    inline static DomainVector rotateCCW ( const DomainVector &v )
    {
      return DomainVector{ -v[ 1 ], v[ 0 ] };
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
