#ifndef DUNE_VOF_CURVATURE_HH
#define DUNE_VOF_CURVATURE_HH

#include <dune/vof/curvature/cartesianheightfunctioncurvature.hh>
#include <dune/vof/curvature/swartzcurvature.hh>

namespace Dune
{

  namespace VoF
  {

    template< class Stencils >
    static inline auto curvature ( Stencils& stencils ) -> SwartzCurvature< typename Stencils::GridView, Stencils >
    {
      return SwartzCurvature< typename Stencils::GridView, Stencils >( stencils );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURE_HH
