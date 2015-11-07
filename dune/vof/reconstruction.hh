#ifndef DUNE_VOF_RECONSTRUCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HH

#include <dune/vof/reconstruction/modifiedswartz.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>
#include <dune/vof/reconstructionSet.hh>

namespace Dune
{
  namespace VoF
  {

    template< class GridView, class Stencils, class ColorFunction >
    auto reconstruction ( const GridView&, const ColorFunction&, Stencils &stencils )
     -> decltype( ModifiedSwartzReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils,
                                                ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils >
                                              >( stencils ) )
    {
      return ModifiedSwartzReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils,
                                           ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils >
                                         >( stencils );
    }

  } // namespace VoF

} // namespace Dune


#endif // #ifndef DUNE_VOF_RECONSTRUCTION_HH
