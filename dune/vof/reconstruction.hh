#ifndef DUNE_VOF_RECONSTRUCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HH

#include <dune/vof/reconstruction/modifiedswartz.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>
#include <dune/vof/reconstructionSet.hh>

namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Method
     * \brief   generate modified Swartz reconstruction operator with a modifed Youngs reconstruction as inital guess
     * \details \see Reconstruction
     *
     * \tparam  GridView
     * \tparam  Stencils
     * \tparam  ColorFunction
     * \param   gv        grid view
     * \param   cf        discrete function
     * \param   stencils  set of stencils
     * \return [description]
     */
    /*
    template< class GridView, class Stencils, class ColorFunction >
    static inline auto reconstruction ( const GridView&, const ColorFunction&, Stencils &stencils )
     -> ModifiedSwartzReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils,
                                      ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils > >
    {
      return ModifiedSwartzReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils,
                                           ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils >
                                         >( stencils );
    }
    */
    template< class GridView, class Stencils, class ColorFunction >
    static inline auto reconstruction ( const GridView&, const ColorFunction&, Stencils &stencils )
     -> ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils >
    {
      return ModifiedYoungsReconstruction< ColorFunction, ReconstructionSet< GridView >, Stencils >( stencils );
    }


  } // namespace VoF

} // namespace Dune


#endif // #ifndef DUNE_VOF_RECONSTRUCTION_HH
