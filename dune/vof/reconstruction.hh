#ifndef DUNE_VOF_RECONSTRUCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HH

#include <dune/vof/reconstruction/heightfunction.hh>
#include <dune/vof/reconstruction/modifiedswartz.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>
#include <dune/vof/reconstruction/swartz.hh>
#include <dune/vof/reconstructionset.hh>

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
     * \param   gv        grid view
     * \param   stencils  set of stencils
     * \return [description]
     */


    template< class Stencils >
    static inline auto reconstruction ( Stencils &stencils )
     -> HeightFunctionReconstruction< typename Stencils::GridView, Stencils,
                                      ModifiedYoungsReconstruction< typename Stencils::GridView, Stencils > >
    {
      return HeightFunctionReconstruction< typename Stencils::GridView, Stencils,
                                           ModifiedYoungsReconstruction< typename Stencils::GridView, Stencils >
                                         >( stencils );
    }

    /*
    template< class GridView, class Stencils >
    static inline auto reconstruction ( const GridView&, Stencils &stencils )
     ->  ModifiedYoungsReconstruction< GridView, Stencils >
    {
      return ModifiedYoungsReconstruction< GridView, Stencils >( stencils );
    }
    */
  } // namespace VoF

} // namespace Dune


#endif // #ifndef DUNE_VOF_RECONSTRUCTION_HH
