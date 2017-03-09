#ifndef DUNE_VOF_RECONSTRUCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HH

#include <dune/vof/reconstruction/heightfunction.hh>
#include <dune/vof/reconstruction/modifiedyoungs.hh>
#include <dune/vof/reconstructionset.hh>

#include <dune/vof/reconstruction/modifiedswartz.hh>

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

    /*
    template< class GridView, class Stencils >
    static inline auto reconstruction ( const GridView&, Stencils &stencils )
     -> typename std::enable_if< GridView::dimension == 2,
                                 ModifiedSwartzReconstruction< GridView, Stencils,
                                  ModifiedYoungsReconstruction< GridView, Stencils > > >::type
    {
      return ModifiedSwartzReconstruction< GridView, Stencils,
                                           ModifiedYoungsReconstruction< GridView, Stencils >
                                         >( stencils );
    }
    */

    template< class GridView, class Stencils >
    static inline auto reconstruction ( const GridView&, Stencils &stencils )
     ->  HeightFunctionReconstruction< GridView, Stencils, ModifiedYoungsReconstruction< GridView, Stencils > >
    {
      return HeightFunctionReconstruction< GridView, Stencils, ModifiedYoungsReconstruction< GridView, Stencils > >( stencils );
    }

  } // namespace VoF

} // namespace Dune


#endif // #ifndef DUNE_VOF_RECONSTRUCTION_HH
