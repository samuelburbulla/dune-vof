#ifndef DUNE_VOF_CURVATURESET_HH
#define DUNE_VOF_CURVATURESET_HH

#include <dune/common/typetraits.hh>

#include <dune/vof/dataset.hh>


namespace Dune
{

  namespace VoF
  {

    // CurvatureSet
    // -----------------

    /**
     * \ingroup Other
     * \brief set of curvature values
     *
     * \tparam  GridView  grid view
     */
    template< class GridView >
    using CurvatureSet = DataSet< GridView, real_t< typename GridView::Grid::ctype > >;


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURESET_HH
