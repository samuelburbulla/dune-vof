#ifndef DUNE_VOF_CURVATURESET_HH
#define DUNE_VOF_CURVATURESET_HH

#include <algorithm>
#include <vector>

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
    struct CurvatureSet
     : public DataSet< GridView, double >
    {
      using BaseType = DataSet< GridView, double >;

      explicit CurvatureSet ( const GridView &gridView )
        : BaseType( gridView )
      {}
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURESET_HH
