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
    using CurvatureSet = DataSet< GridView, real_t< typename __impl::Entity_t< GridView >::Geometry::ctype > >;

    // template< class GridView >
    // struct CurvatureSet
    //  : public DataSet< GridView, double >
    // {
    //   using ThisType = CurvatureSet< GridView >
    //   using BaseType = DataSet< GridView, double >;

    //   explicit CurvatureSet ( const GridView &gridView )
    //     : BaseType( gridView )
    //   {}
    // };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURESET_HH
