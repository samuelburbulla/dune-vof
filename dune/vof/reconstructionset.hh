#ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
#define DUNE_VOF_RECONSTRUCTIONSET_HH

#include <dune/vof/dataset.hh>
#include <dune/vof/geometry/halfspace.hh>


namespace Dune
{

  namespace VoF
  {

    // ReconstructionSet
    // -----------------

    /**
     * \ingroup Other
     * \brief set of reconstructions
     *
     * \tparam  GridView  grid view
     */
    template< class GridView >
    using ReconstructionSet = DataSet< GridView, HalfSpace< typename GridView::template Codim< 0 >::Geometry::GlobalCoordinate > >;

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
