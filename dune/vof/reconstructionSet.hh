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
    using ReconstructionSet = DataSet< GridView, HalfSpace< typename __impl::Entity_t< GridView >::Geometry::GlobalCoordinate > >;

    // template< class GridView >
    // struct ReconstructionSet
    //   : public DataSet< GridView, HalfSpace< typename __impl::Entity_t< GridView >::Geometry::GlobalCoordinate > >
    // {
    //   using Reconstruction = HalfSpace< typename decltype(std::declval< GridView >().template begin< 0 >())::Entity::Geometry::GlobalCoordinate >;
    //   using BaseType = DataSet< GridView, Reconstruction >;

    //   explicit ReconstructionSet ( const GridView &gridView )
    //     : BaseType( gridView )
    //   {}
    // };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
