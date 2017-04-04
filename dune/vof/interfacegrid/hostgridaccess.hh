#ifndef DUNE_VOF_INTERFACEGRID_HOSTGRIDACCESS_HH
#define DUNE_VOF_INTERFACEGRID_HOSTGRIDACCESS_HH

#include <type_traits>

#include <dune/grid/common/intersection.hh>

#include <dune/vof/interfacegrid/declaration.hh>
#include <dune/vof/interfacegrid/intersection.hh>

namespace Dune
{

  // HostGridAccess
  // --------------

  template< class Grid >
  struct HostGridAccess;



  /**
   * \class HostGridAccess
   * \brief provides access to host grid objects
   *
   * \tparam  Grid  meta grid, whose host grid shall be accessed
   *
   * \nosubgrouping
   **/
  template< class R >
  struct HostGridAccess< VoF::InterfaceGrid< R > >
  {
    /** \name Exported Types
     * \{
     **/

    typedef VoF::InterfaceGrid< R > Grid;

    static const int dimension = Grid::dimension;

    //! type of HostGrid
    typedef typename R::GridView::Grid HostGrid;

     /** \} **/

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     **/
    template< int codim, class T = void >
    struct Codim;

    template< int codim >
    struct Codim< codim, std::enable_if_t< (codim == 0) > >
    {
      //! type of the InterfaceGrid entity
      typedef Dune::Entity< codim, dimension, const Grid, VoF::InterfaceGridEntity > Entity;

      //! type of the host entity
      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
    };

    template< int codim >
    struct Codim< codim, std::enable_if_t< (codim > 0) && (codim <= dimension) > >
    {};

    /**
     * \brief Get underlying HostGrid.
     * \param[in] grid  InterfaceGrid
     * \returns HostGrid
     */
    static const HostGrid &hostGrid ( const Grid &grid )
    {
      return grid.reconstruction().gridView().grid();
    }

    template< int codim >
    static const typename Codim< codim >::HostEntity &
    hostEntity ( const Dune::Entity< codim, dimension, const Grid, VoF::InterfaceGridEntity > &entity )
    {
      return Grid::getRealImplementation( entity ).hostElement();
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_HOSTGRIDACCESS_HH
