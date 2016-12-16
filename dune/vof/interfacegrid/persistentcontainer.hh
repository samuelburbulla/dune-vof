#ifndef DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH
#define DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainerwrapper.hh>

#include <dune/vof/interfacegrid/declaration.hh>

namespace Dune
{

  // PersistentContainer for InterfaceGrid
  // ------------------------------

  template< class HostGrid, class T >
  class PersistentContainer< InterfaceGrid< HostGrid >, T >
  : public PersistentContainerWrapper< InterfaceGrid< HostGrid >, T >
  {
    typedef PersistentContainerWrapper< InterfaceGrid< HostGrid >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid, codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH
