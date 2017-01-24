#ifndef DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH
#define DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH

#include <vector>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainervector.hh>

#include <dune/vof/interfacegrid/declaration.hh>
#include <dune/vof/interfacegrid/indexset.hh>

namespace Dune
{

  // PersistentContainer for InterfaceGrid
  // -------------------------------------

  template< class Reconstruction, class T >
  class PersistentContainer< VoF::InterfaceGrid< Reconstruction >, T >
    : public PersistentContainerVector< VoF::InterfaceGrid< Reconstruction >, VoF::InterfaceGridIndexSet< const VoF::InterfaceGrid< Reconstruction > >, std::vector< T > >
  {
    typedef PersistentContainerVector< VoF::InterfaceGrid< Reconstruction >, VoF::InterfaceGridIndexSet< const VoF::InterfaceGrid< Reconstruction > >, std::vector< T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid.leafIndexSet(), codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_PERSISTENTCONTAINER_HH
