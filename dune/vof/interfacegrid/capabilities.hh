#ifndef DUNE_VOF_INTERFACEGRID_CAPABILITIES_HH
#define DUNE_VOF_INTERFACEGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/vof/interfacegrid/declaration.hh>

namespace Dune
{

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< class HostGrid >
    struct hasSingleGeometryType< InterfaceGrid< HostGrid > >
    {
      static const bool v = hasSingleGeometryType< HostGrid >::v;
      static const unsigned int topologyId = hasSingleGeometryType< HostGrid >::topologyId;
    };


    template< class HostGrid >
    struct isCartesian< InterfaceGrid< HostGrid > >
    {
      static const bool v = isCartesian< HostGrid >::v;
    };


    template< class HostGrid, int codim >
    struct hasEntity< InterfaceGrid< HostGrid >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim >::v;
    };

    template< class HostGrid, int codim >
    struct canCommunicate< InterfaceGrid< HostGrid >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v;
    };


    template< class HostGrid >
    struct hasBackupRestoreFacilities< InterfaceGrid< HostGrid > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLevelwiseConforming< InterfaceGrid< HostGrid > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLeafwiseConforming< InterfaceGrid< HostGrid > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct threadSafe< InterfaceGrid< HostGrid > >
    {
      static const bool v = false;
    };

    template< class HostGrid >
    struct viewThreadSafe< InterfaceGrid< HostGrid > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_CAPABILITIES_HH
