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

    template< class Reconstruction >
    struct hasSingleGeometryType< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = false;
      static const unsigned int topologyId = ~static_cast< unsigned int >( 0u );
    };


    template< class Reconstruction >
    struct isCartesian< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = false;
    };


    template< class Reconstruction, int codim >
    struct hasEntity< VoF::InterfaceGrid< Reconstruction >, codim >
    {
      static const bool v = true;
    };

    template< class Reconstruction, int codim >
    struct canCommunicate< VoF::InterfaceGrid< Reconstruction >, codim >
    {
      static const bool v = canCommunicate< Reconstruction, 0 >::v;
    };


    template< class Reconstruction >
    struct hasBackupRestoreFacilities< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = false;
    };

    template< class Reconstruction >
    struct isLevelwiseConforming< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = true;
    };

    template< class Reconstruction >
    struct isLeafwiseConforming< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = true;
    };

    template< class Reconstruction >
    struct threadSafe< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = false;
    };

    template< class Reconstruction >
    struct viewThreadSafe< VoF::InterfaceGrid< Reconstruction > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_CAPABILITIES_HH
