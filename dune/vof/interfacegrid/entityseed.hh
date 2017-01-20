#ifndef DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH
#define DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/entityseed.hh>

namespace Dune
{

  namespace VoF
  {

    template< int codim, class Grd >
    class InterfaceGridEntitySeed
    {
      typedef typename remove_const< Grd >::type::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef typename Traits::Grid Grid;
      typedef typename Traits::template Codim< codim >::Entity Entity;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename HostGrid::template Codim< codim >::EntitySeed HostEntitySeed;

      explicit InterfaceGridEntitySeed ( const HostEntitySeed &hostEntitySeed )
        : hostEntitySeed_( hostEntitySeed )
      {}

      InterfaceGridEntitySeed ()
        : hostEntitySeed_()
      {}

      bool isValid() const { return hostEntitySeed_.isValid(); }

      const HostEntitySeed &hostEntitySeed () const { return hostEntitySeed_; }

    private:
      HostEntitySeed hostEntitySeed_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH
