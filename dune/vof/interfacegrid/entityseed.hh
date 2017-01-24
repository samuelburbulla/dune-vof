#ifndef DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH
#define DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH

#include <type_traits>

#include <dune/grid/common/entityseed.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridEntitySeed
    // -----------------------

    template< int codim, class Grid >
    class InterfaceGridEntitySeed
    {
      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int codimension = codim;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;

      typedef typename Traits::Reconstruction::GridView::Grid::template Codim< 0 >::EntitySeed HostElementSeed;

      InterfaceGridEntitySeed () = default;

      InterfaceGridEntitySeed ( const HostElementSeed &hostElementSeed, int subEntity )
        : hostElementSeed_( hostElementSeed ), subEntity_( subEntity )
      {}

      bool isValid () const { return hostElementSeed_.isValid(); }

      const HostElementSeed &hostElementSeed () const { return hostElementSeed_; }
      int subEntity () const { return subEntity_; }

    private:
      HostElementSeed hostElementSeed_;
      int subEntity_;
    };



    // InterfaceGridEntitySeed for codimension 0
    // -----------------------------------------

    template< class Grid >
    class InterfaceGridEntitySeed< 0, Grid >
    {
      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      static const int codimension = 0;
      static const int dimension = Traits::dimension;
      static const int mydimension = dimension - codimension;

      typedef typename Traits::Reconstruction::GridView::Grid::template Codim< 0 >::EntitySeed HostElementSeed;

      InterfaceGridEntitySeed () = default;

      explicit InterfaceGridEntitySeed ( const HostElementSeed &hostElementSeed )
        : hostElementSeed_( hostElementSeed )
      {}

      bool isValid () const { return hostElementSeed_.isValid(); }

      const HostElementSeed &hostElementSeed () const { return hostElementSeed_; }

    private:
      HostElementSeed hostElementSeed_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_ENTITYSEED_HH
