#ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
#define DUNE_VOF_EDGENEIGHBORSSTENCIL_HH

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // EdgeNeighborsStencil
    // ----------------------

    template< class GV >
    struct EdgeNeighborsStencil
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using EntitySeed = typename Entity::EntitySeed;
      using Stencil = std::vector< Entity >;

    public:
      explicit EdgeNeighborsStencil ( const GridView &gridView ) : gridView_( gridView ) {}

      Stencil operator[] ( const Entity &entity ) const
      {
        Stencil stencil;
        for ( const auto &intersection : intersections( gridView(), entity ) )
          if ( intersection.neighbor() )
          {
            const auto &neighbor = intersection.outside();
            stencil.push_back( neighbor );
          }
        return stencil;
      }

    private:
      const GridView &gridView () const { return gridView_; }

      GridView gridView_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
