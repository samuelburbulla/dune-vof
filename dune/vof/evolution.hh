#ifndef DUNE_VOF_EVOLUTION_HH
#define DUNE_VOF_EVOLUTION_HH

#include <functional>
#include <type_traits>

//- local includes
#include <dune/vof/evolution/evolution.hh>
#include <dune/vof/evolution/characteristicsevolution.hh>

namespace Dune
{
  namespace VoF
  {
    // evolution
    // --------

    /**
     * \ingroup Method
     * \brief generate time evolution operator
     *
     * \tparam  GridView
     * \return [description]
     */
    template< class GridView >
    static inline auto evolution ( const GridView& gv )
     -> decltype( Evolution< GridView >( gv ) )
    {
      return Evolution< GridView >( gv );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EVOLUTION_HH
