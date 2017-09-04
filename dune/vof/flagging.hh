#ifndef DUNE_VOF_FLAGGING_HH
#define DUNE_VOF_FLAGGING_HH

#include <algorithm>
#include <cmath>
#include <utility>

#include <dune/grid/common/partitionset.hh>

#include <dune/vof/flagset.hh>

namespace Dune
{
  namespace VoF
  {

      /**
       * \brief update set of flags
       *
       * \tparam  GV  grid view
       * \param eps   marker tolerance
       */

    template< class GV >
    struct FlagOperator
    {
      using GridView = GV;

    public:
      explicit FlagOperator ( double eps )
       : eps_( eps )
      {}

      template< class ColorFunction, class FlagSet >
      void operator() ( const ColorFunction& color, FlagSet &flags, bool communicate = true ) const
      {
        // TODO
        //if ( communicate )
        //  partition = Partitions::interiorBorder;

        for ( const auto &entity : elements( color.gridView(), Partitions::all ) )
        {
          Flag &flag = flags[ entity ];
          const auto colorEn = color[ entity ];

          if ( colorEn < eps_ )
            flag = Flag::empty;
          else if ( colorEn <= ( 1 - eps_ ) )
            flag = Flag::mixed;
          else
          {
            if ( std::isnan( colorEn ) )
              flag = Flag::nan;
            else
            {
              flag = Flag::full;

              for ( const auto &intersection : intersections( color.gridView(), entity ) )
                if ( intersection.neighbor() && color[ intersection.outside() ] < eps_ )
                {
                  flag = Flag::mixedfull;
                  break;
                }
            }
          }
        }

        if ( communicate )
          flags.communicate();
      }

    private:
      const double eps_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGGING_HH
