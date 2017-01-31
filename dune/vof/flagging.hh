#ifndef DUNE_VOF_FLAGGING_HH
#define DUNE_VOF_FLAGGING_HH

#include <algorithm>
#include <cmath>
#include <utility>

#include <dune/vof/flagset.hh>

namespace Dune
{
  namespace VoF
  {

      /**
       * \brief update set of flags
       *
       * \tparam  DF  discrete function type
       * \param color discrete function
       * \param eps   marker tolerance
       */

    template< class DF, class FS >
    struct FlagOperator
    {
      using ColorFunction = DF;
      using FlagSet = FS;
      using GridView = typename ColorFunction::GridView;

    public:
      explicit FlagOperator ( double eps )
       : eps_( eps )
      {}

      void operator() ( const ColorFunction& color, FlagSet &flags ) const
      {
        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
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
        flags.communicate();

        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( flags[ entity ] == Flag::mixed || flags[ entity ] == Flag::mixedfull )
            for ( const auto &intersection : intersections( color.gridView(), entity ) )
              if ( intersection.neighbor() )
              {
                Flag &flag = flags[ intersection.outside() ];

                if ( flag == Flag::empty )
                  flag = Flag::activeempty;
                else if ( flag == Flag::full )
                  flag = Flag::activefull;
              }
        }
        flags.communicate( Dune::All_All_Interface, [] ( Flag a, Flag b  ) { if ( a == Flag::activefull || b == Flag::activefull ) return Flag::activefull; else return std::max( a, b ); } );
      }

    private:
      const double eps_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGGING_HH
