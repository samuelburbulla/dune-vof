#ifndef DUNE_VOF_FLAGSET_HH
#define DUNE_VOF_FLAGSET_HH

#include <algorithm>
#include <cmath>
#include <utility>

#include <dune/vof/dataset.hh>

namespace Dune
{
  namespace VoF
  {

    // Flag
    // ----

    enum class Flag : std::size_t {
      empty       = 0,
      mixed       = 2,
      full        = 5,
      mixedfull   = 3,
      activeempty = 1,
      activefull  = 4,
      nan         = 99
    };


    // FlagSet
    // -------

    /**
     * \ingroup Method
     * \brief set of flags
     *
     * \tparam  GridView  grid view
     */
    template< class GridView >
    class FlagSet : public DataSet< GridView, Flag >
    {
    private:
      using This = FlagSet< GridView >;
      using Base =  DataSet< GridView, Flag >;

      using Entity = typename Base::Entity;

      template< std::size_t Lower, std::size_t Upper >
      struct Range
      {
        constexpr Range () = default;
        static constexpr bool contains ( std::size_t val ) { return ( val - Lower <= Upper - Lower ); }

        template< class T, class = std::enable_if_t< std::is_enum< T >{} >,
                           class = std::enable_if_t< std::is_integral< std::underlying_type_t< T > >{} > >
        static constexpr bool contains ( T val )
        {
          return contains( static_cast< std::underlying_type_t< T > >( val ) );
        }
      };

    public:
      FlagSet ( GridView gridView )
      : Base( gridView )
      {}

      using Empty = Range< static_cast< std::size_t >( Flag::empty ),
                           static_cast< std::size_t >( Flag::activeempty ) >;

      using Mixed = Range< static_cast< std::size_t >( Flag::mixed ),
                           static_cast< std::size_t >( Flag::mixedfull ) >;

      using Full = Range< static_cast< std::size_t >( Flag::activefull ),
                          static_cast< std::size_t >( Flag::full ) >;

      using Active = Range< static_cast< std::size_t >( Flag::activeempty ),
                            static_cast< std::size_t >( Flag::activefull ) >;

      bool isEmpty  ( const Entity& entity ) const { return inRange( entity, Empty{} ); }
      bool isMixed  ( const Entity& entity ) const { return inRange( entity, Mixed{} ); }
      bool isFull   ( const Entity& entity ) const { return inRange( entity, Full{} ); }
      bool isActive ( const Entity& entity ) const { return inRange( entity, Active{} );; }

    private:
      template< class _Range >
      bool inRange ( const Entity& en, _Range = {} ) const
      {
        return _Range::contains( this->operator[]( en ) );
      }
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGSET_HH
