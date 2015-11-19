#ifndef DUNE_VOF_FLAGS_HH
#define DUNE_VOF_FLAGS_HH

#include <utility>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // Flag the entities of a gridView as mixed or active

    template< class GV >
    struct Flags
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;

    private:
      enum class Flag {
        empty       = 0,
        mixed       = 1,
        full        = 2,
        mixedfull   = 3,
        activeempty = 4,
        activefull  = 5
      };

      using Mapper = MCMGMapper< GridView, MCMGElementLayout >;

    public:
      Flags ( const GridView &gridView )
       : gridView_ ( gridView ), mapper_( gridView ), flags_( mapper_.size(), Flag::empty )
      {}

      const bool isMixed ( const Entity& entity ) const { return flags_[ index( entity ) ] == Flag::mixed; }
      const bool isFullAndMixed ( const Entity& entity ) const { return flags_[ index( entity ) ] == Flag::mixedfull; }

      const bool isActive ( const Entity& entity ) const
      {
        const Flag &flag = flags_[ index( entity ) ];
        return ( flag == Flag::activeempty ) || ( flag == Flag::activefull ); }

      const std::size_t size() const { return mapper_.size(); }

      const int operator[] ( const int i ) const { return static_cast< int >( flags_[ i ] ); }

      template< class DF >
      void reflag ( const DF& color, const double eps )
      {
        for ( const auto &entity : elements( gridView() ) )
        {
          const auto idx = index( entity );
          Flag &flag = flags_[ idx ];
          const auto colorEn = color[ idx ];

          if ( colorEn < eps )
            flag = Flag::empty;
          else if ( colorEn <= (1-eps) )
            flag = Flag::mixed;
          else
          {
            flag = Flag::full;

            for ( const auto &intersection : intersections( gridView(), entity ) )
              if ( intersection.neighbor() && color[ index( intersection.outside() ) ] < eps )
              {
                flag = Flag::mixedfull;
                break;
              }
          }
        }

        for ( const auto &entity : elements( gridView() ) )
        {
          const auto idx = index( entity );

          if ( flags_[ idx ] == Flag::mixed || flags_[ idx ] == Flag::mixedfull )
            for ( const auto &intersection : intersections( gridView(), entity ) )
              if ( intersection.neighbor() )
              {
                Flag &flag = flags_[ index( intersection.outside() ) ];

                if ( flag == Flag::empty )
                  flag = Flag::activeempty;
                else if ( flag == Flag::full )
                  flag = Flag::activefull;
              }
        }
      }

    private:
      const GridView &gridView () const { return gridView_; }

      auto index ( const Entity &entity ) const -> decltype( std::declval< Mapper >().index( entity ) )
      {
        return mapper_.index( entity );
      }

      GridView gridView_;
      Mapper mapper_;
      std::vector< Flag > flags_;
    };


    template< class GV >
    inline static Flags< GV > flags ( const GV &gridView )
    {
      return Flags< GV > ( gridView );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGS_HH
