#ifndef DUNE_VOF_FLAGS_HH
#define DUNE_VOF_FLAGS_HH

#include <algorithm>
#include <utility>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // Flag the entities of a gridView as mixed or active

    /**
     * \ingroup Method
     * \brief set of flags
     *
     * \tparam  GV  grid view
     */
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
        activefull  = 5,
        nan         = 6
      };

      using Mapper = MCMGMapper< GridView, MCMGElementLayout >;


    private:
      /**
       * \brief MPI communication handler
       *
       * \tparam  Reduce  reduction functor
       */
      template< class Reduce >
      struct Exchange;

    public:
      explicit Flags ( const GridView &gridView )
       : gridView_ ( gridView ), mapper_( gridView ), flags_( mapper_.size(), Flag::empty )
      {}

      const bool isMixed ( const Entity& entity ) const { return flags_[ index( entity ) ] == Flag::mixed; }
      const bool isFullAndMixed ( const Entity& entity ) const { return flags_[ index( entity ) ] == Flag::mixedfull; }

      const bool isActive ( const Entity& entity ) const
      {
        const Flag &flag = flags_[ index( entity ) ];
        return ( flag == Flag::activeempty ) || ( flag == Flag::activefull );
      }

      const Flag& operator[] ( const Entity& entity ) const { return flags_[ index( entity ) ]; }
      Flag& operator[] ( const Entity& entity ) { return flags_[ index( entity ) ]; }

      const std::size_t size() const { return mapper_.size(); }

      /**
       * \brief update set of flags
       *
       * \tparam  DF  discrete function type
       * \param color discrete function
       * \param eps   marker tolerance
       */
      template< class DF >
      double reflag ( const DF& color, const double eps )
      {
        double elapsedTime = - MPI_Wtime();

        for ( const auto &entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          const auto idx = index( entity );
          Flag &flag = flags_[ idx ];
          const auto colorEn = color[ entity ];

          if ( colorEn < eps )
            flag = Flag::empty;
          else if ( colorEn <= (1-eps) )
            flag = Flag::mixed;
          else
          {
            if ( colorEn != colorEn )
              flag = Flag::nan;
            else
            {
              flag = Flag::full;

              for ( const auto &intersection : intersections( gridView(), entity ) )
                if ( intersection.neighbor() && color[ intersection.outside() ] < eps )
                {
                  flag = Flag::mixedfull;
                  break;
                }
            }
          }
        }
        elapsedTime += MPI_Wtime();

        auto exchange1 = makeExchange( []( Flag a, Flag b ){ return b; } );
        color.gridView().grid().communicate( exchange1, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );

        elapsedTime -= MPI_Wtime();
        for ( const auto &entity : elements( gridView(), Partitions::interiorBorder ) )
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
        elapsedTime += MPI_Wtime();

        auto exchange2 = makeExchange( [] ( Flag a, Flag b  ) { return std::max( a, b ); } );
        color.gridView().grid().communicate( exchange2, Dune::All_All_Interface, Dune::ForwardCommunication );

        return elapsedTime;
      }

    private:
      template< class Reduce >
      Exchange< Reduce > makeExchange ( Reduce reduce ) { return Exchange< Reduce >( *this, std::move( reduce ) ); }

      const GridView &gridView () const { return gridView_; }

      auto index ( const Entity &entity ) const -> decltype( std::declval< Mapper >().index( entity ) )
      {
        return mapper_.index( entity );
      }

      GridView gridView_;
      Mapper mapper_;
      std::vector< Flag > flags_;
    };


    /**
     * \ingroup Method
     * \brief generate set of flags
     *
     * \tparam GV grid view
     * \param gridView
     */
    template< class GV >
    inline static Flags< GV > flags ( const GV &gridView )
    {
      return Flags< GV > ( gridView );
    }


    // Exchange class for MPI
    template< class GV >
    template < class Reduce >
    class Flags< GV >::Exchange : public Dune::CommDataHandleIF < Exchange < Reduce >, Flag >
    {
      public:

        Exchange ( Flags &flags, Reduce reduce ) : flags_ ( flags ), reduce_( std::move( reduce ) ) {}

        typedef typename Flags::Flag FlagType;

        const bool contains ( const int dim, const int codim ) const { return ( codim == 0 ); }

        const bool fixedsize ( const int dim, const int codim ) const { return true; }

        template < class Entity >
        const size_t size ( const Entity &e ) const { return 1; }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {
          buff.write( flags_[ e ] );
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {}

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
        {
          assert( n == 1 );
          FlagType x ;
          buff.read( x );
          FlagType &y = flags_[ e ];
          y = reduce_( y, x );
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
        {}

      private:
        Flags &flags_;
        Reduce reduce_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGS_HH
