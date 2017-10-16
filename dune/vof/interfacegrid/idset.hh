#ifndef DUNE_VOF_INTERFACEGRID_IDSET_HH
#define DUNE_VOF_INTERFACEGRID_IDSET_HH

#include <cassert>
#include <cstdint>

#include <type_traits>
#include <utility>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridIdSet
    // ------------------

    template< class Grid, class HostIdSet >
    class InterfaceGridIdSet
      : public IdSet< Grid, InterfaceGridIdSet< Grid, HostIdSet >, std::pair< typename HostIdSet::IdType, unsigned int > >
    {
      typedef InterfaceGridIdSet< Grid, HostIdSet > This;
      typedef IdSet< Grid, InterfaceGridIdSet< Grid, HostIdSet >, std::pair< typename HostIdSet::IdType, unsigned int > > Base;

      typedef typename std::remove_const_t< Grid >::Traits Traits;

    public:
      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;

      struct IdType
      {
        typedef typename HostIdSet::IdType HostId;

        IdType () : hostId_( 0 ), codimAndSubEntity_( 0 ) {};

        IdType ( HostId hostId, int codim, int subEntity )
          : hostId_( std::move( hostId ) ),
            codimAndSubEntity_( (static_cast< std::uint32_t >( codim ) << 24) | static_cast< std::uint32_t >( subEntity ) )
        {}

        bool operator== ( const IdType &other ) const { return (codimAndSubEntity_ == other.codimAndSubEntity_) && (hostId_ == other.hostId_); }
        bool operator!= ( const IdType &other ) const { return (codimAndSubEntity_ != other.codimAndSubEntity_) || (hostId_ != other.hostId_); }

        bool operator< ( const IdType &other ) const
        {
          return (codimAndSubEntity_ < other.codimAndSubEntity_) || ((codimAndSubEntity_ == other.codimAndSubEntity_) && (hostId_ < other.hostId_));
        }

        template< class C, class T >
        friend std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &out, const IdType &id )
        {
          return out << "<" << id.hostId_ << ", " << id.codimAndSubEntity_ << ">";
        }

      private:
        HostId hostId_;
        std::uint32_t codimAndSubEntity_;
      };

      InterfaceGridIdSet () = default;

      explicit InterfaceGridIdSet ( const HostIdSet &hostIdSet ) : hostIdSet_( &hostIdSet ) {}

      InterfaceGridIdSet ( const This &other ) : hostIdSet_( other.hostIdSet_ ) {}

      const This &operator= ( const This &other ) { hostIdSet_ = other.hostIdSet_; return *this; }

      explicit operator bool () const { return bool( hostIdSet_ ); }

      template< int cd >
      IdType id ( const typename Traits::template Codim< cd >::Entity &entity ) const
      {
        return id( entity, Dune::Codim< cd >() );
      }

      template< class Entity >
      IdType id ( const Entity &entity ) const
      {
        return id( entity, Dune::Codim< Entity::codimension >() );
      }

      template< int cd >
      IdType subId ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
      {
        return subId( entity, i, codim, Dune::Codim< cd >() );
      }

      //! subId method for all entities
      template< class Entity >
      IdType subId ( const Entity &entity, int i, unsigned int codim ) const
      {
        return subId( entity, i, codim, Dune::Codim< Entity::codimension >() );
      }

    protected:
      IdType id ( const typename Traits::template Codim< 0 >::Entity &entity, Dune::Codim< 0 > ) const
      {
        return IdType( hostIdSet().id( Grid::getRealImplementation( entity ).hostElement() ), 0, 0 );
      }

      template< int cd >
      IdType id ( const typename Traits::template Codim< cd >::Entity &entity, Dune::Codim< cd > ) const
      {
        return IdType( hostIdSet().id( Grid::getRealImplementation( entity ).hostElement() ), cd, Grid::getRealImplementation( entity ).subEntity() );
      }

      IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim, Dune::Codim< 0 > ) const
      {
        assert( codim <= static_cast< unsigned int >( dimension ) );
        return IdType( hostIdSet().id( Grid::getRealImplementation( entity ).hostElement() ), codim, i );
      }

      template< int cd >
      IdType subId ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim, Dune::Codim< cd > ) const
      {
        assert( (codim >= static_cast< unsigned int >( cd )) && (codim <= static_cast< unsigned int >( dimension )) );
        const auto &refElement = ReferenceElements< ctype, dimension >::general( Grid::getRealImplementation( entity ).hostElement().type() );
        return IdType( hostIdSet().id( Grid::getRealImplementation( entity ).hostElement() ), codim, refElement.subEntity( Grid::getRealImplementation( entity ).subEntity(), cd, i, codim ) );
      }

      const HostIdSet &hostIdSet () const { assert( *this ); return *hostIdSet_; }

      const HostIdSet *hostIdSet_ = nullptr;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_IDSET_HH
