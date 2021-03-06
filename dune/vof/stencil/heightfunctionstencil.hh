#ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCILS_HH
#define DUNE_VOF_HEIGHTFUNCTIONSTENCILS_HH

#include <dune/common/power.hh>
#include <dune/grid/spgrid/entity.hh>

namespace Dune
{
  namespace VoF
  {

    // HeightFunctionStencil
    // ---------------------

     /**
     * \ingroup Method
     * \brief  a single stencil for the height function method
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct HeightFunctionStencil
    {
      using GridView = GV;
      static constexpr int dim = GridView::dimension;

      using Entity = typename GridView::template Codim< 0 >::Entity;
      using EntityImpl = SPEntity< 0, dim, const typename GridView::Grid >;
      using EntityInfo = typename EntityImpl::EntityInfo;
      using MultiIndex = typename EntityInfo::MultiIndex;

      using Orientation = std::tuple< int, int >;

      static constexpr int noc = StaticPower< 3, dim-1 >::power;

    public:
      explicit HeightFunctionStencil ( const GridView &gridView, const EntityInfo &entityInfo, const Orientation &orientation, const double tup = 3 )
       : gridView_( gridView ), tup_( tup ), entityInfo_( entityInfo ), orientation_( orientation )
      {}

      std::size_t columns() const { return noc; }

      int tdown() const { return -tup_; }

      int tup() const { return tup_; }

      int effectiveTdown() const
      {
        for ( int i = -1; i >= tdown(); --i )
          if ( !valid( ( noc - 1 ) / 2, i ) )
            return  - 1 - i;
        return tup_;
      }

      Entity operator() ( const std::size_t c, const int t ) const
      {
        EntityInfo entityInfo( entityInfo_ );
        entityInfo.id() = getMultiIndex< dim >( c, t );
        entityInfo.update();
        return EntityImpl( entityInfo );
      }

      bool valid ( const std::size_t c, const int t ) const
      {
        return entityInfo_.gridLevel().template partition< All_Partition >().contains( getMultiIndex< dim >( c, t ), entityInfo_.partitionNumber() );
      }

      template < int dimension >
      auto getMultiIndex( const std::size_t c, const int t ) const
       -> typename std::enable_if< dimension == 2, MultiIndex >::type
      {
        MultiIndex m;

        int i = std::get< 0 >( orientation_ );
        int j = std::get< 1 >( orientation_ );

        m[ 1 - i ] += ( c - 1 ) * ( i == 0 ? -1 : 1 );
        m[ i ] += t;
        m *= j;

        m *= 2;
        return m + entityInfo_.id();
      }

      template < int dimension >
      auto getMultiIndex( const std::size_t c, const int t ) const
       -> typename std::enable_if< dimension == 3, MultiIndex >::type
      {
        MultiIndex m;

        int i = std::get< 0 >( orientation_ );
        int j = std::get< 1 >( orientation_ );

        m[ i ] += j * t;
        m[ (i+1)%3 ] += ( ( c % 3 ) - 1 ) * -j;
        m[ (i+2)%3 ] += ( ( c / 3 ) - 1 );

        m *= 2;
        return m + entityInfo_.id();
      }

    private:
      const GridView &gridView_;
      std::size_t tup_;
      const EntityInfo entityInfo_;
      const Orientation orientation_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCIL_HH
