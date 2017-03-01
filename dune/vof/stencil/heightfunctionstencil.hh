#ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCILS_HH
#define DUNE_VOF_HEIGHTFUNCTIONSTENCILS_HH

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

      static constexpr int noc = std::pow( 3, dim-1 );

    public:
      explicit HeightFunctionStencil ( const GridView &gridView, const EntityInfo &entityInfo, const Orientation &orientation )
       : gridView_( gridView ), tup_( 3 ), entityInfo_( entityInfo ), orientation_( orientation )
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
        entityInfo.id() = getMultiIndex( c, t );
        entityInfo.update();
        return EntityImpl( entityInfo );
      }

      bool valid ( const std::size_t c, const int t ) const
      {
        return entityInfo_.gridLevel().template partition< All_Partition >().contains( getMultiIndex( c, t ), entityInfo_.partitionNumber() );
      }

      MultiIndex getMultiIndex( const std::size_t c, const int t ) const
      {
        MultiIndex m;

      #if GRIDDIM == 2
        m[ 0 ] += c - 1;
        m[ 1 ] += t;

      #elif GRIDDIM == 3
        m[ 0 ] += ( c % 3 ) - 1;
        m[ 1 ] += ( t / 3 ) - 1;
        m[ 2 ] += t;

      #endif

        int i = std::get< 0 >( orientation_ );
        int j = std::get< 1 >( orientation_ );

        m *= 2;
        m *= j;

      #if GRIDDIM == 2
        if ( i == 0 )
        {
          std::swap( m[0], m[1] );
          m[1] *= -1;
        }

      #elif GRIDDIM == 3
        if ( i == 0 && j == -1 ) m[2] *= -1;
        if ( i == 0 && j == 1 ) m[1] *= -1;
        if ( i == 1 && j == -1 ) m[2] *= -1;
        if ( i == 2 && j == -1 ) m[1] *= -1;
        if ( i == 2 && j == 1 ) m[1] *= -1;

      #endif

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
