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

      static constexpr int noc = std::pow( 3, dim-1 );

    public:
      explicit HeightFunctionStencil ( const GridView &gridView )
       : gridView_( gridView ), tup_( 3 ), columns_( noc )
      {
        for ( int i = 0; i < noc; ++i )
          columns_[ i ].resize( 2 * tup_ + 1 );
      }

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

      Entity operator() ( std::size_t i, int t ) const
      {
        assert( valid( i, t ) );
        const EntityInfo entityInfo = columns_[ i ][ t + tup_ ];
        return EntityImpl( entityInfo );
      }

      void add ( const MultiIndex &m, const Entity &entity )
      {
        int columnId = m[ 1 ] + 1;

        if ( dim == 3 )
          columnId += 3 * ( m[ 2 ] + 1 );

        columns_[ columnId ][ m[ 0 ] + tup_ ] = GridView::Grid::getRealImplementation( entity ).entityInfo();
      }

      bool valid ( std::size_t i, int t ) const
      {
        auto multiIndex = columns_[ i ][ t + tup_ ].id();
        return multiIndex != MultiIndex::zero();
      }

    private:
      const GridView &gridView_;
      std::size_t tup_;
      std::vector< std::vector< EntityInfo > > columns_;
    };

    // HeightFunctionStencils
    // ----------------------

     /**
     * \ingroup Method
     * \brief  set of height function stencils
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct HeightFunctionStencils
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using MultiIndex = typename SPEntity< 0, GV::dimension, typename GV::Grid >::MultiIndex;
      using Stencil = HeightFunctionStencil< GridView >;
      static constexpr int dim = GridView::dimension;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit HeightFunctionStencils ( const GridView& gridView )
       : gridView_( gridView )
      {
        for ( std::size_t i = 0; i < 2 * dim; ++i )
          for ( std::size_t id = 0; id < indexSet().size( 0 ); ++id )
            stencils_[ i ].emplace_back( gridView_ );

        for ( const auto& entity : elements( gridView ) )
        {
          deltaX_ = std::pow( entity.geometry().volume(), 1.0 / dim );
          break;
        }

        initialize();
      }

      const Stencil& operator() ( const Coordinate &normal, const Entity& entity ) const
      {
        std::size_t dir = 0;
        double max = std::numeric_limits< double >::min();
        for ( std::size_t i = 0; i < dim; ++i )
         if ( std::abs( normal[ i ] ) > max )
         {
          dir = i;
          max = std::abs( normal[ i ] );
         }

        dir += dim * ( ( normal[ dir ] > 0 ) ? 0 : 1 );

        return stencils_[ dir ][ indexSet().index( entity ) ];
      }

      double deltaX() const { return deltaX_; }

    private:
      const GridView& gridView() const { return gridView_; }
      const IndexSet& indexSet() const { return gridView().indexSet(); }

      void initialize()
      {
        for ( const Entity& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          Index id = indexSet().index( entity );
          auto multiIndex = GridView::Grid::getRealImplementation( entity ).entityInfo().id();

          for ( std::size_t i = 0; i < dim; ++i )
            for ( int j = -1; j <= 1; j = j + 2 )
            {
              MultiIndex direction = MultiIndex::zero();
              direction[ i ] = j;

              Stencil &stencil = stencils_[ i + 0.5 * ( j + 1 ) * dim ][ id ];

              for ( const auto& other : elements( gridView(), Partitions::all ) )
              {
                auto multiIndexOther = GridView::Grid::getRealImplementation( other ).entityInfo().id();

                MultiIndex m = multiIndexOther - multiIndex;
                m /= 2;
                m *= j;

                if ( dim == 2 )
                {
                  if ( i == 0 ) m[1] *= -1;
                }

                if ( dim == 3 )
                {
                  if ( i == 0 && j == -1 ) m[2] *= -1;
                  if ( i == 0 && j == 1 ) m[1] *= -1;
                  if ( i == 1 && j == -1 ) m[2] *= -1;
                  if ( i == 2 && j == -1 ) m[1] *= -1;
                  if ( i == 2 && j == 1 ) m[1] *= -1;
                }

                if ( m[ i ] < stencil.tdown() || m[ i ] > stencil.tup() )
                  continue;

                bool valid = true;
                for ( std::size_t k = 0; k < dim; ++k )
                {
                  if ( k == i )
                    continue;

                  valid &= ( std::abs( m[ k ] ) <= 1 );
                }
                if ( !valid )
                  continue;

                if ( i == 1 )
                  std::swap( m[0], m[1] );

                if ( i == 2 )
                {
                  std::swap( m[1], m[2] );
                  std::swap( m[0], m[1] );
                }

                stencil.add( m, other );
              }
            }

        }
      }

      GridView gridView_;
      std::array< std::vector< Stencil >, 2 * dim > stencils_;
      double deltaX_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCIL_HH
