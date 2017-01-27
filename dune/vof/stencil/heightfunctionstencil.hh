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
     * \tparam  Entity  entity
     */
    template< class E, class M >
    struct HeightFunctionStencil
    {
      using Entity = E;
      using MultiIndex = M;
      static constexpr int noc = 3;

    public:
      explicit HeightFunctionStencil ()
       : tdown_( 3 ), tup_( 3 ), offset_( tdown_ ), maxTdown_( 0 ), maxTup_( 0 ), nocL_( -1 ), nocR_( 1 ), nocOffset_( -nocL_ ), minC_( 0 ), maxC_( 0 ), columns_( columns() )
      {
        for ( int i = cMin(); i < cMax(); ++i )
          columns_[ nocOffset_ + i ].resize( tdown() + tup() + 1 );
      }

      void add ( const MultiIndex &m, const Entity &entity )
      {
        columns_[ m[ 1 ] + 1 ][ m[ 0 ] + offset_ ] = entity;

        if ( m[ 0 ] < -maxTdown_ )
          maxTdown_ = -m[ 0 ];
        else if ( m[ 0 ] > maxTup_ )
          maxTup_ = m[ 0 ];

        if ( m[ 1 ] < minC_ )
          minC_ = m[ 1 ];
        if ( m[ 1 ] > maxC_ )
          maxC_ = m[ 1 ];
      }

      void fixBounds ()
      {
        tdown_ = maxTdown_;
        tup_ = maxTup_;
        nocL_ = minC_;
        nocR_ = maxC_;
      }

      int cMin() const { return nocL_; }
      int cMax() const { return columns() + nocL_; }
      std::size_t columns() const { return 1 - nocL_ + nocR_; }
      int tdown() const { return tdown_; }
      int tup() const { return tup_; }
      const Entity &operator() ( std::size_t i, int t ) const { return columns_[ nocOffset_ + i ][ t + offset_ ]; }

    private:
      std::size_t tdown_, tup_, offset_;
      int maxTdown_, maxTup_, nocL_, nocR_, nocOffset_, minC_, maxC_;
      std::vector< std::vector< Entity > > columns_;
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
      using Stencil = HeightFunctionStencil< Entity, MultiIndex >;
      static constexpr int dim = GridView::dimension;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit HeightFunctionStencils ( const GridView& gridView )
       : gridView_( gridView )
      {
        std::size_t size = indexSet().size( 0 );
        for ( std::size_t i = 0; i < 2 * dim; ++i )
          stencils_[ i ].resize( size );

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
                m[ i ] *= j;

                if ( m[ i ] < -stencil.tdown() || m[ i ] > stencil.tup() )
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

                std::swap( m[0], m[i] );
                stencil.add( m, other );
              }

              stencil.fixBounds();
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
