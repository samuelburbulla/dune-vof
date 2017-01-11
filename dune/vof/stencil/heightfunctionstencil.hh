#ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCIL_HH
#define DUNE_VOF_HEIGHTFUNCTIONSTENCIL_HH

namespace Dune
{
  namespace VoF
  {

    // HeightFunctionStencil
    // ---------------------

     /**
     * \ingroup Method
     * \brief  set of height function stencils
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct HeightFunctionStencil
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Stencil = std::vector< Entity >;
      static constexpr int dim = GridView::dimension;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      explicit HeightFunctionStencil ( const GridView& gridView )
       : gridView_( gridView )
      {
        std::size_t size = indexSet().size( 0 );
        for ( std::size_t i = 0; i < dim; ++i )
          stencils_[ i ].resize( size );
        initialize();
      }

      const Stencil& operator() ( std::size_t direction, const Entity& entity ) const
      {
        return stencils_[ direction ][ indexSet().index( entity ) ];
      }

    private:
      const GridView& gridView() const { return gridView_; }
      const IndexSet& indexSet() const { return gridView().indexSet(); }

      void initialize()
      {
        for ( const Entity& entity : elements( gridView(), Partitions::interiorBorder ) )
        {
          Index id = indexSet().index( entity );
          Coordinate center = entity.geometry().center();
          const double deltaX = std::sqrt( entity.geometry().volume() );

          for ( const Entity& other : elements( gridView(), Partitions::all ) )
          {
            if ( other == entity )
              continue;

            Coordinate d = other.geometry().center() - center;

            for ( std::size_t i = 0; i < dim; ++i )
            {
              Coordinate direction ( 0.0 );
              direction[ i ] = 1.0;

              double dDir = d * direction;
              Coordinate dOrth = d;
              dOrth.axpy( - dDir, direction );

              if ( std::abs( dDir ) <= 3.0 * deltaX && dOrth.two_norm() <= std::sqrt( dim - 1.0 ) * deltaX )
                stencils_[ i ][ id ].push_back( other );
            }
          }
        }
      }

      GridView gridView_;
      std::array< std::vector< Stencil >, dim > stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_HEIGHTFUNCTIONSTENCIL_HH
