#ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
#define DUNE_VOF_EDGENEIGHBORSSTENCIL_HH

namespace Dune
{
  namespace VoF
  {

    // EdgeNeighborsStencil
    // ----------------------

    /**
     * \ingroup Method
     * \brief  set of edge neighbors stencils
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct EdgeNeighborsStencil
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Stencil = std::vector< Entity >;
    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );
    public:
      explicit EdgeNeighborsStencil ( const GridView& gridView )
       : gridView_( gridView ), stencils_( indexSet().size( 0 ) )
      {
        initialize();
      }

      const Stencil& operator[] ( const Entity& entity ) const
      {
        return stencils_[ indexSet().index( entity ) ];
      }

    private:
      const GridView& gridView () const { return gridView_; }
      const IndexSet& indexSet() const { return gridView().indexSet(); }

      void initialize()
      {
        for ( const auto& entity : elements( gridView(), Partitions::all ) )
        {
          Stencil& stencil = stencils_[ indexSet().index( entity ) ];
          for ( const auto& intersection : intersections( gridView(), entity ) )
            if ( intersection.neighbor() )
              stencil.push_back( intersection.outside() );

        }
      }

      GridView gridView_;
      std::vector< Stencil > stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
