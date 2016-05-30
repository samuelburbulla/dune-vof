#ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
#define DUNE_VOF_EDGENEIGHBORSSTENCIL_HH

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

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
      using Mapper = Dune::VoF::MCMGMapper< GridView, Dune::MCMGElementLayout >;

    public:
      explicit EdgeNeighborsStencil ( const GridView& gridView )
       : gridView_( gridView ), mapper_( gridView_ ), stencils_( mapper().size() )
      {
        initialize();
      }

      const Stencil& operator[] ( const Entity& entity ) const
      {
        return stencils_[ mapper().index( entity ) ];
      }

    private:
      const GridView& gridView () const { return gridView_; }
      const Mapper& mapper() const { return mapper_; }

      void initialize()
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          Stencil& stencil = stencils_[ mapper().index( entity ) ];
          for ( const auto& intersection : intersections( gridView(), entity ) )
            if ( intersection.neighbor() )
            {
              const auto &neighbor = intersection.outside();
              stencil.push_back( neighbor );
            }
        }
      }

      GridView gridView_;
      Mapper mapper_;
      std::vector< Stencil > stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_EDGENEIGHBORSSTENCIL_HH
