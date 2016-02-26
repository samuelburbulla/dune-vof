#ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
#define DUNE_VOF_RECONSTRUCTIONSET_HH

#include <vector>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>

//- dune-vof includes
#include <dune/vof/geometry/halfspace.hh>

namespace Dune
{
  namespace VoF
  {

    // ReconstructionSet
    // -----------------

    template< class GV >
    struct ReconstructionSet
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;

    private:
      using Mapper = Dune::VoF::MCMGMapper< GridView, Dune::MCMGElementLayout >;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      using Reconstruction = HalfSpace< Coordinate >;
      using Intersections = std::vector< Coordinate >;

      using iterator = typename std::vector< Reconstruction >::iterator;
      using const_iterator = typename std::vector< Reconstruction >::const_iterator;

      explicit ReconstructionSet ( const GridView &gridView )
       : mapper_( gridView ), reconstructionSet_( mapper().size() ), intersectionsSet_( mapper().size() )
       {}

      const Reconstruction& operator[] ( const Entity &entity ) const { return reconstructionSet_[ mapper().index( entity ) ]; }
      Reconstruction& operator[] ( const Entity &entity ) { return reconstructionSet_[ mapper().index( entity ) ]; }

      iterator begin () { return reconstructionSet_.begin(); }
      const_iterator begin () const { return reconstructionSet_.begin(); }

      iterator end () { return reconstructionSet_.end(); }
      const_iterator end () const { return reconstructionSet_.end(); }

      const Intersections& intersections ( const Entity &entity ) const { return intersectionsSet_[ mapper().index( entity ) ]; }
      Intersections& intersections ( const Entity &entity ) { return intersectionsSet_[ mapper().index( entity ) ]; }

      const std::vector< Intersections >& intersectionsSet () const { return intersectionsSet_; }

      const void clear()
      {
        std::fill( reconstructionSet_.begin(), reconstructionSet_.end(), Reconstruction() );
        std::fill( intersectionsSet_.begin(), intersectionsSet_.end(), Intersections() );
      }

      struct Exchange;

    private:
      const Mapper &mapper () const { return mapper_; }

      Mapper mapper_;
      std::vector< Reconstruction > reconstructionSet_;
      std::vector< Intersections > intersectionsSet_;
    };


    // Exchange class for MPI
    template< class GV >
    class ReconstructionSet< GV >::Exchange : public Dune::CommDataHandleIF < Exchange, ReconstructionSet::Reconstruction >
    {
      public:
        Exchange ( ReconstructionSet &reconstructionSet ) : reconstructionSet_ ( reconstructionSet ) {}

        typedef typename ReconstructionSet::Reconstruction ReconstructionType;

        const bool contains ( const int dim, const int codim ) const { return ( codim == 0 ); }

        const bool fixedsize ( const int dim, const int codim ) const { return true; }

        template < class Entity >
        const size_t size ( const Entity &e ) const { return 1; }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {
          buff.write( reconstructionSet_[ e ] );
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {}

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, size_t n )
        {
          ReconstructionType x ;
          buff.read( x );
          reconstructionSet_[ e ] = x;
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
        {}

      private:
        ReconstructionSet &reconstructionSet_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
