#ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
#define DUNE_VOF_RECONSTRUCTIONSET_HH

#include <algorithm>
#include <vector>

#include <dune/common/deprecated.hh>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // ReconstructionSet
    // -----------------

    /**
     * \ingroup Other
     * \brief set of reconstructions
     *
     * \tparam  GV  grid view
     */
    template< class GV >
    struct ReconstructionSet
    {
      using GridView = GV;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

    public:
      using Reconstruction = HalfSpace< Coordinate >;

      using iterator = typename std::vector< Reconstruction >::iterator;
      using const_iterator = typename std::vector< Reconstruction >::const_iterator;
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );

      /**
       * \brief MPI communication handler
       */
      struct Exchange;

      explicit ReconstructionSet ( const GridView &gridView )
        : gridView_( gridView ),
          reconstructionSet_( indexSet().size( 0 ) )
       {}

      const Reconstruction& operator[] ( const Entity &entity ) const /*DUNE_DEPRECATED_MSG( "Use access via index instead." )*/ { return reconstructionSet_[ indexSet().index( entity ) ]; }
      Reconstruction& operator[] ( const Entity &entity ) /*DUNE_DEPRECATED_MSG( "Use access via index instead." )*/ { return reconstructionSet_[ indexSet().index( entity ) ]; }

      const Reconstruction& operator[] ( const Index &index ) const { return reconstructionSet_[ index ]; }
      Reconstruction& operator[] ( const Index &index ) { return reconstructionSet_[ index ]; }

      iterator begin () { return reconstructionSet_.begin(); }
      iterator end () { return reconstructionSet_.end(); }

      const_iterator begin () const { return reconstructionSet_.begin(); }
      const_iterator end () const { return reconstructionSet_.end(); }

      void clear() { std::fill( reconstructionSet_.begin(), reconstructionSet_.end(), Reconstruction() ); }

    private:
      const IndexSet &indexSet () const { return gridView_.indexSet(); }

      GridView gridView_;
      std::vector< Reconstruction > reconstructionSet_;
    };


    // Exchange class for MPI
    template< class GV >
    struct ReconstructionSet< GV >::Exchange : public Dune::CommDataHandleIF < Exchange, ReconstructionSet::Reconstruction >
    {
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
