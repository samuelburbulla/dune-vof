#ifndef DUNE_VOF_DATASET__HH
#define DUNE_VOF_DATASET__HH

#include <algorithm>
#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

namespace Dune
{
  namespace VoF
  {

    // DataSet
    // -------

    /**
     * \ingroup Other
     * \brief storage of grid data
     *
     * \tparam  GV  grid view
     * \tparam  T   data type
     */
    template< class GV, class T >
    struct DataSet
    {
      using GridView = GV;
      using DataType = T;
      using Entity = typename GridView::template Codim< 0 >::Entity;

    private:
      using IndexSet = typename GridView::IndexSet;

    public:
      using iterator = typename std::vector< DataType >::iterator;
      using const_iterator = typename std::vector< DataType >::const_iterator;
      using Index = typename IndexSet::IndexType;

      /**
       * \brief MPI communication handler
       */
      template < class Reduce >
      struct Exchange;

      explicit DataSet ( GridView gridView )
        : gridView_( gridView ),
          dataSet_( indexSet().size( 0 ) )
       {}

      const DataType& operator[] ( const Entity &entity ) const { return dataSet_[ indexSet().index( entity ) ]; }
      DataType& operator[] ( const Entity &entity ) { return dataSet_[ indexSet().index( entity ) ]; }

      const DataType& operator[] ( const Index &index ) const { return dataSet_[ index ]; }
      DataType& operator[] ( const Index &index ) { return dataSet_[ index ]; }

      iterator begin () { return dataSet_.begin(); }
      iterator end () { return dataSet_.end(); }

      const_iterator begin () const { return dataSet_.begin(); }
      const_iterator end () const { return dataSet_.end(); }

      void clear() { std::fill( dataSet_.begin(), dataSet_.end(), DataType() ); }

      std::size_t size() const { return dataSet_.size(); }

      const GridView &gridView () const { return gridView_; }

      template< class Reduce >
      void communicate ( Dune::InterfaceType interface, Reduce reduce )
      {
        auto exchange = Exchange< Reduce > ( *this, std::move( reduce ) );
        gridView_.communicate( exchange, interface, Dune::ForwardCommunication );
      }

      void communicate ()
      {
        auto reduce = [] ( DataType a, DataType b ) { return a; };
        auto exchange = Exchange< decltype( reduce ) > ( *this, std::move( reduce ) );
        gridView_.communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );
      }

    private:
      const IndexSet &indexSet () const { return gridView_.indexSet(); }
      GridView gridView_;
      std::vector< DataType > dataSet_;
    };



    // Exchange class for MPI
    template< class GV, class T >
    template < class Reduce >
    struct DataSet< GV, T >::Exchange
     : public Dune::CommDataHandleIF < Exchange< Reduce >, T >
    {
        Exchange ( DataSet &dataSet, Reduce reduce ) : dataSet_ ( dataSet ), reduce_( std::move( reduce ) ) {}

        const bool contains ( const int dim, const int codim ) const { return ( codim == 0 ); }

        const bool fixedsize ( const int dim, const int codim ) const { return true; }

        template < class Entity >
        const size_t size ( const Entity &e ) const { return 1; }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {
          buff.write( dataSet_[ e ] );
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void gather ( MessageBuffer &buff, const Entity &e ) const
        {}

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
        {
          T x ;
          buff.read( x );
          T &y = dataSet_[ e ];
          y = reduce_( x, y );
        }

        template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
        void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
        {}

      private:
        DataSet &dataSet_;
        Reduce reduce_;
    };


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_DATASET__HH
