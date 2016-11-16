#ifndef COLORFUNCTION_HH
#define COLORFUNCTION_HH

//- dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>


// ColorFunction
// -------------

template< class GV >
struct ColorFunction
{
  typedef GV GridView;
  typedef typename GridView::template Codim< 0 >::Entity Entity;
  typedef double ctype;
  struct Exchange;

private:
  typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > Mapper;

public:
  ColorFunction ( const GridView &gridView )
   : gridView_( gridView ), mapper_( gridView ), color_( mapper_.size(), 0.0 )
   {}

  ctype& operator[] ( const Entity& entity ) { return color_[ mapper_.index( entity ) ]; }
  const ctype& operator[] ( const Entity& entity ) const { return color_[ mapper_.index( entity ) ]; }

  ctype& operator[] ( const int i ) { return color_[ i ]; }
  const ctype& operator[] ( const int i ) const {return color_[ i ]; }

  std::size_t size () const { return color_.size(); }

  const GridView &gridView () const { return gridView_; }

  void axpy ( const ctype a, ColorFunction &x )
  {
    assert( x.size() == color_.size() );

    for ( std::size_t i = 0; i < color_.size();  ++i )
      color_[ i ] += a * x[ i ];
  }

  void clear () { std::fill( color_.begin(), color_.end(), 0.0 ); }

  void communicate ()
  {
    auto exchange = Exchange ( *this );
    gridView_.communicate( exchange, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );
  }

private:
  GridView gridView_;
  Mapper mapper_;
  std::vector< ctype > color_;
};

// Exchange class for MPI
template< class GV >
struct ColorFunction< GV >::Exchange : public Dune::CommDataHandleIF < Exchange, ColorFunction::ctype >
{
    Exchange ( ColorFunction &color ) : color_ ( color ) {}

    typedef typename ColorFunction::ctype ctype;

    const bool contains ( const int dim, const int codim ) const { return ( codim == 0 ); }

    const bool fixedsize ( const int dim, const int codim ) const { return true; }

    template < class Entity >
    const size_t size ( const Entity &e ) const { return 1; }

    template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
    void gather ( MessageBuffer &buff, const Entity &e ) const
    {
      buff.write( color_[ e ] );
    }

    template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
    void gather ( MessageBuffer &buff, const Entity &e ) const
    {}

    template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension == 0, int >::type = 0 >
    void scatter ( MessageBuffer &buff, const Entity &e, size_t n )
    {
      ctype x ;
      buff.read( x );
      color_[ e ] = x;
    }

    template < class MessageBuffer, class Entity, typename std::enable_if< Entity::codimension != 0, int >::type = 0 >
    void scatter ( MessageBuffer &buff, const Entity &e, std::size_t n )
    {}

  private:
    ColorFunction &color_;
};

#endif // #ifndef COLORFUNCTION_HH
