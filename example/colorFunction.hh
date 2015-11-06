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

private:
  typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > Mapper;

public:
  ColorFunction ( const GridView &gridView )
   : gridView_( gridView ), mapper_( gridView ), color_( mapper_.size(), 0.0 )
   {}

  double& operator[] ( const Entity& entity ) { return color_[ mapper_.index( entity ) ]; }
  const double& operator[] ( const Entity& entity ) const { return color_[ mapper_.index( entity ) ]; }

  double& operator[] ( const int i ) { return color_[ i ]; }
  const double& operator[] ( const int i ) const {return color_[ i ]; }

  const double size () const { return color_.size(); }

  const GridView &gridView () const { return gridView_; }

  void axpy ( const double a, ColorFunction &x )
  {
    assert( x.size() == color_.size() );

    for ( std::size_t i = 0; i < color_.size();  ++i )
      color_[ i ] += a * x[ i ];
  }

  void clear () { std::fill( color_.begin(), color_.end(), 0.0 ); }

private:
  GridView gridView_;
  Mapper mapper_;
  std::vector< double > color_;
};

#endif // #ifndef COLORFUNCTION_HH
