#ifndef COLORFUNCTION_HH
#define COLORFUNCTION_HH

//- dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>

#include "../dataset.hh"


// ColorFunction
// -------------

template< class GV >
struct ColorFunction
  : public Dune::VoF::DataSet < GV, double >
{
  using GridView = GV;
  using ThisType = ColorFunction< GridView >;
  using BaseType = typename Dune::VoF::DataSet< GridView, double >;
  using ctype = double;

public:
  ColorFunction ( const GridView &gridView ) : BaseType( gridView ) {}

  void axpy ( const ctype a, ThisType &x )
  {
    assert( x.size() == this->size() );

    for ( const auto& entity : elements( this->gridView() ) )
      this->operator[]( entity ) += a * x[ entity ];
  }

  template< class BinaryInStream >
  void write ( BinaryInStream &in ) const
  {
    for ( const auto &entity : elements( this->gridView() ) )
      in.writeDouble( this->operator[]( entity ) );
  }

  template< class BinaryOutStream >
  void read ( BinaryOutStream &out )
  {
    for ( const auto &entity : elements( this->gridView() ) )
      out.readDouble( this->operator[]( entity ) );
  }
};

#endif // #ifndef COLORFUNCTION_HH
