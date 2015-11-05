#ifndef DUNE_VOF_COLORFUNCTION_HH
#define DUNE_VOF_COLORFUNCTION_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // ColorFunction for a given gridView

    template< class GridView >
    class ColorFunction 
    {

      public:

        typedef typename GridView::template Codim< 0 >::Entity Entity;

        ColorFunction ( const GridView &gridView ) 
         : _gridView( gridView ), _mapper ( gridView ), _color ( _mapper.size(), 0.0 )
         {}


        double& operator[] ( const Entity& entity )
        {
          assert ( _mapper.index( entity ) < _color.size() );
          return _color[ _mapper.index( entity ) ];
        }

        const double& operator[] ( const Entity& entity ) const
        {
          assert ( _mapper.index( entity ) < _color.size() );
          return _color[ _mapper.index( entity ) ];
        }

        double& operator[] ( const std::size_t i )
        {
          assert ( i < _color.size() );
          return _color[ i ];
        }

        const double& operator[] ( const std::size_t i ) const
        {
          assert ( i < _color.size() );
          return _color[ i ];
        }

        const double size() const
        {
          return _color.size();
        }

        const GridView& gridView() const
        {
          return _gridView;
        }

        void axpy ( const double a, ColorFunction &x )
        {
          assert( x.size() == _color.size() );

          for ( std::size_t i = 0; i < _color.size();  i++ )
          {
            double tmp = x[ i ];
            tmp *= a;
            _color[ i ] += tmp;
          }
        }

        void clear()
        {
          std::fill( _color.begin(), _color.end(), 0.0 );
        }


      private:

        GridView _gridView;
        Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > _mapper;
        std::vector< double > _color;
    };

  } // end of namespace VoF
} // end of namespace Dune

#endif
