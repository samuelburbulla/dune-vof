#ifndef DUNE_VOF_FLAGS_HH
#define DUNE_VOF_FLAGS_HH

//- dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // Flag the entities of a gridView as mixed or active

    template< class GridView, class Domain >
    struct Flags
    {
      const int dimworld = GridView::dimensionworld;
      typedef typename GridView::template Codim< 0 >::Entity Entity;
      //typedef typename Dune::MCMGElementLayout< 2 > MapperLayout;

      Flags ( const GridView &gridView, const Domain &domain )
       : _gridView ( gridView ), _domain ( domain ), _mapper ( gridView ), _mixed ( _mapper.size(), false ),
         _fullandmixed ( _mapper.size(), false ), _active ( _mapper.size(), false )
      {}

      const bool isMixed ( const Entity& entity ) const { return _mixed[ _mapper.index( entity ) ]; }
      const bool isFullAndMixed ( const Entity& entity ) const { return _fullandmixed[ _mapper.index( entity ) ]; }
      const bool isActive ( const Entity& entity ) const { return _active[ _mapper.index( entity ) ]; }

      const std::size_t size() const { return _mapper.size(); }

      const int operator[] ( const int i ) const
      {
        if ( _fullandmixed[ i ] ) return 3;
        if ( _mixed[ i ] ) return 2;
        if ( _active[ i ] ) return 1;
        return 0;
      }

      template< class ColorFunction >
      void reflag ( const ColorFunction& colorFunction, const double eps )
      {
        // mixed cells
        for( const auto &entity : elements( _gridView ) )
        {
          _active[ _mapper.index( entity ) ] = false;
          _fullandmixed[ _mapper.index( entity ) ] = false;
          _mixed[ _mapper.index( entity ) ] = ( colorFunction[ entity ] >= eps && colorFunction[ entity ] <= 1 - eps );


          if( colorFunction[ entity ] > 1 - eps )
            for( auto&& intersection : intersections( _gridView, entity ) )
            {
              if ( intersection.neighbor() )
              {
                const Entity &neighbor = intersection.outside();

                if ( colorFunction[ neighbor ] < eps )
                  _fullandmixed[ _mapper.index( entity ) ] = true;
              }
            }
        }

        // active cells
        for( const auto &entity : elements( _gridView ))
        {
          if ( _mixed[ _mapper.index( entity ) ] )
          {
            for( auto&& is : intersections( _gridView,  entity ) )
            {
              if( is.neighbor() )
              {
                auto neighbor = is.outside();
                if( !_mixed[ _mapper.index( neighbor ) ] )
                  _active[  _mapper.index( neighbor ) ] = true;
              }
            }
          }
        }
      }

    private:
      GridView _gridView;
      Domain _domain;
      MultipleCodimMultipleGeomTypeMapper< GridView, MCMGElementLayout > _mapper;
      std::vector< bool > _mixed;
      std::vector< bool > _fullandmixed;
      std::vector< bool > _active;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_FLAGS_HH
