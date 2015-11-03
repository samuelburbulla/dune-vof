#ifndef DUNE_VOF_FLAGS_HH
#define DUNE_VOF_FLAGS_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune
{
  namespace VoF
  {

    // Flag the entities of a gridView as mixed or active

    template< class GridView, class Domain >
    class Flags 
    {

      public:

        const int dimworld = GridView::dimensionworld;
        typedef typename GridView::template Codim< 0 >::Entity Entity;
        //typedef typename Dune::MCMGElementLayout< 2 > MapperLayout;

        Flags ( const GridView &gridView, const Domain &domain ) 
         : _gridView ( gridView ), _domain ( domain ), _mapper ( gridView ), _mixed ( _mapper.size(), false ), _active ( _mapper.size(), false )
         {}

        bool isMixed ( const Entity& entity ) const
        {
          return _mixed[ _mapper.index( entity ) ];
        }

        bool isActive ( const Entity& entity ) const
        {
          return _active[ _mapper.index( entity ) ];
        }

        const std::size_t size() const
        {
          return _mapper.size();
        }

        int operator[] ( int i )
        {
          return _mixed[ i ] ? 2 : (_active[ i ] ? 1 : 0);
        }



        template< class ColorFunction >
        void reflag ( const ColorFunction& colorFunction, const double eps )
        {
 
          // mixed cells
          for( auto&& entity : elements( _gridView ) )
          {
            _active[  _mapper.index( entity ) ] = false;  

            if( colorFunction[ entity ] >= eps && colorFunction[ entity ] <= 1 - eps )
            {
              _mixed[ _mapper.index( entity ) ] = true;
            }
            else
            {
              _mixed[ _mapper.index( entity ) ] = false;
            }
          }
                  /*
        else if( colorFunction[ entity ] > 1 - eps )
          for( auto&& intersection : intersections( gridView, entity ) )
          {
            if ( intersection.neighbor() )
            {
              const Entity &neighbor = intersection.outside();
              
              if( colorFunction[ neighbor ] < eps )
              {
                flags.addMixedCell[ entity ];

                // the reconstruction is the edge
                const IntersectionGeometry isGeo = intersection.geometry();
                auto n = intersection.centerUnitOuterNormal();
                n *= -1.0;

                reconstruction[ entity ] = std::array< fvector, 3 >( {{ isGeo.corner( 0 ), isGeo.corner( 1 ), n }} );
              }
            }
          }
          */


          // active cells
          for( auto&& entity : elements( _gridView ))
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


      private:
      
        GridView _gridView;
        Domain _domain;
        Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > _mapper;
        std::vector< bool > _mixed;
        std::vector< bool > _active;

    };

  } // end of namespace VoF
} // end of namespace Dune

#endif
