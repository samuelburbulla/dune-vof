#ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
#define DUNE_VOF_RECONSTRUCTIONSET_HH

//-dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>

#include "hyperplane.hh"

namespace Dune
{
  namespace VoF
  {


    // ReconstructionSet
    // -----------------


    template< class GridView >
    class ReconstructionSet
    {

      public:

        static constexpr int dimworld = GridView::dimensionworld;
        typedef typename GridView::ctype ctype;
        typedef typename Dune::FieldVector< ctype, dimworld > fvector;
        typedef typename GridView::template Codim< 0 >::Entity Entity;
        typedef typename std::vector< Hyperplane< fvector > >::iterator iterator;
        typedef typename std::vector< Hyperplane< fvector > >::const_iterator const_iterator;

        ReconstructionSet ( const GridView &gridView )
         : _mapper ( gridView ), _reconstructionSet ( _mapper.size() ), _reconstructionIntersections ( _mapper.size() )
         {}

        const Dune::VoF::Hyperplane< fvector >& operator[] ( const Entity &entity ) const
        {
          assert ( _mapper.index( entity ) < _reconstructionSet.size() );
          return _reconstructionSet[ _mapper.index( entity ) ];
        }

        Dune::VoF::Hyperplane< fvector >& operator[] ( const Entity &entity ) { return _reconstructionSet[ _mapper.index( entity ) ]; }

        iterator begin () { return _reconstructionSet.begin(); }
        const_iterator begin () const { return _reconstructionSet.begin(); }

        iterator end () { return _reconstructionSet.end(); }
        const_iterator end () const { return _reconstructionSet.end(); }

        const std::vector< fvector >& intersections ( const Entity &entity ) const
        {
          assert ( _mapper.index( entity ) < _reconstructionIntersections.size() );
          return _reconstructionIntersections[ _mapper.index( entity ) ];
        }

        std::vector< fvector >& intersections ( const Entity &entity )
        {
          assert ( _mapper.index( entity ) < _reconstructionIntersections.size() );
          return _reconstructionIntersections[ _mapper.index( entity ) ];
        }

        const std::vector< std::vector< fvector > >& intersections ( ) const
        {
          return _reconstructionIntersections;
        }


        /* right class for definition?
        const bool isInner ( const Entity &entity, const fvector &x ) const
        {
          auto h = _reconstructionSet[ _mapper.index( entity ) ];
          return h.n * x + h.p >= 0;
        }
        */

        const void clear()
        {
          Dune::VoF::Hyperplane< fvector > h;
          std::fill( _reconstructionSet.begin(), _reconstructionSet.end(), h );

          std::vector< fvector > v;
          std::fill( _reconstructionIntersections.begin(), _reconstructionIntersections.end(), v );
        }


      private:

        Dune::MultipleCodimMultipleGeomTypeMapper< GridView, Dune::MCMGElementLayout > _mapper;
        std::vector< Dune::VoF::Hyperplane< fvector > > _reconstructionSet; // should be private with range-based iterator
        std::vector< std::vector< fvector > > _reconstructionIntersections; // should be private with range-based iterator
    };



  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTIONSET_HH
