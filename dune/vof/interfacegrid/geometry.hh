#ifndef DUNE_VOF_INTERFACEGRID_GEOMETRY_HH
#define DUNE_VOF_INTERFACEGRID_GEOMETRY_HH

#include <cassert>
#include <cstddef>

#include <type_traits>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  namespace VoF
  {

    // PolygonGeometry
    // ---------------

    template< class ct, int cdim >
    struct PolygonGeometry;

    template< class ct, int cdim >
    struct PolygonGeometry< ct, 3 >
    {
      typedef ct ctype;
      
      static const int mydimension = 2;
      static const int coorddimension = 3;
      
      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
      
      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

      InterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize )
        : normal_( normal ), cbegin_( cbegin ), csize_( csize )
      {
        assert( csize >= 3u );
      }

      int corners () const { return csize_; }

      const GlobalCoordinate &corner ( int i ) const { assert( (i >= 0) && (i < corners) ); return cbegin_[ i ]; }

      GlobalCoordinate center () const
      {
        GlobalCoordinate center( 0 );
        ctype volume = 0;
        for( std::size_t i = 0; i < csize_; ++i )
        {
          const GlobalCoordinate &x = cbegin_[ i ];
          const GlobalCoordinate &y = cbegin_[ (i+1) % csize_ ];
          const ctype weight = det( x, y, normal );
          center.axpy( weight, x + y );
          volume += weight;
        }
        return center *= ctype( 1 ) / (ctype( 3 ) * volume);
      }

      ctype volume () const
      {
        using std::abs;

        ctype volume = 0;
        for( std::size_t i = 0; i < csize_; ++i )
          volume += det( cbegin_[ i ], cbegin_[ (i+1) % csize_ ], normal );
        return abs( volume / ctype( 2 ) );
      }

      GeometryType type () const noexcept { return GeometryType( GeometryType::none, mydimension ); }

      bool affine () const
      {
        DUNE_THROW( InvalidStateException, "Geometry::affine does not make for arbitrary polygons." );
      }

      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        DUNE_THROW( InvalidStateException, "Geometry::global does not make for arbitrary polygons." );
      }

      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        DUNE_THROW( InvalidStateException, "Geometry::local does not make for arbitrary polygons." );
      }

      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        DUNE_THROW( InvalidStateException, "Geometry::integrationElement does not make for arbitrary polygons." );
      }

      JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
      {
        DUNE_THROW( InvalidStateException, "Geometry::jacobianTransposed does not make for arbitrary polygons." );
      }

      JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        DUNE_THROW( InvalidStateException, "Geometry::jacobianInverseTransposed does not make for arbitrary polygons." );
      }

    private:
      static ctype det ( const GlobalCoordinate &x, const GlobalCoordinate &y, const GlobalCoordinate &z )
      {
        return x[ 0 ]*y[ 1 ]*z[ 2 ] + x[ 1 ]*y[ 2 ]*z[ 0 ] + x[ 2 ]*y[ 0 ]*z[ 1 ] - x[ 2 ]*y[ 1 ]*z[ 0 ] - x[ 0 ]*y[ 2 ]*z[ 1 ] - x[ 1 ]*y[ 0 ]*z[ 2 ];
      }

      GlobalCoordinate normal_;
      const GlobalCoordinate *cbegin_;
      std::size_t csize_;
    }



    // InterfaceGridGeometry
    // ---------------------

    template< int mydim, int cdim, class Grid >
    class InterfaceGridGeometry;

    template< int cdim, class Grid >
    class InterfaceGridGeometry< 0, cdim, Grid >
      : public AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 0, cdim >
    {
      typedef InterfaceGridGeometry< 0, cdim, Grid > This;
      typedef AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 0, cdim > Base;

    public:
      typedef typename Base::ctype ctype;

      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      InterfaceGridGeometry ( const GlobalCoordinate *cbegin, std::size_t csize )
        : Base( ReferenceElements< ctype, 0 >::cube(), cbegin[ 0 ], {} )
      {
        assert( csize == 1u );
      }
    };

    template< int cdim, class Grid >
    class InterfaceGridGeometry< 1, cdim, Grid >
      : public AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 1, cdim >
    {
      typedef InterfaceGridGeometry< 1, cdim, Grid > This;
      typedef AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 1, cdim > Base;

    public:
      typedef typename Base::ctype ctype;

      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      InterfaceGridGeometry ( const GlobalCoordinate *cbegin, std::size_t csize )
        : Base( ReferenceElements< ctype, 1 >::cube(), cbegin[ 0 ], { cbegin[ 1 ] - cbegin[ 0 ] } ) 
      {
        assert( csize == 2u );
      }

      InterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize )
        : This( cbegin, csize )
      {}
    };

    template< int cdim, class Grid >
    class InterfaceGridGeometry< 2, cdim, Grid >
      : public PolygonGeometry< typename std::remove_const_t< Grid >::Traits::ctype, cdim >
    {
      typedef InterfaceGridGeometry< 2, cdim, Grid > This;
      typedef PolygonGeometry< typename std::remove_const_t< Grid >::Traits::ctype, cdim > Base;

    public:
      typedef typename Base::ctype ctype;

      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      InterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize )
        : Base( normal, cbegin, csize )
      {}
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GEOMETRY_HH
