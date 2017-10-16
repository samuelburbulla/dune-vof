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

    // BasicPolygonGeometry
    // --------------------

    template< class ct, int cdim, class SignedTrapezoidArea >
    struct BasicPolygonGeometry
    {
      typedef ct ctype;

      static const int mydimension = 2;
      static const int coorddimension = cdim;

      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

      template< class... Args >
      BasicPolygonGeometry ( const GlobalCoordinate *cbegin, std::size_t csize, Args &&... args )
        : cbegin_( cbegin ), csize_( csize ), signedTrapezoidArea_( std::forward< Args >( args )... )
      {
        assert( csize >= 3u );
      }

      int corners () const { return csize_; }

      const GlobalCoordinate &corner ( int i ) const { assert( (i >= 0) && (i < corners()) ); return cbegin_[ i ]; }

      GlobalCoordinate center () const
      {
        GlobalCoordinate center( 0 );
        ctype volume = 0;
        for( std::size_t i = 0; i < csize_; ++i )
        {
          const GlobalCoordinate &x = cbegin_[ i ];
          const GlobalCoordinate &y = cbegin_[ (i+1) % csize_ ];
          const ctype weight = signedTrapezoidArea_( x, y );
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
          volume += signedTrapezoidArea_( cbegin_[ i ], cbegin_[ (i+1) % csize_ ] );
        return abs( volume / ctype( 2 ) );
      }

      GeometryType type () const noexcept { return GeometryTypes::none( mydimension ); }

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
      const GlobalCoordinate *cbegin_;
      std::size_t csize_;
      SignedTrapezoidArea signedTrapezoidArea_;
    };



    // SignedTrapezoidArea
    // -------------------

    template< class ctype, int dim >
    struct SignedTrapezoidArea;

    template< class ctype >
    struct SignedTrapezoidArea< ctype, 2 >
    {
      ctype operator() ( const FieldVector< ctype, 2 > &x, const FieldVector< ctype, 2 > &y ) const { return det( x, y ); }

    private:
      static ctype det ( const FieldVector< ctype, 2 > &x, const FieldVector< ctype, 2 > &y )
      {
        return x[ 0 ]*y[ 1 ] - x[ 1 ]*y[ 0 ];
      }
    };

    template< class ctype >
    struct SignedTrapezoidArea< ctype, 3 >
    {
      explicit SignedTrapezoidArea ( const FieldVector< ctype, 3 > &normal ) : normal_( normal ) {}

      ctype operator() ( const FieldVector< ctype, 3 > &x, const FieldVector< ctype, 3 > &y ) const { return det( x, y, normal_ ); }

    private:
      static ctype det ( const FieldVector< ctype, 3 > &x, const FieldVector< ctype, 3 > &y, const FieldVector< ctype, 3 > &z )
      {
        return x[ 0 ]*y[ 1 ]*z[ 2 ] + x[ 1 ]*y[ 2 ]*z[ 0 ] + x[ 2 ]*y[ 0 ]*z[ 1 ] - x[ 2 ]*y[ 1 ]*z[ 0 ] - x[ 0 ]*y[ 2 ]*z[ 1 ] - x[ 1 ]*y[ 0 ]*z[ 2 ];
      }

      FieldVector< ctype, 3 > normal_;
    };



    // PolygonGeometry
    // ---------------

    template< class ct, int cdim >
    using PolygonGeometry = BasicPolygonGeometry< ct, cdim, SignedTrapezoidArea< ct, cdim > >;



    // BasicInterfaceGridGeometry
    // --------------------------

    template< class ctype, int mydim, int cdim >
    class BasicInterfaceGridGeometry;

    template< class ctype, int cdim >
    class BasicInterfaceGridGeometry< ctype, 0, cdim >
      : public AffineGeometry< ctype, 0, cdim >
    {
      typedef BasicInterfaceGridGeometry< ctype, 0, cdim > This;
      typedef AffineGeometry< ctype, 0, cdim > Base;

    public:
      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      explicit BasicInterfaceGridGeometry ( const GlobalCoordinate &x )
        : Base( ReferenceElements< ctype, 0 >::cube(), x, {} )
      {}

      explicit BasicInterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize  )
        : Base( ReferenceElements< ctype, 0 >::cube(), cbegin[ 0 ], {} )
      {}
    };

    template< class ctype, int cdim >
    class BasicInterfaceGridGeometry< ctype, 1, cdim >
      : public AffineGeometry< ctype, 1, cdim >
    {
      typedef BasicInterfaceGridGeometry< ctype, 1, cdim > This;
      typedef AffineGeometry< ctype, 1, cdim > Base;

    public:
      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      BasicInterfaceGridGeometry ( const GlobalCoordinate &x, const GlobalCoordinate &y )
        : Base( ReferenceElements< ctype, 1 >::cube(), x, { y - x } )
      {}

      BasicInterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize )
        : This( cbegin[ 0 ], cbegin[ 1 ] )
      {
        assert( csize == 2u );
      }
    };

    template< class ctype, int cdim >
    class BasicInterfaceGridGeometry< ctype, 2, cdim >
      : public PolygonGeometry< ctype, cdim >
    {
      typedef BasicInterfaceGridGeometry< ctype, 2, cdim > This;
      typedef PolygonGeometry< ctype, cdim > Base;

    public:
      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      BasicInterfaceGridGeometry ( const GlobalCoordinate *cbegin, std::size_t csize )
        : Base( cbegin, csize )
      {}

      BasicInterfaceGridGeometry ( const GlobalCoordinate &normal, const GlobalCoordinate *cbegin, std::size_t csize )
        : Base( cbegin, csize, normal )
      {}
    };



    // InterfaceGridGeometry
    // ---------------------

    template< int mydim, int cdim, class Grid >
    using InterfaceGridGeometry = BasicInterfaceGridGeometry< typename std::remove_const_t< Grid >::Traits::ctype, mydim, cdim >;

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GEOMETRY_HH
