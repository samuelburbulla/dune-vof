#ifndef DUNE_VOF_INTERFACEGRID_GEOMETRY_HH
#define DUNE_VOF_INTERFACEGRID_GEOMETRY_HH

#include <type_traits>

#include <dune/geometry/affinegeometry.hh>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  namespace VoF
  {

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
      InterfaceGridGeometry ( const HostGeometry &hostGeometry )
      : Base( hostGeometry )
      {}
    };

    template< int cdim, class Grid >
    class InterfaceGridGeometry< 1, cdim, Grid >
      : public AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 1, cdim >
    {
      typedef InterfaceGridGeometry< 1, cdim, Grid > This;
      typedef AffineGeometry< typename std::remove_const_t< Grid >::Traits::ctype, 1, cdim > Base;

    public:
      InterfaceGridGeometry ( const HostGeometry &hostGeometry )
      : Base( hostGeometry )
      {}
    };

    template< int cdim, class Grid >
    class InterfaceGridGeometry< 2, cdim, Grid >
    {
      typedef InterfaceGridGeometry< 2, cdim, Grid > This;

    public:
      typedef ct ctype;
      
      static const int mydimension = 2;
      static const int coorddimension = cdim;
      
      typedef FieldVector< ctype, mydimension > LocalCoordinate;
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;
      
      typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;
      typedef FieldMatrix< ctype, coorddimension, mydimension > JacobianInverseTransposed;

      InterfaceGridGeometry ( ... )
      {
        // TODO: Please implement me
      }

      int corners () const
      {
        // TODO: Please implement me
      }

      GlobalCoordinate corner ( int i ) const
      {
        // TODO: Please implement me
      }

      GlobalCoordinate center () const
      {
        // TODO: Please implement me
      }

      ctype volume () const
      {
        // TODO: Please implement me
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
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_GEOMETRY_HH
