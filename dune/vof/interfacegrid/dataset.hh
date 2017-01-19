#ifndef DUNE_VOF_INTERFACEGRID_DATASET_HH
#define DUNE_VOF_INTERFACEGRID_DATASET_HH

#include <cassert>
#include <cstddef>

#include <utility>
#include <vector>

#include <dune/geometry/dimension.hh>

#include <dune/vof/flags.hh>
#include <dune/vof/interfacegrid/geometry.hh>
#include <dune/vof/mixedcellmapper.hh>
#include <dune/vof/utility.hh>

namespace Dune
{

  namespace VoF
  {

    // BasicInterfaceGridDataSet
    // -------------------------

    template< class R >
    struct BasicInterfaceGridDataSet
    {
      typedef R Reconstruction;

      typedef typename Reconstruction::GridView GridView;

      typedef VoF::Flags< GridView > Flags;

      typedef typename Reconstruction::ColorFunction ColorFunction;
      typedef typename Reconstruction::ReconstructionSet ReconstructionSet;

      typedef typename GridView::template Codim< 0 >::Entity Element;

      typedef typename ReconstructionSet::DataType::Coordinate GlobalCoordinate;

      template< class... Args >
      explicit BasicInterfaceGridDataSet ( const ColorFunction &colorFunction, Args &&... args )
        : reconstruction_( std::forward< Args >( args )... ),
          flags_( colorFunction.gridView() ),
          reconstructionSet_( colorFunction.gridView() )
      {
        update( colorFunction );
      }

      const Reconstruction &reconstruction () { return reconstruction_; }
      const Flags &flags () const { return flags_; }
      const ReconstructionSet &reconstructionSet () const { return reconstructionSet_; }

      const GridView &gridView () const { return flags().gridView(); }

      void update ( const ColorFunction &colorFunction )
      {
        flags_.reflag( colorFunction, 1e-6 );
        reconstruction_( colorFunction, reconstructionSet_, flags_ );
      }

      const GlobalCoordinate &normal ( const Element &element ) const { return reconstructionSet()[ element ].innerNormal(); }

    protected:
      Reconstruction reconstruction_;
      Flags flags_;
      ReconstructionSet reconstructionSet_;
    };



    // InterfaceGridDataSet
    // --------------------

    template< class R >
    class InterfaceGridDataSet
      : public BasicInterfaceGridDataSet< R >
    {
      typedef InterfaceGridDataSet< R > This;
      typedef BasicInterfaceGridDataSet< R > Base;

    public:
      using Base::flags;
      using Base::reconstructionSet;
      using Base::normal;

      typedef typename Base::ColorFunction ColorFunction;
      typedef typename Base::Element Element;
      typedef typename Base::GridView GridView;
      typedef typename Base::GlobalCoordinate GlobalCoordinate;

      typedef MixedCellMapper< GridView > Indices;
      typedef std::vector< GlobalCoordinate > Vertices;
      typedef std::vector< std::size_t > Offsets;

      typedef typename GridView::ctype ctype;

      template< int mydim >
      using Geometry = BasicInterfaceGridGeometry< ctype, mydim, GridView::dimensionworld >;

      template< class... Args >
      explicit InterfaceGridDataSet ( const ColorFunction &colorFunction, Args &&... args )
        : Base( colorFunction, std::forward< Args >( args )... ), indices_( flags() )
      {
        getInterfaceVertices( reconstructionSet(), flags(), vertices_, offsets_ );
      }

      void update ( const ColorFunction &colorFunction )
      {
        Base::update( colorFunction );
        indices_.update( flags() );
        getInterfaceVertices( reconstructionSet(), flags(), vertices_, offsets_ );
      }

      void covariantOuterNormal ( const Element &element, std::size_t i, FieldVector< ctype, 2 > &n ) const
      {
        const auto elementIndex = indices().index( element );
        const std::size_t index = offsets()[ elementIndex ];
        const std::size_t size = offsets()[ elementIndex + 1 ] - index;
        assert( i < size );
        n = vertices()[ index + i ] - vertices()[ index + (i + 1) % size ];
      }

      void covariantOuterNormal ( const Element &element, std::size_t i, FieldVector< ctype, 3 > &n ) const
      {
        const auto elementIndex = indices().index( element );
        const std::size_t index = offsets()[ elementIndex ];
        const std::size_t size = offsets()[ elementIndex + 1 ] - index;
        assert( i < size );
        n = generalizedCrossProduct( vertices()[ index + (i + 1) % size ] - vertices()[ index + i ], normal( element ) );
      }

      Geometry< 0 > geometry ( const Element &element, std::size_t i, Dune::Dim< 0 > ) const
      {
        const auto elementIndex = indices().index( element );
        const std::size_t index = offsets()[ elementIndex ];
        assert( index + i < offsets()[ elementIndex+1 ] );
        return Geometry< 0 >( vertices()[ index + i ] );
      }

      Geometry< 1 > geometry ( const Element &element, std::size_t i, Dune::Dim< 1 > ) const
      {
        const auto elementIndex = indices().index( element );
        const std::size_t index = offsets()[ elementIndex ];
        const std::size_t size = offsets()[ elementIndex + 1 ] - index;
        assert( i < size );
        return Geometry< 1 >( vertices()[ index + i ], vertices()[ index + (i + 1) % size ] );
      }

      std::size_t numVertices ( const Element &element ) const
      {
        const auto elementIndex = indices().index( element );
        return (offsets()[ elementIndex + 1 ] - offsets()[ elementIndex ]);
      }

      const Indices &indices () const { return indices_; }
      const Vertices &vertices () const { return vertices_; }
      const Offsets &offsets () const { return offsets_; }

    private:
      Indices indices_;
      Vertices vertices_;
      Offsets offsets_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_DATASET_HH
