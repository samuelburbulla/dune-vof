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

      typedef VoF::Flags< typename Reconstruction::GridView > Flags;

      typedef typename Reconstruction::ColorFunction ColorFunction;
      typedef typename Reconstruction::ReconstructionSet ReconstructionSet;

      template< class... Args >
      explicit BasicInterfaceGridDataSet ( const ColorFunction &colorFunction, Args &&... args )
        : reconstruction_( std::forward< Args >( args )... ),
          flags_( reconstruction_.gridView(), colorFunction, 1e-6 ),
          reconstructionSet_( reconstruction_.gridView() )
      {
        update( colorFunction );
      }

      const Reconstruction &reconstruction () { return reconstruction_; }
      const Flags &flags () const { return flags_; }
      const ReconstructionSet &reconstructionSet () const { return reconstructionSet_; }

      void update ( const ColorFunction &colorFunction )
      {
        flags_.reflag( colorFunction, 1e-6 );
        reconstruction_( colorFunction, reconstructionSet_, flags_ );
      }

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

      typedef typename Base::ColorFunction ColorFunction;

      typedef MixedCellMapper< typename R::GridView > Indices;
      typedef std::vector< typename ReconstructionSet::DataType::Coordinate > Vertices;
      typedef typename R::GridView::template Codim< 0 >::Entity Element;
      typedef std::vector< std::size_t > Offsets;

      template< int mydim >
      using Geometry = BasicInterfaceGridGeometry< typename R::GridView::ctype, mydim, R::dimensionworld >;

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
