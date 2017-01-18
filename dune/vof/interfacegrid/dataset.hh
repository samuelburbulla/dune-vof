#ifndef DUNE_VOF_INTERFACEGRID_DATASET_HH
#define DUNE_VOF_INTERFACEGRID_DATASET_HH

#include <cstddef>
#include <utility>
#include <vector>

#include <dune/vof/flags.hh>
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

      typedef MixedCellMapper< typename R::GridView > Indices;
      typedef std::vector< typename ReconstructionSet::DataType::Coordinate > Vertices;
      typedef std::vector< std::size_t > Offsets;

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
