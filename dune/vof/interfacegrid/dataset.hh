#ifndef DUNE_VOF_INTERFACEGRID_DATASET_HH
#define DUNE_VOF_INTERFACEGRID_DATASET_HH

#include <cstddef>
#include <utility>
#include <vector>

#include <dune/vof/flags.hh>
#include <dune/vof/utility.hh>

namespace Dune
{

  namespace VoF
  {

    // InterfaceGridDataSet
    // --------------------

    template< class R >
    struct InterfaceGridDataSet
    {
      typedef R Reconstruction;

      typedef VoF::Flags< typename Reconstruction::GridView > Flags;

      typedef typename Reconstruction::ColorFunction ColorFunction;

      template< class... Args >
      explicit InterfaceGridDataSet ( const ColorFunction &colorFunction, Args &&... args )
        : reconstruction_( std::forward< Args >( args )... ),
          flags_( reconstruction_.gridView(), colorFunction, 1e-6 ),
          reconstructionSet_( reconstruction_.gridView() )
      {
        update( colorFunction );
      }

      const Reconstruction &reconstruction () { return reconstruction_; }
      const Flags &flags () const { return flags_; }

      void update ( const ColorFunction &colorFunction )
      {
        flags_.reflag( colorFunction, 1e-6 );
        reconstruction_( colorFunction, reconstructionSet_, flags_ );
        getInterfaceVertices( reconstructionSet_, flags_, vertices_, offsets_ );
      }

    protected:
      Reconstruction reconstruction_;
      Flags flags_;
      typename Reconstruction::ReconstructionSet reconstructionSet_;
      std::vector< typename ReconstructionSet::DataType::Coordinate > vertices_;
      std::vector< std::size_t > offsets_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_INTERFACEGRID_DATASET_HH
