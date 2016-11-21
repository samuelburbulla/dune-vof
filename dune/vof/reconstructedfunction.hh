#ifndef DUNE_VOF_RECONSTRUCTEDFUNCTION_HH
#define DUNE_VOF_RECONSTRUCTEDFUNCTION_HH

#include <dune/vof/femdfwrapper.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>


namespace Dune
{
  namespace VoF
  {

    // ReconstructionFunction
    // ----------------------

    /**
     * \ingroup Other
     * \brief reconstructed function
     *
     * \tparam  DF  discrete function
     */
    template< class DF >
    struct ReconstructedFunction
    {
      using GridPartType = typename DF::GridPartType;
      using DomainType = typename DF::DomainType;
      using RangeType = typename DF::RangeType;
      using Flags = typename Dune::VoF::Flags< GridPartType >;
      using ReconstructionSet = typename Dune::VoF::ReconstructionSet< GridPartType >;
      using Stencils = Dune::VoF::VertexNeighborsStencil< GridPartType >;

    public:
      ReconstructedFunction ( const DF &uh )
       : uh_( uh ), flags_( Dune::VoF::flags( uh.gridPart() ) ), reconstructions_( uh.gridPart() ), stencils_ ( uh.gridPart() )
      {
        auto cuh = Dune::VoF::discreteFunctionWrapper( uh );
        auto reconstruction = Dune::VoF::reconstruction( uh.gridPart(), cuh, stencils_ );
        const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );

        flags_.reflag( cuh, eps );
        reconstruction( cuh, reconstructions_, flags_ );
      }

      template< class Entity >
      void evaluate ( const Entity &entity, const DomainType &x, RangeType &u ) const
      {
        if ( flags_.isFull( entity ) )
          u = 1.0;
        else if ( flags_.isMixed( entity ) )
          u = ( reconstructions_[ entity ].levelSet( x ) > 0.0 );
        else
          u = 0.0;
      }

    private:
      const DF uh_;
      Flags flags_;
      ReconstructionSet reconstructions_;
      Stencils stencils_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTEDFUNCTION_HH
