#ifndef DUNE_VOF_ALGORITHM_HH
#define DUNE_VOF_ALGORITHM_HH

// C++ includes
#include <utility>

// dune-vof includes
#include "test/errors.hh"
#include <dune/vof/evolution.hh>
#include <dune/vof/flagset.hh>
#include <dune/vof/flagging.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionset.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/velocity.hh>

namespace Dune
{
  namespace VoF
  {

    // Algorithm
    // ---------

    /**
     * \ingroup   Method
     * \brief     volume of fluid evolution algorithm
     * \details   Rider, W.J., Kothe, D.B., Reconstructing Volume Tracking, p. 15ff
     *
     * \tparam  GV  grid view type
     * \tparam  PR  problem type
     * \tparam  DW  data writer type
     */
    template< class GV, class PR, class DW >
    struct Algorithm
    {
      using GridView = GV;
      using Problem = PR;
      using DataWriter = DW;
      using Stencils = VertexNeighborsStencil< GridView >;
      using Reconstructions = ReconstructionSet< GridView >;
      using Flags = FlagSet< GridView >;
      using VelocityField = Velocity< Problem, GridView >;


      Algorithm ( const GridView &gridView, const Problem& problem, DataWriter& dataWriter, double cfl, double eps, const bool verbose = false )
       : gridView_( gridView ), problem_( problem ), dataWriter_( dataWriter ), cfl_( cfl ), eps_( eps ), verbose_( verbose ),
         stencils_( gridView ), reconstructions_( gridView ), flags_( gridView )
      {}

      template< class ColorFunction >
      double operator() ( ColorFunction& uh, double start, double end, int level = 0 )
      {
        // Create operators
        auto reconstructionOperator = reconstruction( gridView_, uh, stencils_ );
        auto flagOperator = FlagOperator< ColorFunction, Flags >( eps_ );
        auto evolutionOperator = evolution( gridView_ );

        double time = start, dt = 0.0, dtEst = 0.0;
        double error = 0.0;
        ColorFunction update( gridView_ );

        // Time Iteration
        do
        {
          if ( verbose_ )
            std::cerr << "time = " << time << ", " << "dt = " << dt << std::endl;

          VelocityField velocity( problem_, time );

          flagOperator( uh, flags_ );
          reconstructionOperator( uh, reconstructions_, flags_ );

          double dtEstimate = evolutionOperator( reconstructions_, flags_, velocity, dt, update );
          dtEst = gridView_.comm().min( dtEstimate );

          update.communicate();
          uh.axpy( 1.0, update );

          if ( dt > 0.0 )
          {
            error += dt * Dune::VoF::l1error( gridView_, reconstructions(), flags(), problem_, time, level );
            time += dt;
          }

          dataWriter_.write( time );

          dt = dtEst * cfl_;
        }
        while( time < end );

        if ( time == start )
          error += Dune::VoF::l1error( gridView_, reconstructions(), flags(), problem_, time );

        return error;
      }

      const Reconstructions& reconstructions() const { return reconstructions_; }

      const Flags& flags() const { return flags_; }


    private:
      const GridView& gridView_;
      const Problem& problem_;
      DataWriter& dataWriter_;
      const double cfl_, eps_;
      const bool verbose_;
      const Stencils stencils_;
      Reconstructions reconstructions_;
      Flags flags_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_ALGORITHM_HH
