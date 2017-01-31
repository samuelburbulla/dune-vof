#ifndef DUNE_VOF_ALGORITHM_HH
#define DUNE_VOF_ALGORITHM_HH

// C++ includes
#include <utility>

// dune-vof includes
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
      void operator() ( ColorFunction& uh, double start, double &end )
      {
        // Create operators
        auto reconstructionOperator = reconstruction( gridView_, uh, stencils_ );
        auto flagOperator = FlagOperator< ColorFunction, Flags >( eps_ );
        auto evolutionOperator = evolution( gridView_ );

        double time = start, dt = 0.0, dtEst = 0.0;
        ColorFunction update( gridView_ );

        // Time Iteration
        for( ; time <= end; time += dt, dt = dtEst * cfl_)
        {
          if ( verbose_ )
            std::cerr << "time = " << time << ", " << "dt = " << dt << std::endl;

          VelocityField velocity( problem_, time );

          flagOperator( uh, flags_ );
          reconstructionOperator( uh, reconstructions_, flags_ );
          dtEst = evolutionOperator( reconstructions_, flags_, velocity, dt, update );

          gridView_.comm().min( dtEst );

          update.communicate();
          uh.axpy( 1.0, update );

          dataWriter_.write( time );
        }

        end = time;
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