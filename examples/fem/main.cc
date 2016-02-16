#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>
#include <utility>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-grid include
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/common/instationary.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/operator/projection/l2projection.hh>

// dune-vof includes
#include <dune/vof/femdfwrapper.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

// local includes
#include "binarywriter.hh"
#include "polygon.hh"
#include "reconstructionwriter.hh"
#include "../problems/rotatingcircle.hh"



// Problem
// -------
template < class Base, class FunctionSpace >
struct Problem : public Base
{
  using FunctionSpaceType = FunctionSpace;
};


// algorithm
// ---------

template< class Grid, class DF, class P >
const double algorithm ( Grid &grid, DF& uh, P& problem, int level, double start, double end, double cfl, double eps )
{
  using GridType = Grid;
  using ProblemType = P;
  using GridPartType = Dune::Fem::LeafGridPart< GridType >;

  using DomainType = typename Dune::FieldVector < double, GridType::dimensionworld >;

  using FunctionSpaceType = Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;
  using DiscreteFunctionSpaceType = Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;

  using TimeProviderType = Dune::Fem::FixedStepTimeProvider< typename GridType::CollectiveCommunication >;

  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridPartType >;
  using Polygon = Polygon< typename ReconstructionSet::Reconstruction::Coordinate >;

  using DataOutputType = BinaryWriter;
  using RecOutputType = ReconstructionWriter< ReconstructionSet, Polygon >;

  using Stencils = Dune::VoF::VertexNeighborsStencil< GridPartType >;

  GridPartType gridPart( grid );

  Stencils stencils( gridPart );
  ReconstructionSet reconstructions( gridPart );

  DiscreteFunctionSpaceType space( gridPart );
  DiscreteFunctionType update( "update", space );

  auto cuh = Dune::VoF::discreteFunctionWrapper( uh );
  auto cupdate = Dune::VoF::discreteFunctionWrapper( update );

  auto evolution = Dune::VoF::evolution( reconstructions, cuh, eps );
  auto reconstruction = Dune::VoF::reconstruction( gridPart, cuh, stencils );
  auto flags = Dune::VoF::flags( gridPart );

  double timeStep = std::pow( 2, - ( 3 + level ) );
  timeStep *= cfl;
  timeStep /= ProblemType::maxVelocity();
  TimeProviderType timeProvider( start, timeStep, gridPart.comm() );

  DataOutputType dataOutput( level, timeProvider );
  RecOutputType recOutput ( level );

  auto velocity = [ &timeProvider, &problem ] ( const auto &x ) { DomainType u; problem.velocityField( x, timeProvider.time(), u ); return u; };

  // Time Iteration
  // ==============
  flags.reflag( cuh, eps );
  reconstruction( cuh, reconstructions, flags );

  recOutput.write ( reconstructions );
  dataOutput.write( uh );

  std::size_t count = 0;
  for( ; timeProvider.time() <= end; )
  {
    if( Dune::Fem::Parameter::verbose() )
      std::cerr << "time step = " << timeProvider.timeStep() << ", "
                << "time = " << timeProvider.time() << ", "
                << "dt = " << timeProvider.deltaT() << std::endl;

    flags.reflag( cuh, eps );
    reconstruction( cuh, reconstructions, flags );

    evolution( cuh, reconstructions, velocity, timeProvider.deltaT(), cupdate, flags );
    update.communicate();
    uh.axpy( 1.0, update );

    timeProvider.next();

    if ( dataOutput.willWrite( timeProvider ) )
    {
      recOutput.write ( reconstructions );
      dataOutput.write( uh );

      if( Dune::Fem::Parameter::verbose() )
          std::cout << "written reconstructions count=" << count << std::endl;
    }
  }

  return timeProvider.time();
}


int main(int argc, char** argv)
try {

  Dune::Fem::MPIManager::initialize( argc, argv );


  // Read Parameter File
  // ===================
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  const int level = Dune::Fem::Parameter::getValue< int >( "level", 0 );
  const int repeats = Dune::Fem::Parameter::getValue< int >( "repeats", 2 );
  const double startTime = Dune::Fem::Parameter::getValue< double >( "start", 0.0 );
  const double endTime = Dune::Fem::Parameter::getValue< double >( "end", 2.5 );
  const double cfl = Dune::Fem::Parameter::getValue< double >( "cfl", 1.0 );
  const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );
  const std::string path = Dune::Fem::Parameter::getValue< typename std::string >( "fem.io.path", "./data/" );
  const int restartStep = Dune::Fem::Parameter::getValue< int >( "restartStep", -1 );

  using GridType = Dune::GridSelector::GridType;
  using GridPartType = Dune::Fem::LeafGridPart< GridType >;
  using FunctionSpaceType = Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;
  using DiscreteFunctionSpaceType = Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;

  // Testproblem
  using ProblemType = Problem< RotatingCircle< double, GridPartType::dimensionworld >, FunctionSpaceType >;
  using SolutionType = Dune::Fem::InstationaryFunction< ProblemType >;
  using GridSolutionType = Dune::Fem::GridFunctionAdapter< SolutionType, GridPartType >;
  using L1NormType = Dune::Fem::L1Norm< GridPartType >;


  // Create Grid
  // ===========
  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;
  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );


  // EOC Calculation
  // ===============
  double oldL1Error;
  for( int step = level; step <= level+repeats; ++step )
  {
    // Initialize Data
    // ===============
    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType space( gridPart );
    ProblemType problem;
    DiscreteFunctionType uh( "uh", space );
    L1NormType l1norm( gridPart );

    SolutionType solution( problem, startTime );
    GridSolutionType u( "solution", solution, gridPart, 9 );
    Dune::Fem::L2Projection < GridSolutionType, DiscreteFunctionType > l2projection ( 15 );

    if ( restartStep == -1 )
    {
      // Use initial data of problem.
      l2projection( u, uh );
    }
    else
    {
      // Use given data in binary file.
      std::stringstream namedata;
      namedata.fill('0');
      namedata << "vof-fem-" << std::to_string( step ) << "-" << std::setw(5) << restartStep << ".bin";
      const auto filename = Dune::concatPaths( path, namedata.str() );
      if ( !Dune::Fem::fileExists( filename ) )
      {
        std::cout << "Restart error: A step is given to load binary file but this file does not exist." << std::endl;
        return 1;
      }
      else
      {
        Dune::Fem::BinaryFileInStream binaryStream ( filename );
        uh.read( binaryStream );
      }
    }

    uh.communicate();


    // Run Algorithm
    // ===============
    const double actualEndTime = algorithm( grid, uh, problem, step, startTime, endTime, cfl, eps );


    // Calculate Error
    // ===============
    solution.setTime( actualEndTime );
    DiscreteFunctionType uhEnd( "", space );
    l2projection( u, uhEnd );
    const double newL1Error = l1norm.distance( uh, uhEnd );

    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      const double l1eoc = log( oldL1Error / newL1Error ) / M_LN2;
      std::cout << "L1 Error: " << std::setw( 16 ) << newL1Error << std::endl;

      if ( step > level )
        std::cout << "L1 EOC( " << std::setw( 2 ) << step << " ) = " << std::setw( 11 ) << l1eoc << std::endl;
    }

    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );
    oldL1Error = newL1Error;
  }

  return 0;
}
catch ( Dune::Exception &e ) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...) {
  std::cerr << "Unknown exception thrown!" << std::endl;
}
