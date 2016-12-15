#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>
#include <utility>
#include <chrono>

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
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/operator/projection/l2projection.hh>

// dune-vof includes
#include <dune/vof/femdfwrapper.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstructedfunction.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

// local includes
#include "binarywriter.hh"
#include "polygon.hh"
#include "reconstructionwriter.hh"
#include "velocity.hh"
#include "../problems/linearwall.hh"
#include "../problems/rotatingcircle.hh"
#include "../problems/sflow.hh"
#include "../problems/slottedcylinder.hh"



// Problem
// -------
template < class Base, class FunctionSpace >
struct Problem : public Base
{
  using FunctionSpaceType = FunctionSpace;
};

template < class GridPart, class Velocity, class Flags >
double initTimeStep( const GridPart& gridPart, Velocity &velocity, const Flags &flags )
{
  double dtMin = std::numeric_limits< double >::max();
  for( const auto &entity : elements( gridPart ) )
  {
    if( !flags.isActive( entity ) )
    continue;

    const auto geoEn = entity.geometry();

    for ( const auto &intersection : intersections( gridPart, entity ) )
    {
      if ( !intersection.neighbor() )
        continue;

      const auto geoIs = intersection.geometry();
      velocity.bind( intersection );
      auto v = velocity( geoIs.local( geoIs.center() ) );
      dtMin = std::min( dtMin, geoEn.volume() / std::abs( intersection.integrationOuterNormal( typename decltype( geoIs )::LocalCoordinate( 0 ) ) * v ) );
    }
  }
  return dtMin;
}



// Algorithm
// ---------

template< class Grid, class DF, class P >
double algorithm ( Grid &grid, DF& uh, P& problem, int level, double start, double end, double cfl, double eps, const bool writeData = true )
{
  using GridType = Grid;
  using GridPartType = Dune::Fem::LeafGridPart< GridType >;
  GridPartType gridPart( grid );

  // Create stencils
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridPartType >;
  Stencils stencils( gridPart );

  // Create discrete functions
  using FunctionSpaceType = Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;
  using DiscreteFunctionSpaceType = Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  DiscreteFunctionSpaceType space( gridPart );

  using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;
  DiscreteFunctionType update( "update", space );

  auto cuh = Dune::VoF::discreteFunctionWrapper( uh );
  auto cupdate = Dune::VoF::discreteFunctionWrapper( update );

  // Create reconstruction set
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridPartType >;
  ReconstructionSet reconstructions( gridPart );

  // Create operators
  auto reconstruction = Dune::VoF::reconstruction( gridPart, cuh, stencils );
  auto flags = Dune::VoF::flags( gridPart );
  auto evolution = Dune::VoF::evolution( gridPart );

  // Calculate initial data
  flags.reflag( cuh, eps );
  reconstruction( cuh, reconstructions, flags );

  using Velocity = Velocity< P, typename GridPartType::GridViewType >;
  Velocity velocity( problem, start );

  // Create and initialize time provider
  using TimeProviderType = Dune::Fem::TimeProvider< typename GridType::CollectiveCommunication >;
  TimeProviderType timeProvider( start, cfl, gridPart.comm() );
  timeProvider.init( initTimeStep( gridPart, velocity, flags ) );

  // Create data output
  using DataOutputType = BinaryWriter;
  DataOutputType dataOutput( level, timeProvider );

  // Calculate and write initial data
  if ( writeData ) dataOutput.write( grid, uh, timeProvider );



  // Time Iteration
  for( ; timeProvider.time() <= end; )
  {
    if( Dune::Fem::Parameter::verbose() )
      std::cerr << "time step = " << timeProvider.timeStep() << ", "
                << "time = " << timeProvider.time() << ", "
                << "dt = " << timeProvider.deltaT() << std::endl;

                // Create velocity object
    Velocity velocity( problem, timeProvider.time() );

    flags.reflag( cuh, eps );
    reconstruction( cuh, reconstructions, flags );
    double dtEst = evolution( reconstructions, flags, velocity, timeProvider.deltaT(), cupdate );

    update.communicate();
    uh.axpy( 1.0, update );

    timeProvider.provideTimeStepEstimate( dtEst );
    timeProvider.next();

    if ( writeData ) dataOutput.write( grid, uh, timeProvider );
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
  double startTime = Dune::Fem::Parameter::getValue< double >( "start", 0.0 );
  const double endTime = Dune::Fem::Parameter::getValue< double >( "end", 2.5 );
  const double cfl = Dune::Fem::Parameter::getValue< double >( "cfl", 1.0 );
  const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );
  const std::string path = Dune::Fem::Parameter::getValue< std::string >( "fem.io.path", "data" );
  const int restartStep = Dune::Fem::Parameter::getValue< int >( "restartStep", -1 );

  using GridType = Dune::GridSelector::GridType;
  using GridPartType = Dune::Fem::LeafGridPart< GridType >;
  using FunctionSpaceType = Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;
  using DiscreteFunctionSpaceType = Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;

  // Testproblem
  using ProblemType = Problem< SFlow< double, GridPartType::dimensionworld >, FunctionSpaceType >;
  using SolutionType = Dune::Fem::InstationaryFunction< ProblemType >;
  using GridSolutionType = Dune::Fem::GridFunctionAdapter< SolutionType, GridPartType >;


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
  double oldError = 0;
  for( int step = level; step <= level+repeats; ++step )
  {
    // Initialize Data
    // ===============
    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType space( gridPart );
    ProblemType problem;
    DiscreteFunctionType uh( "uh", space );

    SolutionType solution( problem, startTime );
    GridSolutionType u( "solution", solution, gridPart, 19 );
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
      namedata << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-p" << std::setw(4) << Dune::Fem::MPIManager::rank()
        << "-vof-fem-" << std::to_string( step ) << "-" << std::setw(5) << restartStep << ".bin";
      const auto filename = Dune::concatPaths( path, namedata.str() );
      if ( !Dune::Fem::fileExists( filename ) )
      {
        std::cout << "Restart error: A step is given to load binary file but this file does not exist." << std::endl;
        return 1;
      }
      else
      {
        Dune::Fem::BinaryFileInStream binaryStream ( filename );
        double timestamp;
        binaryStream >> timestamp;
        startTime = timestamp;
        uh.read( binaryStream );
        std::cout << "Restarted in file " << filename << std::endl;
      }
    }

    uh.communicate();


    // Run Algorithm
    double actualEndTime = algorithm( grid, uh, problem, step, startTime, endTime, cfl, eps );


    // Calculate Error
    solution.setTime( actualEndTime );
    Dune::VoF::ReconstructedFunction< DiscreteFunctionType > uRec ( uh );

    using Quadrature = Dune::Fem::CachingQuadrature < GridPartType, 0 >;
    double error = 0.0;
    for ( const auto& entity : elements( gridPart, Dune::Partitions::interior ) )
    {
      Quadrature quad ( entity, 19 );
      for ( const auto& qp : quad )
      {
        const auto x = entity.geometry().global( qp.position() );
        DiscreteFunctionType::RangeType v, w;
        solution.evaluate( x, v );
        uRec.evaluate( entity, x, w );
        error += std::abs( v - w ) * qp.weight() * entity.geometry().integrationElement( qp.position() );
      }
    }

    if ( Dune::Fem::MPIManager::rank() == 0 )
    {
      const double l1eoc = log( oldError / error ) / M_LN2;
      std::cout << "L1-Error =" << std::setw( 16 ) << error << std::endl;

      if ( step > level )
        std::cout << "L1 EOC( " << std::setw( 2 ) << step << " ) = " << std::setw( 11 ) << l1eoc << std::endl;
    }

    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );
    oldError = error;
  }

  return 0;
}
catch ( Dune::Exception &e ) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...) {
  std::cerr << "Unknown exception thrown!" << std::endl;
}
