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
#include <dune/fem/misc/mpimanager.hh>

// dune-vof includes
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

#include "../dune/vof/test/colorfunction.hh"
#include "../dune/vof/test/average.hh"
#include "../dune/vof/test/errors.hh"
#include "../dune/vof/test/velocity.hh"
#include "../dune/vof/test/binarywriter.hh"
#include "../dune/vof/test/polygon.hh"
#include "../dune/vof/test/reconstructionwriter.hh"
#include "../dune/vof/test/problems/linearwall.hh"
#include "../dune/vof/test/problems/rotatingcircle.hh"
#include "../dune/vof/test/problems/sflow.hh"
#include "../dune/vof/test/problems/slottedcylinder.hh"


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

template< class Grid, class ColorFunction, class P >
double algorithm ( Grid &grid, ColorFunction& uh, P& problem, int level, double start, double end, double cfl, double eps, const bool writeData = true )
{
  using GridView = typename Grid::LeafGridView;
  GridView gridView( grid.leafGridView() );

  // Create stencils
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
  Stencils stencils( gridView );

  // Create reconstruction set
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
  ReconstructionSet reconstructions( gridView );

  // Create operators
  auto reconstruction = Dune::VoF::reconstruction( gridView, uh, stencils );
  auto flags = Dune::VoF::flags( gridView );
  auto evolution = Dune::VoF::evolution( gridView );

  // Calculate initial data
  flags.reflag( uh, eps );
  reconstruction( uh, reconstructions, flags );

  using Velocity = Velocity< P, GridView >;
  Velocity velocity( problem, start );

  // Create and initialize time provider
  double time = start;
  double deltaT = initTimeStep( gridView, velocity, flags );

  // Create data output
  using DataOutputType = BinaryWriter;
  DataOutputType dataOutput( level, time );

  // Calculate and write initial data
  if ( writeData )
    dataOutput.write( grid, uh, time );

  ColorFunction update( gridView );

  // Time Iteration
  for( ; time <= end; )
  {
    if( Dune::Fem::Parameter::verbose() )
      std::cerr << "time = " << time<< ", "
                << "dt = " << deltaT << std::endl;

                // Create velocity object
    Velocity velocity( problem, time );

    flags.reflag( uh, eps );
    reconstruction( uh, reconstructions, flags );
    double dtEst = evolution( reconstructions, flags, velocity, deltaT, update );

    update.communicate();
    uh.axpy( 1.0, update );

    time += deltaT;
    deltaT = dtEst * cfl;

    if ( writeData )
      dataOutput.write( grid, uh, time );
  }

  return Dune::VoF::l1error( gridView, reconstructions, flags, problem, time );
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

  using Grid = Dune::GridSelector::GridType;
  using GridView = typename Grid::LeafGridView;


  // Create Grid
  // ===========
  std::stringstream gridFile;
  gridFile << "../dune/vof/test/" << Grid::dimension << "dgrid.dgf";

  Dune::GridPtr< Grid > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  Grid& grid = *gridPtr;
  const int refineStepsForHalf = Dune::DGFGridInfo< Grid >::refineStepsForHalf();
  grid.globalRefine( level * refineStepsForHalf );

  GridView gridView( grid.leafGridView() );

  using ColorFunction = ColorFunction< GridView >;

  // Testproblem
  using ProblemType = SFlow< double, GridView::dimensionworld >;
  ProblemType problem;


  // EOC Calculation
  // ===============
  double oldError = 0;
  for( int step = level; step <= level+repeats; ++step )
  {
    // Initialize Data
    // ===============
    ColorFunction uh( gridView );


    if ( restartStep == -1 )
    {
      // Use initial data of problem.
      Dune::VoF::average( uh, problem, 0.0 );
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
    double error = algorithm( grid, uh, problem, step, startTime, endTime, cfl, eps );


    if ( Dune::Fem::MPIManager::rank() == 0 )
    {
      const double l1eoc = log( oldError / error ) / M_LN2;
      std::cout << "L1-Error =" << std::setw( 16 ) << error << std::endl;

      if ( step > level )
        std::cout << "L1 EOC( " << std::setw( 2 ) << step << " ) = " << std::setw( 11 ) << l1eoc << std::endl;
    }

    grid.globalRefine( refineStepsForHalf );
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
