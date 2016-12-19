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
#include <dune/common/parametertreeparser.hh>

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
double algorithm ( Grid &grid, ColorFunction& uh, P& problem, const Dune::ParameterTree &parameters, int level, double start, double end, double cfl, double eps, const bool verbose = false, const bool writeData = true )
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
  DataOutputType dataOutput( parameters, level, time );

  // Calculate and write initial data
  if ( writeData )
    dataOutput.write( grid, uh, time );

  ColorFunction update( gridView );

  // Time Iteration
  for( ; time <= end; )
  {
    if ( verbose )
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

  Dune::MPIHelper::instance( argc, argv );

  // Read Parameter File
  // ===================
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  int level = parameters.get< int >( "grid.level", 0 );
  int repeats = parameters.get< int >( "grid.repeats", 0 );
  double start = parameters.get< double >( "scheme.start", 0.0 );
  double end = parameters.get< double >( "scheme.end", 2.5 );
  double cfl = parameters.get< double >( "scheme.cfl", 1.0 );
  double eps = parameters.get< double >( "scheme.eps", 1e-9 );
  std::string path = parameters.get< std::string >( "fem.io.path", "data" );
  int restartStep = parameters.get< int >( "io.restartStep", -1 );
  int verboserank = parameters.get< int >( "io.verboserank", -1 );
  bool writeData = parameters.get< bool >( "io.writeData", true );

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
      namedata << "s" << std::setw(4) << grid.comm().size() << "-p" << std::setw(4) << grid.comm().rank()
        << "-vof-" << std::to_string( step ) << "-" << std::setw(5) << restartStep << ".bin";
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
        start = timestamp;
        uh.read( binaryStream );
        std::cout << "Restarted in file " << filename << std::endl;
      }
    }

    uh.communicate();

    bool verbose = ( verboserank == grid.comm().rank() );

    // Run Algorithm
    double error = algorithm( grid, uh, problem, parameters, step, start, end, cfl, eps, verbose, writeData );


    if ( grid.comm().rank() == 0 )
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
