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

// dune-vof includes
#include <dune/vof/colorfunction.hh>
#include <dune/vof/eoc.hh>
#include "../dune/vof/test/average.hh"
#include "../dune/vof/test/problems/linearwall.hh"
#include "../dune/vof/test/problems/rotatingcircle.hh"
#include "../dune/vof/test/problems/sflow.hh"
#include "../dune/vof/test/problems/slottedcylinder.hh"
#include <dune/vof/algorithm.hh>

#include "binarywriter.hh"

int main(int argc, char** argv)
try {

  Dune::MPIHelper::instance( argc, argv );

  // Read Parameter File
  // ===================
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  int level0 = parameters.get< int >( "grid.level", 0 );
  int repeats = parameters.get< int >( "grid.repeats", 0 );
  double start = parameters.get< double >( "scheme.start", 0.0 );
  double end = parameters.get< double >( "scheme.end", 2.5 );
  double cfl = parameters.get< double >( "scheme.cfl", 1.0 );
  double eps = parameters.get< double >( "scheme.eps", 1e-9 );
  std::string path = parameters.get< std::string >( "io.path", "data" );
  int restartStep = parameters.get< int >( "io.restartStep", -1 );
  int verboserank = parameters.get< int >( "io.verboserank", -1 );

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
  grid.globalRefine( level0 * refineStepsForHalf );

  GridView gridView( grid.leafGridView() );

  using ColorFunction = Dune::VoF::ColorFunction< GridView >;

  // Testproblem
  using ProblemType = RotatingCircle< double, GridView::dimensionworld >;
  ProblemType problem;

  const double dx0 = 1.0 / 8.0 * std::pow( 2.0, -level0 );
  Dune::EocOutput eocOutput ( "errors", dx0, level0, 2 );

  // EOC Calculation
  // ===============
  for( int level = level0; level <= level0+repeats; ++level )
  {
    // Initialize Data
    // ===============
    ColorFunction uh( gridView );


    if ( restartStep == -1 )
    {
      // Use initial data of problem.
      Dune::VoF::Average< ProblemType > average ( problem );
      average( uh, start, ( gridView.comm().rank() == 0 && level == level0 ) );
    }
    else
    {
      // Use given data in binary file.
      std::stringstream namedata;
      namedata.fill('0');
      namedata << "s" << std::setw(4) << grid.comm().size() << "-p" << std::setw(4) << grid.comm().rank()
        << "-vof-" << std::to_string( level ) << "-" << std::setw(5) << restartStep << ".bin";
      const auto filename = Dune::concatPaths( path, namedata.str() );
      if ( !Dune::Fem::fileExists( filename ) )
      {
        std::cout << "Restart error: A time step is given to load binary file but this file does not exist." << std::endl;
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
      uh.communicate();
    }


    bool verbose = ( verboserank == grid.comm().rank() );

    using DataOutputType = BinaryWriter< GridView, ColorFunction >;
    DataOutputType dataOutput( gridView, uh, parameters, level );

    // Run Algorithm
    Dune::VoF::Algorithm< GridView, ProblemType, DataOutputType > algorithm( gridView, problem, dataOutput, cfl, eps, verbose );

    double partError = algorithm( uh, start, end, level );
    double error = grid.comm().sum( partError );

    if ( grid.comm().rank() == 0 )
      eocOutput.add( error );

    grid.globalRefine( refineStepsForHalf );
  }

  return 0;
}
catch ( Dune::Exception &e ) {
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...) {
  std::cerr << "Unknown exception thrown!" << std::endl;
}
