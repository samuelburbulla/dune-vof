#include "config.h"

//- C++ includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <vector>

//- dune-common includes
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/path.hh>
#include <dune/common/parametertreeparser.hh>

//- dune-grid includes
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

//- dune-vof includes
#include <dune/vof/algorithm.hh>
#include <dune/vof/colorfunction.hh>

//- local includes
#include "average.hh"
#include "errors.hh"
#include "io.hh"
#include "problems/rotatingcircle.hh"
#include "problems/linearwall.hh"
#include "problems/slope.hh"

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  double cfl = 0.25;
  double eps = 1e-6;
  double start = 0.0;
  double end = 1.0;
  int numRuns = 2;
  int level = 0;

  //  create grid
  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;
  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();

  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  grid.globalRefine( refineStepsForHalf * level );

  // Testproblem
  using ProblemType = PROBLEM < double, GridView::dimensionworld >;
  ProblemType problem;

  double lastL1Error;
  std::ofstream errorsFile;
  errorsFile.open( "errors" );
  errorsFile << "#dx \terror " << std::endl;

  for ( int i = 0; i < numRuns; ++i )
  {
    ColorFunction< GridView > uh( gridView );
    Dune::VoF::average( uh, problem );
    uh.communicate();

    using DataWriter = Dune::VTKSequenceWriter< GridView >;
    std::stringstream path;
    path << "./" << parameters.get< std::string >( "io.path" ) << "/vof-" << std::to_string( level );
    createDirectory( path.str() );

    DataWriter vtkwriter ( gridView, "vof", path.str(), "" );
    vtkwriter.addCellData ( uh, "celldata" );

    // start time integration
    Dune::VoF::Algorithm< GridType::LeafGridView, ProblemType, DataWriter > algorithm( gridView, problem, vtkwriter, cfl, eps );
    vtkwriter.addCellData( algorithm.flags(), "flags" );
    algorithm( uh, start, end );

    auto partL1Error = Dune::VoF::l1error( gridView, algorithm.reconstructions(), algorithm.flags(), problem, end );
    double L1Error = grid.comm().sum( partL1Error );

    errorsFile << 1.0 / 8.0 * std::pow( 2, -level ) << " \t" << L1Error << std::endl;

    // print errors and eoc
    if ( i > 0 && grid.comm().rank() == 0 )
    {
      const double eoc = log( lastL1Error / L1Error ) / M_LN2;

      std::cout << "EOC " << i << ": " << eoc << std::endl;
    }

    lastL1Error = L1Error;

    // refine
    parameters[ "grid.level" ] = std::to_string( ++level );
    grid.globalRefine( refineStepsForHalf );
  }

  return 0;
}
catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch ( std::exception &e )
{
  std::cerr << "STD::EXCEPTION THROWN: \"" << e.what() << "\"" << std::endl;
}
catch (...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
