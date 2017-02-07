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
  double end = 0.0;
  int numRuns = 10;
  int level = 0;

  //  create grid
  Dune::GridPtr< GridType > gridPtr( std::to_string( GridType::dimension ) + "dgrid.dgf" );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;
  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();

  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  grid.globalRefine( refineStepsForHalf * level );

  // Testproblem
  using ProblemType = PROBLEM < double, GridView::dimensionworld >;
  ProblemType problem;

  double lastL1Error = 0.0;
  std::ofstream errorsFile;
  errorsFile.open( "errors" );
  errorsFile << "#dx \terror " << std::endl;

  std::ofstream eocFile;
  eocFile.open( "eoc" );

  for ( int i = level; i < level + numRuns; ++i )
  {
    ColorFunction< GridView > uh( gridView );
    Dune::VoF::averageRecursive( uh, problem, start, i-level );
    uh.communicate();

    using DataWriter = Dune::VTKSequenceWriter< GridView >;
    const std::string path = "./" + parameters.get< std::string >( "io.path" ) + "/vof-" + std::to_string( i );
    createDirectory( path );

    DataWriter vtkwriter ( gridView, "vof", path, "" );
    vtkwriter.addCellData ( uh, "celldata" );

    // start time integration
    Dune::VoF::Algorithm< GridType::LeafGridView, ProblemType, DataWriter > algorithm( gridView, problem, vtkwriter, cfl, eps );
    vtkwriter.addCellData( algorithm.flags(), "flags" );
    double realEnd = end;
    algorithm( uh, start, realEnd );

    double partL1Error = Dune::VoF::l1error( gridView, algorithm.reconstructions(), algorithm.flags(), problem, realEnd, i-level );
    //double partL1Error = Dune::VoF::cellwiseL1error( uh, problem, realEnd, i-level );
    double L1Error = grid.comm().sum( partL1Error );

    // print errors and eoc
    if ( grid.comm().rank() == 0 )
    {
      const double eoc = log( lastL1Error / L1Error ) / M_LN2;

      errorsFile << 1.0 / 8.0 * std::pow( 2, -i ) << " \t" << L1Error << std::endl;
      eocFile << std::setprecision(0) << "    $" << 8 * std::pow( 2, i ) << "^2$ \t& " << std::scientific << std::setprecision(2) << L1Error << " & " << std::fixed << eoc << " \\\\" << std::endl;
      std::cout << "L1-Error( " << i << " ) =\t" << L1Error << std::endl;

      if ( i > level )
        std::cout << "EOC " << i << ": " << eoc << std::endl;
    }

    lastL1Error = L1Error;

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
