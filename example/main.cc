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
#include <dune/grid/common/mcmgmapper.hh>

//- dune-vof includes
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/vertexneighborsstencil.hh>

//- local includes
#include "average.hh"
#include "colorFunction.hh"
#include "errors.hh"
#include "io.hh"
#include "polygon.hh"
#include "problem.hh"
#include "vtu.hh"




template < class GridView >
std::tuple< double, double > algorithm ( const GridView& gridView, const Dune::ParameterTree &parameters )
{
  // build domain references for each cell
  Dune::VoF::VertexNeighborsStencil< GridView > stencil( gridView );

  // allocate and initialize objects for data representation
  ColorFunction< GridView > colorFunction( gridView );
  ColorFunction< GridView > update( gridView );
  Dune::VoF::ReconstructionSet< GridView > reconstructionSet ( gridView );
  Dune::VoF::Flags< GridView, Dune::VoF::VertexNeighborsStencil< GridView > > flags ( gridView, stencil );

  // VTK Writer
  int numCells = parameters.get< int >( "grid.numCells" );
  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( numCells );
  createDirectory( path.str() );

  Dune::VTKSequenceWriter< GridView > vtkwriter ( gridView, "vof", path.str(), "~/dune" );
  vtkwriter.addCellData ( colorFunction, "celldata" );
  vtkwriter.addCellData ( flags, "flags" );

  std::vector< Polygon > recIO;
  VTUWriter< std::vector< Polygon > > vtuwriter( recIO );
  std::stringstream name;
  name.fill('0');
  name << "vof-rec-" << std::setw(5) << 0 << ".vtu";


  // calculate dt
  double dt = parameters.get< double >( "scheme.cflFactor" )  * ( 1.0 / numCells ) / psiMax();
  const double endTime = parameters.get< double >( "scheme.end", 10 );
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );

  int saveNumber = 1;
  const double saveInterval = std::max( parameters.get< double >( "io.saveInterval", 1 ), dt );
  double nextSaveTime = saveInterval;

  // Initial reconstruction
  average( colorFunction, [ ] ( const auto &x ) { return f( x, 0.0 ); } );

  flags.reflag( colorFunction, eps );
  Dune::VoF::reconstruct( gridView, colorFunction, reconstructionSet, domain, flags, eps );
  filterReconstruction( reconstructionSet, recIO );

  vtkwriter.write( 0 );
  vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );

  double t = 0;
  int k = 0;

  auto psit = [ &t ] ( const auto &x ) { return psi( x, t ); };

  while ( t < endTime )
  {
    ++k;

    flags.reflag( colorFunction, eps );

    Dune::VoF::reconstruct( gridView, colorFunction, reconstructionSet, domain, flags, eps );
    Dune::VoF::evolve( gridView, colorFunction, reconstructionSet, t, dt, flags, psit, eps, update );
    colorFunction.axpy( 1.0, update );

    t += dt;

    if ( 2.0 * std::abs( t - nextSaveTime ) < dt )
    {
      vtkwriter.write( t );

      filterReconstruction( reconstructionSet, recIO );
      std::stringstream name_;
      name_.fill('0');
      name_ << "vof-rec-" << std::setw(5) << saveNumber << ".vtu" ;
      vtuwriter.write( Dune::concatPaths( path.str(), name_.str() ) );


      nextSaveTime += saveInterval;
      ++saveNumber;
    }


    if ( (int)(( t / endTime ) * 100) > (int)(( (t-1) / endTime ) * 100) )
    	std::cerr << "\r" << numCells << " [" << (int)(( t / endTime ) * 100) << "%]";

    //std::cerr << "s=" << grid.size(0) << " k=" << k << " t=" << t << " dt=" << dt << " saved=" << saveNumber-1 << std::endl
  }

  auto ft = [ t ] ( const auto &x ) { return f( x, t ); };

  return std::make_tuple( l1error( colorFunction, ft ), l2error( colorFunction, ft ) );
}

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  std::cerr << std::endl << "### Parameters ###" << std::endl;
  parameters.report( std::cerr );
  std::cerr << "#################" << std::endl << std::endl;

  std::tuple<double, double> lastErrorTuple;

  // open errors file
  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" );
  createDirectory( path.str() );
  path << "/eoc.txt";
  std::ofstream eocFile( path.str(), std::ios_base::out );

  //  create grid
  std::stringstream gridFile;
  gridFile << "cube2d.dgf";

  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;

  int numCells = parameters.get< int >( "grid.numCells" );
  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  grid.globalRefine( refineStepsForHalf * log( numCells ) / log( 2 ) );

  int numRuns = parameters.get<int>( "runs", 1 );

  for ( int i = 0; i < numRuns; ++i )
  {
    int numCells = parameters.get< int >( "grid.numCells" );

    // start time integration
    auto errorTuple = algorithm( grid.leafGridView(), parameters );

      // print errors and eoc
    if ( i > 0 )
    {
      const double eocL1 = log( std::get< 0 >( lastErrorTuple )  / std::get< 0 >( errorTuple ) ) / M_LN2;
      const double eocL2 = log( std::get< 1 >( lastErrorTuple )  / std::get< 1 >( errorTuple ) ) / M_LN2;

      eocFile << "             " << eocL1 << "        " << eocL2 << std::endl;
    }

    eocFile << numCells << "    " << std::get<0> ( errorTuple ) << "       " << std::get<1> ( errorTuple ) << std::endl;

    lastErrorTuple = errorTuple;

    // refine
    parameters[ "grid.numCells" ] = std::to_string( numCells * 2 );
    std::cout << std::endl;

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
