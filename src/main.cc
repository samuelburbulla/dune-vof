#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/path.hh>
#include <fstream>
#include <iostream>
#include <vector>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/yaspgrid.hh>
#include <cmath>
#include <functional>
#include <stdexcept>

#include <dune/common/parametertreeparser.hh>


#include <dune/vof/flagCells.hh>
#include <dune/vof/errors.hh>
#include <dune/vof/io.hh>
#include <dune/vof/domain.hh>
#include <dune/vof/initialize.hh>
#include <dune/vof/evolve.hh>
#include <dune/vof/reconstruct.hh>


#include "polygon.hh"
#include "vtu.hh"

using fvector =
 Dune::FieldVector<double,2>;
using polygon =
  Polygon< fvector >;

void filterReconstruction( const std::vector< std::array<fvector,3 > > &rec, std::vector< polygon > &io )
{
  io.clear();
  for( auto && it: rec)
    if( it[ 2 ] != fvector( 0.0 ) )
      io.push_back( polygon{ it[ 0 ], it[ 1 ] } );
}



template < class Grid >
std::tuple<double, double> algorithm ( const Grid& grid, const Dune::ParameterTree &parameters )
{
    // build domain references for each cell
  Dune::VoF::DomainOfCells<Grid> domain ( grid );

  auto gridView = grid.leafGridView();

  int n = grid.leafGridView().size(0); 

  // allocate and initialize vectors for data representation
  std::vector<double> concentration ( n );
  std::vector< std::array<fvector,3> > reconstruction( n );
  std::vector< polygon > recIO;
  std::vector<bool> cellIsMixed ( n );
  std::vector<bool> cellIsActive ( n );
  std::vector<fvector> velocityField ( n );
  std::vector<double> overundershoots ( n );


  Dune::VoF::initialize( grid, concentration, Dune::VoF::f0<fvector> );

  // calculate dt
  int numCells = parameters.get< int >( "grid.numCells" );
  double dt = parameters.get< double >( "scheme.cflFactor" )  * ( 1.0 / numCells ) / Dune::VoF::psiMax();
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );


  double t = 0;
  int k = 0;


  Dune::VoF::flagCells( gridView, concentration, reconstruction, domain, cellIsMixed, cellIsActive, eps );

  Dune::VoF::clearReconstruction( reconstruction );
  Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, eps );
  filterReconstruction( reconstruction, recIO );

  int saveNumber = 1;
  const double saveInterval = parameters.get< double >( "io.saveInterval", 1 );
  double nextSaveTime = saveInterval;


  // VTK Writer
  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( numCells );

  Dune::VoF::createDirectory( path.str() );

  Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkwriter ( gridView, "vof", path.str(), "~/dune" );

  vtkwriter.addCellData ( concentration, "celldata" );
  vtkwriter.addCellData ( cellIsMixed, "cellmixed" );
  vtkwriter.addCellData ( cellIsActive, "cellactive" );
  vtkwriter.addCellData ( overundershoots, "overundershoots" );

  VTUWriter< std::vector< polygon > > vtuwriter( recIO );

  vtkwriter.write( 0 );

  std::stringstream name;
  name.fill('0');
  name << "vof-rec-" << std::setw(5) << 0 << ".vtu";

  vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );


  auto psit = std::bind( Dune::VoF::psi<fvector>, std::placeholders::_1, 0);
  Dune::VoF::L1projection( grid, velocityField, psit );

  const double endTime = parameters.get< double >( "scheme.end", 10 );

  while ( t < endTime )
  {
    ++k;

    Dune::VoF::clearReconstruction( reconstruction );

    Dune::VoF::flagCells( grid.leafGridView(), concentration, reconstruction, domain, cellIsMixed, cellIsActive, eps );
    Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, eps );
    Dune::VoF::evolve( grid, concentration, reconstruction, domain, numCells, t, dt, eps, cellIsMixed, cellIsActive, velocityField, overundershoots );
    filterReconstruction( reconstruction, recIO );

    t += dt;

    if ( std::abs( t - nextSaveTime ) < saveInterval/4.0 )
    {
      vtkwriter.write( t );

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

  auto ft = std::bind( Dune::VoF::f<fvector>, std::placeholders::_1, t);

  double L1 = Dune::VoF::l1error( grid, concentration, ft );
  double L2 = Dune::VoF::l2error( grid, concentration, ft );

  return std::tuple<double, double> ( L1, L2 );
}

int main(int argc, char** argv)
{

    Dune::MPIHelper::instance( argc, argv );

    try {


      // type declarations
    const int dim = 2;
      typedef Dune::YaspGrid<dim> GridType;
    typedef typename Dune::FieldVector<double,dim> fvector;

    // set parameters
    Dune::ParameterTree parameters;
    Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
    Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

    std::cerr << std::endl << "### Parameters ###" << std::endl;
    parameters.report( std::cerr );
    std::cerr << "#################" << std::endl << std::endl;

    std::tuple<double, double> lastErrorTuple;

    std::stringstream errorsPath;
    errorsPath << "./" << parameters.get< std::string >( "io.folderPath" ) << "/errors";

    std::fstream errorFile;
    errorFile.open( errorsPath.str(), std::fstream::out );

    errorFile << "Cells   L1    eoc     L2     eoc    " << std::endl << std::endl;



    int numRuns = parameters.get<int>( "runs", 1 );

    for ( int i = 0; i < numRuns; ++i )
    {

      // build Grid
      fvector upper( 1.0 );
      int numCells = parameters.get< int >( "grid.numCells" );
      Dune::array<int,dim> noc;
      std::fill( noc.begin(), noc.end(), numCells );
      GridType grid( upper, noc );

      // start time integration
      auto errorTuple = algorithm( grid, parameters );

        // print errors and eoc
      if ( i > 0 )
      {
        const double eocL1 = log( std::get< 0 >( lastErrorTuple )  / std::get< 0 >( errorTuple ) ) / M_LN2;
        const double eocL2 = log( std::get< 1 >( lastErrorTuple )  / std::get< 1 >( errorTuple ) ) / M_LN2;

        errorFile << "             " << eocL1 << "        " << eocL2 << std::endl;

      }

      errorFile << numCells << "    " << std::get<0> ( errorTuple ) << "       " << std::get<1> ( errorTuple ) << std::endl;

      lastErrorTuple = errorTuple;

      // refine
      parameters[ "grid.numCells" ] = std::to_string( numCells * 2 );
      std::cerr << std::endl;
    }

    errorFile.close();

    return 0;
  }

  catch ( Dune::Exception &e )
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch( int e )
  {
    if ( e == 10 )
      std::cerr << "Error: No intersection in cell with his reconstruction." << std::endl;
    else
      throw;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }

}
