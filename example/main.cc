#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <math.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh> 
#include <dune/common/path.hh>
#include <fstream>
#include <iostream>
#include <vector>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/common/parametertreeparser.hh>


#include "colorFunction.hh"
#include "domain.hh"
#include "testproblem.hh"

#include "src/flags.hh"
#include "src/errors.hh"
#include "src/io.hh"
#include "src/L1projection.hh"
#include "src/evolve.hh"
#include "src/reconstruct.hh"
#include "src/reconstructionSet.hh"
#include "src/geometricaltoolbox.hh"
#include "src/cputime.hh"
#include "src/polygon.hh"
#include "src/vtu.hh"




template < class GridView >
std::tuple< double, double > algorithm ( const GridView& gridView, const Dune::ParameterTree &parameters )
{
  const int dimworld = GridView::dimensionworld;
  typedef typename GridView::ctype ctype;
  typedef typename Dune::FieldVector< ctype, dimworld > fvector;
  typedf typename Polygon< fvector > Polygon;

  // build domain references for each cell
  Dune::VoF::DomainOfPointNeighbors< GridView > domain ( gridView ); 

  // allocate and initialize objects for data representation
  Dune::VoF::ColorFunction< GridView > colorFunction( gridView );
  Dune::VoF::ColorFunction< GridView > update( gridView );
  Dune::VoF::ReconstructionSet< GridView > reconstructionSet ( gridView );
  Dune::VoF::Flags< GridView, Dune::VoF::DomainOfPointNeighbors<GridView> > flags ( gridView, domain );



  // VTK Writer
  int numCells = parameters.get< int >( "grid.numCells" );
  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( numCells );
  Dune::VoF::createDirectory( path.str() );


  Dune::VTKSequenceWriter< GridView > vtkwriter ( gridView, "vof", path.str(), "~/dune" );
  vtkwriter.addCellData ( colorFunction, "celldata" );
  vtkwriter.addCellData ( flags, "flags" );


  std::vector< Polygon > recIO;
  VTUWriter< std::vector< Polygon > > vtuwriter( recIO );
  std::stringstream name;
  name.fill('0');
  name << "vof-rec-" << std::setw(5) << 0 << ".vtu";


  int saveNumber = 1;
  const double saveInterval = parameters.get< double >( "io.saveInterval", 1 );
  double nextSaveTime = saveInterval;

  // calculate dt
  double dt = std::min( parameters.get< double >( "scheme.cflFactor" )  * ( 1.0 / numCells ) / Dune::VoF::psiMax(), saveInterval );
  const double endTime = parameters.get< double >( "scheme.end", 10 );
  
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );




  // Initial reconstruction
  Dune::VoF::L1projection( colorFunction, Dune::VoF::f0<fvector> );

  flags.reflag( colorFunction, eps );
  Dune::VoF::reconstruct( gridView, colorFunction, reconstructionSet, domain, flags, eps );
  Dune::VoF::filterReconstruction( reconstructionSet, recIO );

  vtkwriter.write( 0 );
  vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );


  double t = 0;
  int k = 0;

  auto psit = [t] ( fvector x ) -> fvector { return Dune::VoF::psi( x, t ); };

  while ( t < endTime )
  {
    ++k;

    flags.reflag( colorFunction, eps );

    Dune::VoF::reconstruct( gridView, colorFunction, reconstructionSet, domain, flags, eps );

    Dune::VoF::evolve( gridView, colorFunction, reconstructionSet, domain, t, dt, flags, psit, eps, update );
    colorFunction.axpy( 1.0, update );

    t += dt;

    if ( std::abs( t - nextSaveTime ) <= saveInterval / 2.0 )
    {

      vtkwriter.write( t );

      Dune::VoF::filterReconstruction( reconstructionSet, recIO );
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

  auto ft = [t] ( fvector x ) -> ctype { return Dune::VoF::f( x, t ); };

  double L1 = Dune::VoF::l1error( gridView, colorFunction, ft );
  double L2 = Dune::VoF::l2error( gridView, colorFunction, ft );

  return std::tuple<double, double> ( L1, L2 );

  return std::tuple<double, double> ( 0 , 0 );
}

int main(int argc, char** argv)
{

    Dune::MPIHelper::instance( argc, argv );

    try {

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
    std::stringstream errorsPath;
    errorsPath << "./" << parameters.get< std::string >( "io.folderPath" ) << "/errors";
    std::fstream errorFile;
    errorFile.open( errorsPath.str(), std::fstream::out );

    errorFile << "Cells   L1    eoc     L2     eoc    " << std::endl << std::endl;


    //  create grid
    std::stringstream gridFile;
    gridFile << "cube2.dgf";

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

        errorFile << "             " << eocL1 << "        " << eocL2 << std::endl;
      }

      errorFile << numCells << "    " << std::get<0> ( errorTuple ) << "       " << std::get<1> ( errorTuple ) << std::endl;

      lastErrorTuple = errorTuple;

      // refine
      parameters[ "grid.numCells" ] = std::to_string( numCells * 2 );
      std::cout << std::endl;
    
      grid.globalRefine( refineStepsForHalf );
    }

    errorFile.close();

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
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }

}
