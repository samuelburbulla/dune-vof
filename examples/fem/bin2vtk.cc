#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/io/streams/binarystreams.hh>

// dune-vof includes
#include <dune/vof/femdfwrapper.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

// local includes
#include "polygon.hh"
#include "reconstructionwriter.hh"


// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int level )
  : level_( level )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : level_( other.level_ )
  {}

  bool willWrite ( ) { return true; }

  std::string prefix () const
  {
    std::stringstream s;
    s << "vof-fem-" << level_ << "-";
    return s.str();
  }

private:
  int level_;
};


// Write Binary data file to VTK file
// ----------------------------------

int main( int argc, char** argv )
try {
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  // check if input file is valid
  if ( argc != 2 )
  {
    std::cout << "No input file given." << std::endl;
    return 1;
  }

  auto file = std::string( argv[1] );

  std::ifstream f;
  f.open( file );
  if ( !f.good() )
  {
    std::cout << "Input file not found." << std::endl;
    return 1;
  }


  // Parse filename
  // --------------

  auto pfile = file.substr( 0, file.find_last_of('.') );
  const int number = std::stoi( pfile.substr( pfile.find_last_of('-')+1 ) );
  auto pfile2 = pfile.substr( 0, pfile.find_last_of('-') );
  const int level = std::stoi( pfile2.substr( pfile2.find_last_of('-')+1 ) );
  auto pfile3 = pfile2.substr( 0, pfile.find_last_of('-') );
  pfile3 = pfile3.substr( 0, pfile.find_last_of('-') );
  pfile3 = pfile3.substr( 0, pfile.find_last_of('-') );
  //const int rank = std::stoi( pfile3.substr( pfile3.find_last_of('p')+1 ) );
  auto pfile4 = pfile3.substr( 0, pfile.find_last_of('-') );
  const int size = std::stoi( pfile4.substr( pfile4.find_last_of('s')+1 ) );


  // supress output
  // --------------
  std::cout.setstate( std::ios_base::failbit );
  std::cerr.setstate( std::ios_base::failbit );


  // Create grid
  // -----------
  using GridType = Dune::GridSelector::GridType;

  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

  // GridType *gridPtr = Dune::BackupRestoreFacility< GridType >::restore( gridfile );
  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;

  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );


  // Create discrete function
  // ------------------------
  using GridPartType =
    Dune::Fem::LeafGridPart< GridType >;
  using FunctionSpaceType =
    Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;
  using DiscreteFunctionSpaceType =
    Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType =
    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;
  using TimeProviderType =
    Dune::Fem::FixedStepTimeProvider< typename GridType::CollectiveCommunication >;

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType space( gridPart );
  DiscreteFunctionType uh( "uh", space );


  // Read data from binary file
  std::stringstream binfile;
  binfile.fill('0');
  binfile << "./data/" << "s" << std::setw(4) << size << "-p" << std::setw(4) << Dune::Fem::MPIManager::rank() << "-vof-fem-"
    << std::setw(1) << level << "-" << std::setw(5) << number << ".bin";

  using BinaryStream = Dune::Fem::BinaryFileInStream;
  BinaryStream binaryStream ( binfile.str() );
  uh.read( binaryStream );




  // Rebuild flags and reconstruction
  // --------------------------------

  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridPartType >;
  ReconstructionSet reconstructions( gridPart );

  using Stencils = Dune::VoF::EdgeNeighborsStencil< GridPartType >;
  Stencils stencils( gridPart );

  const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );
  auto cuh = Dune::VoF::discreteFunctionWrapper( uh );
  auto reconstruction = Dune::VoF::reconstruction( gridPart, cuh, stencils );
  auto flags = Dune::VoF::flags( gridPart );
  DiscreteFunctionType dfFlags ( "flags", space );

  flags.reflag( cuh, eps );
  reconstruction( cuh, reconstructions, flags );


  for ( auto& entity : elements( gridPart ) )
    dfFlags.localFunction( entity )[0] = static_cast< double > ( flags[ entity ] );

  // Write data to vtk file
  // ----------------------
  using DataIOTupleType = std::tuple< DiscreteFunctionType *, DiscreteFunctionType * >;
  using DataOutputType = Dune::Fem::DataOutput< GridType, DataIOTupleType >;

  DataIOTupleType dataIOTuple = std::make_tuple( &uh, &dfFlags );
  DataOutputType dataOutput( grid, dataIOTuple );

  double savestep = Dune::Fem::Parameter::getValue< double >( "fem.io.savestep", 0 );
  savestep *= 1.0001;
  TimeProviderType tp ( 0.0, savestep );
  for ( int i = 0; i < number; ++i ) tp.next();

  std::stringstream outfile;
  outfile << pfile << ".vtk";
  dataOutput.write( tp, outfile.str() );


  // Write reconstruction to vtu file
  // --------------------------------
  using Polygon = Polygon< typename ReconstructionSet::Reconstruction::Coordinate >;
  using RecOutputType = ReconstructionWriter< ReconstructionSet, Polygon >;
  RecOutputType recOutput ( level );

  recOutput.count() = number;
  recOutput.write ( reconstructions );


  return 0;
}
catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
}

