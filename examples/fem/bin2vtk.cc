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
  DataOutputParameters ( const int level, const int savecount )
   : level_ ( level ), savecount_ ( savecount ) {}

  virtual bool willWrite ( bool write ) const { return true; }

  virtual int startcounter () const { return savecount_; }

  virtual std::string prefix () const
  {
    std::stringstream s;
    s << "vof-fem-" << level_ << "-";
    return s.str();
  }

private:
  int level_, savecount_;
};



struct PVDWriter
{
  PVDWriter ( const std::string& filename ) : filename_ ( filename )
  {
    content_ << "<?xml version=\"1.0\"?>" << std::endl << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl << "  <Collection>" << std::endl;
  }

  ~PVDWriter ()
  {
    content_ << "  </Collection>" << std::endl << "</VTKFile>" << std::endl;
    std::ofstream file ( filename_ );
    file << content_.str();
    file.close();
  }

  void addDataSet ( const std::string& prefix, const std::size_t number, const double timeValue )
  {
    std::stringstream pvtu;
    pvtu.fill('0');
    pvtu << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-" << prefix << std::setw(6) << number << ".pvtu";
    content_ << "    <DataSet timestep=\"" << timeValue << "\" group=\"\" part=\"0\" file=\"" << pvtu.str() << "\"/>" << std::endl;
  }

private:
  std::string filename_;
  std::stringstream content_;
};




// Write Binary data file to VTK file
// ----------------------------------

int main( int argc, char** argv )
try {
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  // Create filename
  // ---------------
  const std::string name ( "vof-fem" );
  std::size_t level ( Dune::Fem::Parameter::getValue( "level", 0 ) );
  std::size_t repeats ( Dune::Fem::Parameter::getValue( "repeats", 0 ) );

  if ( argc > 1 )
  {
    level = std::stoi( argv[ 1 ] );
    repeats = 0;
  }

  const std::size_t maxLevel = level + repeats;
  for ( ; level <= maxLevel; ++level )
  {
    std::stringstream seriesName;
    seriesName.fill('0');
    seriesName << "./data/" << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-vof-fem-" << level;

    PVDWriter dataPVDWriter ( seriesName.str() + "-data.pvd" );
    PVDWriter recPVDWriter ( seriesName.str() + "-reconstruction.pvd" );

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

    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType space( gridPart );
    DiscreteFunctionType uh( "uh", space );


    for ( std::size_t number = 0; ; number++ )
    {
      // Open next binary file
      // ---------------------
      std::stringstream filename;
      filename.fill('0');
      filename << "./data/" << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-p" << std::setw(4) << Dune::Fem::MPIManager::rank()
        << "-" << name << "-" << level << "-" << std::setw(5) << number;

      if ( !std::ifstream ( filename.str() + ".bin" ) ) break;

      using BinaryStream = Dune::Fem::BinaryFileInStream;
      BinaryStream binaryStream ( filename.str() + ".bin" );

      double timeValue;
      binaryStream >> timeValue;

      uh.read( binaryStream );
      uh.communicate();



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
      DataOutputParameters dataOutputParameters ( DataOutputParameters( level, number ) );
      DataOutputType dataOutput( grid, dataIOTuple, dataOutputParameters );
      dataOutput.write();

      dataPVDWriter.addDataSet( dataOutputParameters.prefix(), number, timeValue );

      // Write reconstruction to vtu file
      // --------------------------------
      using Polygon = OutputPolygon< typename ReconstructionSet::Reconstruction::Coordinate >;
      using RecOutputType = ReconstructionWriter< GridPartType, ReconstructionSet, Polygon >;
      RecOutputType recOutput ( gridPart, level );

      recOutput.count() = number;
      recOutput.write ( reconstructions );

      recPVDWriter.addDataSet( "vof-rec-" + std::to_string( level ) + "-", number, timeValue );
    }

  }

  return 0;
}
catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
}

