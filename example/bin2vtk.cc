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
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune-fem includes
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/streams/binarystreams.hh>

// dune-vof includes
#include <dune/vof/curvatureSet.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/curvature/generalheightfunctioncurvature.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flags.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

#include "../dune/vof/test/colorfunction.hh"
#include "../dune/vof/test/polygon.hh"
#include "../dune/vof/test/reconstructionwriter.hh"

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
    std::string prefix = Dune::Fem::Parameter::getValue< std::string >( "fem.io.prefix", "vof-data" );
    std::stringstream s;
    s << prefix << "-" << level_ << "-";
    return s.str();
  }

  virtual std::string path() const
  {
    return Dune::Fem::Parameter::getValue< std::string >( "fem.io.path", "data" );
  }

private:
  int level_, savecount_;
};


// RecOutputParameters
// -------------------

struct RecOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, RecOutputParameters >
{
  RecOutputParameters ( const int level, const int savecount )
   : level_ ( level ), savecount_ ( savecount ) {}

  virtual bool willWrite ( bool write ) const { return true; }

  virtual int startcounter () const { return savecount_; }

  virtual std::string prefix () const
  {
    std::string prefix = Dune::Fem::Parameter::getValue< std::string >( "fem.io.recprefix", "vof-rec" );
    std::stringstream s;
    s << prefix << "-" << level_ << "-";
    return s.str();
  }

  virtual std::string path() const
  {
    return Dune::Fem::Parameter::getValue< std::string >( "fem.io.path", "data" );
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
    pvtu << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-" << prefix << std::setw(5) << number << ".pvtu";
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
  std::size_t level ( Dune::Fem::Parameter::getValue< std::size_t >( "level", 0 ) );
  std::size_t repeats ( Dune::Fem::Parameter::getValue< std::size_t >( "repeats", 0 ) );
  const std::string path = Dune::Fem::Parameter::getValue< std::string >( "fem.io.path", "data" );


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
    seriesName << "./" << path << "/" << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-vof-fem-" << level;

    PVDWriter dataPVDWriter ( seriesName.str() + "-data.pvd" );
    PVDWriter recPVDWriter ( seriesName.str() + "-reconstruction.pvd" );

    // Create grid
    // -----------
    using GridType = Dune::GridSelector::GridType;

    std::stringstream gridFile;
    gridFile << "../dune/vof/test/" << GridType::dimension << "dgrid.dgf";

    Dune::GridPtr< GridType > gridPtr( gridFile.str() );
    gridPtr->loadBalance();
    GridType& grid = *gridPtr;

    const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
    grid.globalRefine( level * refineStepsForHalf );


    // Create discrete function
    // ------------------------
    using GridView = typename GridType::LeafGridView;
    using ColorFunction = ColorFunction< GridView >;

    GridView gridView( grid.leafGridView() );
    ColorFunction uh( gridView );

    using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
    ReconstructionSet reconstructions( gridView );

    using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
    Stencils stencils( gridView );

    const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );
    auto reconstruction = Dune::VoF::reconstruction( gridView, uh, stencils );

    using Flags = Dune::VoF::Flags< GridView >;
    Flags flags = Dune::VoF::flags( gridView );
    ColorFunction dfFlags ( gridView );

    using CurvatureOperator = Dune::VoF::GeneralHeightFunctionCurvature< GridView, Stencils, decltype( uh ), ReconstructionSet, Flags >;
    CurvatureOperator curvatureOperator ( gridView, stencils );
    using CurvatureSet = Dune::VoF::CurvatureSet< GridView >;
    CurvatureSet curvatureSet( gridView );

    using DataWriter = Dune::VTKWriter< GridView >;

    DataWriter vtkwriter ( gridView );
    vtkwriter.addCellData ( uh, "celldata" );
    vtkwriter.addCellData ( flags, "flags" );
    vtkwriter.addCellData ( curvatureSet, "curvature" );



    for ( std::size_t number = 0; ; number++ )
    {
      DataOutputParameters dataOutputParameters ( level, number );
      RecOutputParameters recOutputParameters ( level, number );

      // Open next binary file
      // ---------------------
      std::stringstream filename;
      filename.fill('0');
      filename  << "./" << path << "/"  << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-p" << std::setw(4) << Dune::Fem::MPIManager::rank()
        << "-" << dataOutputParameters.prefix() << std::setw(5) << number;

      if ( !std::ifstream ( filename.str() + ".bin" ) ) break;

      using BinaryStream = Dune::Fem::BinaryFileInStream;
      BinaryStream binaryStream ( filename.str() + ".bin" );

      double timeValue;
      binaryStream >> timeValue;

      uh.read( binaryStream );
      uh.communicate();

      // Rebuild flags and reconstruction
      // --------------------------------
      flags.reflag( uh, eps );
      reconstruction( uh, reconstructions, flags );
      curvatureOperator( reconstructions, flags, curvatureSet );

      // Write data to vtk file
      // ----------------------
      std::stringstream vtkfile;
      vtkfile.fill('0');
      vtkfile << dataOutputParameters.prefix() << std::setw(5) << number;
      vtkwriter.pwrite( vtkfile.str(), dataOutputParameters.path(), "" );
      dataPVDWriter.addDataSet( dataOutputParameters.prefix(), number, timeValue );

      // Write reconstruction to vtu file
      // --------------------------------
      using Polygon = OutputPolygon< typename ReconstructionSet::DataType::Coordinate >;
      using RecOutputType = ReconstructionWriter< GridView, ReconstructionSet, Flags, RecOutputParameters, Polygon >;

      RecOutputType recOutput ( gridView, reconstructions, flags, recOutputParameters );
      recOutput.write();

      recPVDWriter.addDataSet( recOutputParameters.prefix(), number, timeValue );
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

