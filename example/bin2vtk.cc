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
#include <dune/common/parametertreeparser.hh>

// dune-fem includes
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/streams/binarystreams.hh>

// dune-vof includes
#include <dune/vof/colorfunction.hh>
#include <dune/vof/curvatureset.hh>
#include <dune/vof/reconstructionset.hh>
#include <dune/vof/curvature/cartesianheightfunctioncurvature.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flagging.hh>
#include <dune/vof/flagset.hh>
#include <dune/vof/interfacegrid/grid.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/stencil/edgeneighborsstencil.hh>

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int level, const int savecount )
   : level_ ( level ), savecount_ ( savecount )
  {
    Dune::ParameterTreeParser::readINITree( "parameter", parameters_ );
  }

  virtual bool willWrite ( bool write ) const { return true; }

  virtual int startcounter () const { return savecount_; }

  virtual std::string prefix () const
  {
    std::string prefix = parameters_.get< std::string >( "io.prefix", "vof-data" );
    std::stringstream s;
    s << prefix << "-" << level_ << "-";
    return s.str();
  }

  virtual std::string path() const
  {
    return parameters_.get< std::string >( "io.path", "data" );
  }

private:
  int level_, savecount_;
  Dune::ParameterTree parameters_;
};


// RecOutputParameters
// -------------------

struct RecOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, RecOutputParameters >
{
  RecOutputParameters ( const int level, const int savecount )
   : level_ ( level ), savecount_ ( savecount )
   {
     Dune::ParameterTreeParser::readINITree( "parameter", parameters_ );
   }

  virtual bool willWrite ( bool write ) const { return true; }

  virtual int startcounter () const { return savecount_; }

  virtual std::string prefix () const
  {
    std::string prefix = parameters_.get< std::string >( "io.recprefix", "vof-rec" );
    std::stringstream s;
    s << prefix << "-" << level_ << "-";
    return s.str();
  }

  virtual std::string path() const
  {
    return parameters_.get< std::string >( "io.path", "data" );
  }

private:
  int level_, savecount_;
  Dune::ParameterTree parameters_;

};



struct PVDWriter
{
  PVDWriter ( const std::string& filename, const std::string ending = "pvtu" ) : filename_ ( filename ), ending_( ending )
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

  void addDataSet ( const std::size_t size, const std::string& prefix, const std::size_t number, const double timeValue )
  {
    std::stringstream pvtu;
    pvtu.fill('0');
    pvtu << "s" << std::setw(4) << size<< "-" << prefix << std::setw(5) << number << "." << ending_;
    content_ << "    <DataSet timestep=\"" << timeValue << "\" group=\"\" part=\"0\" file=\"" << pvtu.str() << "\"/>" << std::endl;
  }

private:
  std::string filename_, ending_;
  std::stringstream content_;
};




// Write Binary data file to VTK file
// ----------------------------------

int main( int argc, char** argv )
try {
  Dune::MPIHelper::instance( argc, argv );

  // read parameter file
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  // Create filename
  // ---------------
  const std::string name ( "vof" );
  std::size_t level ( parameters.get< std::size_t >( "grid.level", 0 ) );
  std::size_t repeats ( parameters.get< std::size_t >( "grid.repeats", 0 ) );
  const std::string path = parameters.get< std::string >( "io.path", "data" );

  if ( argc > 1 )
  {
    level = std::stoi( argv[ 1 ] );
    repeats = 0;
  }

  const std::size_t maxLevel = level + repeats;
  for ( ; level <= maxLevel; ++level )
  {
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


    std::stringstream seriesName;
    seriesName.fill('0');
    seriesName << "./" << path << "/" << "s" << std::setw(4) << grid.comm().size() << "-vof-" << level;

    PVDWriter dataPVDWriter ( seriesName.str() + "-data.pvd", "pvtu" );
    PVDWriter recPVDWriter ( seriesName.str() + "-reconstruction.pvd", "pvtp" );

    // Create discrete function
    // ------------------------
    using GridView = typename GridType::LeafGridView;
    using ColorFunction = ColorFunction< GridView >;

    GridView gridView( grid.leafGridView() );
    ColorFunction uh( gridView );

    using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
    Stencils stencils( gridView );

    using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
    ReconstructionSet reconstructions( gridView );
    auto reconstruction = Dune::VoF::reconstruction( gridView, uh, stencils );

    const double eps = parameters.get< double >( "scheme.eps", 1e-9 );

    using FlagSet = Dune::VoF::FlagSet< GridView >;
    FlagSet flags( gridView );
    auto flagOperator = Dune::VoF::FlagOperator< ColorFunction, FlagSet >( eps );
    ColorFunction dfFlags ( gridView );

    using CurvatureOperator = Dune::VoF::CartesianHeightFunctionCurvature< GridView, Stencils, decltype( uh ), ReconstructionSet, FlagSet >;
    CurvatureOperator curvatureOperator ( gridView, stencils );
    using CurvatureSet = Dune::VoF::CurvatureSet< GridView >;
    CurvatureSet curvatureSet( gridView );

    Dune::VTKWriter< GridView > vtkwriter ( gridView );
    vtkwriter.addCellData ( uh, "celldata" );
    vtkwriter.addCellData ( dfFlags, "flags" );

    for ( std::size_t number = 0; ; number++ )
    {
      DataOutputParameters dataOutputParameters ( level, number );
      RecOutputParameters recOutputParameters ( level, number );

      // Open next binary file
      // ---------------------
      std::stringstream filename;
      filename.fill('0');
      filename  << "./" << path << "/"  << "s" << std::setw(4) << grid.comm().size() << "-p" << std::setw(4) << grid.comm().rank()
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
      flagOperator( uh, flags );
      reconstruction( uh, reconstructions, flags );
      curvatureOperator( uh, reconstructions, flags, curvatureSet );

      // Write data
      // ----------
      for ( const auto &entity : elements( gridView ) )
        dfFlags[ entity ] = static_cast< double > ( flags[ entity ] );

      std::stringstream vtkfile;
      vtkfile.fill('0');
      vtkfile << dataOutputParameters.prefix() << std::setw(5) << number;
      vtkwriter.pwrite( vtkfile.str(), dataOutputParameters.path(), "" );
      dataPVDWriter.addDataSet( grid.comm().size(), dataOutputParameters.prefix(), number, timeValue );

      // Write reconstruction
      // --------------------
      using InterfaceGrid = Dune::VoF::InterfaceGrid< decltype( reconstruction ) >;
      InterfaceGrid interfaceGrid( uh, stencils );

      std::stringstream recfile;
      recfile.fill('0');
      recfile << recOutputParameters.prefix() << std::setw(5) << number;


      Dune::VoF::DataSet< InterfaceGrid::LeafGridView, double > curvatureOnInterface ( interfaceGrid.leafGridView() );
      for ( const auto entity : elements( interfaceGrid.leafGridView() ) )
        curvatureOnInterface[ entity ] = curvatureSet[ entity.impl().hostElement() ];

      Dune::VTKWriter< InterfaceGrid::LeafGridView > interfaceVtkWriter( interfaceGrid.leafGridView() );
      interfaceVtkWriter.addCellData( curvatureOnInterface, "curvature" );
      interfaceVtkWriter.pwrite( recfile.str(), recOutputParameters.path(), "" );
      recPVDWriter.addDataSet( grid.comm().size(), recOutputParameters.prefix(), number, timeValue );
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

