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
#include <dune/vof/colorfunction.hh>
#include <dune/vof/curvatureset.hh>
#include <dune/vof/evolution.hh>
#include <dune/vof/flagset.hh>
#include <dune/vof/flagging.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionset.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/curvature/cartesianheightfunctioncurvature.hh>

//- local includes
#include "average.hh"
#include "errors.hh"
#include "io.hh"
#include "problems/ellipse.hh"
#include "problems/rotatingcircle.hh"
#include "problems/slope.hh"


// initTimeStep
// ---------------------

template < class GridView, class Velocity >
double initTimeStep( const GridView& gridView, const Velocity &velocity )
{
  double dtMin = std::numeric_limits< double >::max();

  for( const auto &entity : elements( gridView, Dune::Partitions::interiorBorder ) )
  {
    const auto geoEn = entity.geometry();
    for ( const auto &intersection : intersections( gridView, entity ) )
    {
      const auto geoIs = intersection.geometry();
      auto v = velocity( geoIs.center() );
      dtMin = std::min( dtMin, ( geoEn.volume() / geoIs.volume() ) / v.two_norm() );
    }
  }
  return dtMin;
}

// algorithm
// ---------

template < class GridView >
double algorithm ( const GridView& gridView, const Dune::ParameterTree &parameters )
{
  using DomainVector = Dune::FieldVector< double, GridView::dimensionworld >;
  using ColorFunction = ColorFunction< GridView >;
  using CurvatureSet = Dune::VoF::CurvatureSet< GridView >;
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
  using Flags = Dune::VoF::FlagSet< GridView >;
  using FlagOperator = Dune::VoF::FlagOperator< ColorFunction, Flags >;

  using DataWriter = Dune::VTKSequenceWriter< GridView >;

  // Testproblem
  using ProblemType = Ellipse< double, GridView::dimensionworld >;
  DomainVector axis ( { 1.0 / std::sqrt( 2 ), 1.0 / std::sqrt( 2 ) } );
  ProblemType problem ( { axis, Dune::VoF::generalizedCrossProduct( axis ) }, { 0.2, 0.5 } );
  //DomainVector axis2 ( { 0.0, 1.0, 0.0 } );
  //DomainVector axis3 ( { 0.0, 0.0, 1.0 } );
  //ProblemType problem ( { axis, axis2, axis3 }, { 0.4, 0.4, 0.4 } );
  //using ProblemType = Slope< double, GridView::dimensionworld >;
  //ProblemType problem ( 0.125 * M_PI );

  // calculate dt
  int level = parameters.get< int >( "grid.level" );

  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );
  const bool writeData = parameters.get< bool >( "io.writeData", 0 );

  // build domain references for each cell
  Stencils stencils( gridView );

  using CurvatureOperator = Dune::VoF::CartesianHeightFunctionCurvature< GridView, Stencils, ColorFunction, ReconstructionSet, Flags >;
  CurvatureOperator curvatureOperator ( gridView, stencils );
  CurvatureSet curvatureSet( gridView );
  ColorFunction curvatureError( gridView );

  // allocate and initialize objects for data representation
  ColorFunction colorFunction( gridView );
  ColorFunction update( gridView );
  ReconstructionSet reconstructionSet( gridView );
  Flags flags ( gridView );
  auto flagOperator = FlagOperator( eps );
  auto reconstruction = Dune::VoF::reconstruction( gridView, colorFunction, stencils );

  ColorFunction normalX( gridView );
  ColorFunction normalY( gridView );
  //ColorFunction normalZ( gridView );

  // File io
  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.path" ) << "/vof-" << std::to_string( level );
  createDirectory( path.str() );

  DataWriter vtkwriter ( gridView, "vof", path.str(), "" );
  vtkwriter.addCellData ( colorFunction, "celldata" );
  vtkwriter.addCellData ( curvatureSet, "curvature" );
  vtkwriter.addCellData ( curvatureError, "curvatureError" );
  vtkwriter.addCellData ( normalX, "nX" );
  vtkwriter.addCellData ( normalY, "nY" );
  //vtkwriter.addCellData ( normalZ, "nZ" );


  // Initial data
  Dune::VoF::average( colorFunction, problem, level );
  colorFunction.communicate();

  // Initial reconstruction
  flagOperator( colorFunction, flags );
  reconstruction( colorFunction, reconstructionSet, flags );
  curvatureOperator( colorFunction, reconstructionSet, flags, curvatureSet );

  double error = Dune::VoF::curvatureError( curvatureSet, flags, reconstructionSet, problem, curvatureError );

  if ( writeData )
  {
    for ( const auto &entity : elements( gridView ) )
    {
      normalX[ entity ] = reconstructionSet[ entity ].innerNormal()[ 0 ];
      normalY[ entity ] = reconstructionSet[ entity ].innerNormal()[ 1 ];
      //normalZ[ entity ] = reconstructionSet[ entity ].innerNormal()[ 2 ];
    }
    vtkwriter.write( 0 );
  }

  return error;
}

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  double lastL1Error = 0.0;

  //  create grid
  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;

  int level = parameters.get< int >( "grid.level" );
  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  grid.globalRefine( refineStepsForHalf * level );

  int maxLevel = parameters.get<int>( "grid.runs", 1 ) + level;

  std::ofstream errorsFile;
  errorsFile.open( "errors" );
  errorsFile << "#dx \terror " << std::endl;
  std::ofstream eocFile;
  eocFile.open( "eoc" );

  for ( ; level < maxLevel; ++level )
  {
    // start time integration
    double singleL1Error = algorithm( grid.leafGridView(), parameters );

    double L1Error = grid.comm().sum( singleL1Error );

    if ( grid.comm().rank() == 0 )
    {
      const double eoc = ( log( lastL1Error ) - log( L1Error ) ) / M_LN2;
      eocFile << std::setprecision(0) << "    $" << 8 * std::pow( 2, level ) << "^2$ \t& " << std::scientific << std::setprecision(2) << L1Error << " & " << std::fixed << eoc << " \\\\" << std::endl;
      errorsFile << 1.0 / 8.0 * std::pow( 2, -level ) << " \t" << L1Error << std::endl;

      if ( level > 0 )
        std::cout << "  EOC " << level << ": " << eoc << std::endl;

      std::cout << "Err " << level << ": " << L1Error << std::endl;
    }
    lastL1Error = L1Error;

    // refine
    parameters[ "grid.level" ] = std::to_string( level+1 );
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
