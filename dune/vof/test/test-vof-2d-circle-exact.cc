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
#include <dune/vof/evolution.hh>
#include <dune/vof/flagset.hh>
#include <dune/vof/flagging.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionset.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/curvature/generalheightfunctioncurvature.hh>
#include <dune/vof/curvature/cartesianheightfunctioncurvature.hh>
#include <dune/vof/interpolation.hh>

//- local includes
#include "average.hh"
#include "../curvatureset.hh"
#include "colorfunction.hh"
#include "errors.hh"
#include "io.hh"
#include "polygon.hh"
#include "problems/rotatingcircle.hh"
#include "velocity.hh"
#include "vtu.hh"

// filterReconstruction
// --------------------

template < class GridView, class ReconstructionSet, class Flags, class Polygon >
void filterReconstruction( const GridView &gridView, const ReconstructionSet &reconstructionSet, const Flags &flags, std::vector< Polygon > &io )
{
  using Coord = typename Polygon::Position;

  io.clear();
  for ( const auto& entity : elements( gridView ) )
  {
    if ( !flags.isMixed( entity ) )
      continue;

    auto is = intersect( Dune::VoF::makePolytope( entity.geometry() ), reconstructionSet[ entity ].boundary() );
    auto intersection = static_cast< typename decltype( is )::Result > ( is );

    if ( intersection == typename decltype( is )::Result() )
      continue;

    std::vector< Coord > vertices;
    for ( std::size_t i = 0; i < intersection.size(); ++i )
      vertices.push_back( intersection.vertex( i ) );

    if ( vertices.size() > 0 )
      io.push_back( Polygon( vertices, reconstructionSet[ entity ].innerNormal() ) );
  }
}


// TimeProvider
// ------------

template < class Comm >
class TimeProvider
{
 public:
  TimeProvider ( const double cfl, const double dt, const double start, const Comm &comm )
   : cfl_ ( cfl ), dt_( dt * cfl_ ), time_( start ), dtEst_ ( std::numeric_limits< double >::max() ), comm_( comm ) {}

  double time() const { return time_; }

  double deltaT() const { return dt_; }

  void next() {
    time_ += dt_;
    step_++;
    dt_ = comm_.min( dtEst_ );
    dtEst_ = std::numeric_limits< double >::max();
  }

  void provideTimeStepEstimate( const double dtEst )
  {
    dtEst_ = std::min( dtEst * cfl_, dtEst_ );
  }

 private:
  double cfl_, dt_, time_, dtEst_;
  int step_ = 0;
  Comm comm_;
};

// initTimeStep
// ---------------------

template < class GridView, class Velocity >
double initTimeStep( const GridView& gridView, Velocity &velocity )
{
  double dtMin = std::numeric_limits< double >::max();
  for( const auto &entity : elements( gridView ) )
  {
    const auto& geoEn = entity.geometry();
    for ( const auto& intersection : intersections( gridView, entity ) )
    {
      const auto& geoIs = intersection.geometry();
      velocity.bind( intersection );
      auto v = velocity( geoIs.local( geoIs.center() ) );
      dtMin = std::min( dtMin, geoEn.volume() / std::abs( intersection.integrationOuterNormal( 0.0 ) * v ) );
    }
  }
  return dtMin;
}

// algorithm
// ---------

template < class GridView >
double algorithm ( const GridView& gridView, const Dune::ParameterTree &parameters )
{
  using ColorFunction = ColorFunction< GridView >;
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
  using FlagSet = Dune::VoF::FlagSet< GridView >;
  using CurvatureSet = Dune::VoF::CurvatureSet< GridView >;


  using DataWriter = Dune::VTKSequenceWriter< GridView >;

  using Polygon = OutputPolygon< typename ReconstructionSet::DataType::Coordinate >;

  // Testproblem
  using ProblemType = RotatingCircle< double, GridView::dimensionworld >;
  using Velocity = Velocity< ProblemType, GridView >;

  ProblemType circle;

  // calculate dt
  int level = parameters.get< int >( "grid.level" );
  const double startTime = parameters.get< double >( "scheme.start", 0.0 );
  const double endTime = parameters.get< double >( "scheme.end", 10 );
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );

  Velocity velocity( circle, startTime );
  double dt = initTimeStep( gridView, velocity );

  const bool writeData = parameters.get< bool >( "io.writeData", dt );

  int saveNumber = 1;
  const double saveInterval = std::max( parameters.get< double >( "io.saveInterval", 1 ), dt );
  double nextSaveTime = saveInterval;

  // build domain references for each cell
  Stencils stencils( gridView );

  // allocate and initialize objects for data representation
  ColorFunction colorFunction( gridView );
  ColorFunction update( gridView );
  ReconstructionSet reconstructionSet( gridView );
  FlagSet flagSet ( gridView );

  auto flagOperator = Dune::VoF::FlagOperator< ColorFunction, FlagSet >( eps );
  auto reconstruction = Dune::VoF::reconstruction( gridView, colorFunction, stencils );
  auto evolution = Dune::VoF::evolution( gridView );


  using CurvatureOperator = Dune::VoF::CartesianHeightFunctionCurvature< GridView, Stencils, ColorFunction, ReconstructionSet, FlagSet >;
  CurvatureOperator curvatureOperator ( gridView, stencils );
  CurvatureSet curvatureSet( gridView );

  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( level );
  createDirectory( path.str() );

  ColorFunction normalX( gridView );
  ColorFunction normalY( gridView );

  DataWriter vtkwriter ( gridView, "vof", path.str(), "" );
  vtkwriter.addCellData ( colorFunction, "celldata" );
  vtkwriter.addCellData ( curvatureSet, "curvature" );
  vtkwriter.addCellData ( normalX, "nX" );
  vtkwriter.addCellData ( normalY, "nY" );

  std::vector< Polygon > recIO;
  VTUWriter< std::vector< Polygon > > vtuwriter( recIO );
  std::stringstream name;
  name.fill('0');
  name << "vof-rec-" << std::setw(5) << 0 << ".vtu";

  // Initial data
  Dune::VoF::circleInterpolation( circle.center( startTime ), circle.radius( startTime ), colorFunction );
  colorFunction.communicate();

  // Initial reconstruction
  flagOperator( colorFunction, flagSet );
  reconstruction( colorFunction, reconstructionSet, flagSet );

  curvatureOperator( colorFunction, reconstructionSet, flagSet, curvatureSet );
  for ( auto entity : elements( gridView ) )
  {
    normalX[ entity ] = reconstructionSet[ entity ].innerNormal()[ 0 ];
    normalY[ entity ] = reconstructionSet[ entity ].innerNormal()[ 1 ];
  }

  if ( writeData )
  {
    filterReconstruction( gridView, reconstructionSet, flagSet, recIO );
    vtkwriter.write( 0 );
    vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );
  }
  const auto comm = gridView.grid().comm();
  TimeProvider< decltype( comm ) > tp ( parameters.get< double >( "scheme.cflFactor" ), dt, startTime, comm );
  double error = 0;

  while ( tp.time() < endTime )
  {
    flagOperator( colorFunction, flagSet );

    reconstruction( colorFunction, reconstructionSet, flagSet );

    curvatureOperator( colorFunction, reconstructionSet, flagSet, curvatureSet );
    for ( auto entity : elements( gridView ) )
    {
      normalX[ entity ] = reconstructionSet[ entity ].innerNormal()[ 0 ];
      normalY[ entity ] = reconstructionSet[ entity ].innerNormal()[ 1 ];
    }

    if ( writeData && 2.0 * tp.time() > nextSaveTime - 0.5 * tp.deltaT() )
    {
      vtkwriter.write( tp.time() );

      filterReconstruction( gridView, reconstructionSet, flagSet, recIO );
      std::stringstream name_;
      name_.fill('0');
      name_ << "vof-rec-" << std::setw(5) << saveNumber << ".vtu" ;
      vtuwriter.write( Dune::concatPaths( path.str(), name_.str() ) );


      nextSaveTime += saveInterval;
      ++saveNumber;
    }

    Velocity velocity( circle, tp.time() );
    double dtEst = evolution( reconstructionSet, flagSet, velocity, tp.deltaT(), update );

    colorFunction.axpy( 1.0, update );
    colorFunction.communicate();

    double dt = tp.deltaT();
    tp.provideTimeStepEstimate( dtEst );
    tp.next();

    error += dt * Dune::VoF::exactL1Error( colorFunction, flagSet, reconstructionSet, circle.center( tp.time() ), circle.radius( tp.time() ) );
  }

  if ( tp.time() == startTime )
    error += Dune::VoF::exactL1Error( colorFunction, flagSet, reconstructionSet, circle.center( tp.time() ), circle.radius( tp.time() ) );

  return error;
}

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
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

  for ( ; level < maxLevel; ++level )
  {
    // start time integration
    double singleL1Error = algorithm( grid.leafGridView(), parameters );

    double L1Error = grid.comm().sum( singleL1Error );

    if ( grid.comm().rank() == 0 )
    {
      // print errors and eoc
      if ( level > 0 )
      {
        const double eoc = ( log( lastL1Error ) - log( L1Error ) ) / M_LN2;

        //if( eoc < 1.5 )
          //DUNE_THROW( Dune::InvalidStateException, "EOC check of 2d rotating circle problem failed.");

        std::cout << "  EOC " << level << ": " << eoc << std::endl;
      }

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
