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
#include <dune/vof/flags.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/geometry/intersect.hh>

//- local includes
#include "average.hh"
#include "colorfunction.hh"
#include "errors.hh"
#include "io.hh"
#include "problem.hh"
#include "utility.hh"
#include "vtu.hh"

// filterReconstruction
// --------------------

template < class GridView, class ReconstructionSet, class Polygon >
void filterReconstruction( const GridView &gridView, const ReconstructionSet &reconstructionSet, std::vector< Polygon > &io )
{
  io.clear();
  for ( const auto& entity : elements( gridView ) )
  {
    Polygon is = intersect( Dune::VoF::makePolytope( entity.geometry() ), reconstructionSet[ entity ] );

    if ( is.size() != 0 )
      io.push_back( Polygon( is ) );
  }
}


// FixedStepTimeProvider
// ---------------------

class FixedStepTimeProvider
{
 public:
  FixedStepTimeProvider ( const double dt, const double start = 0.0 )
   : dt_( dt ), time_( start ) {}

  double time() const { return time_; }

  double deltaT() const { return dt_; }

  void next() {
    time_ += dt_;
    step_++;
  }

  void provideTimeStepEstimate( const double dtEst ) {}

 private:
  double dt_, time_;
  int step_ = 0;
};


// initTimeStep
// ---------------------

template < class GridView, class Velocity >
double initTimeStep( const GridView& gridView, const Velocity &velocity )
{
  double dtMin = std::numeric_limits< double >::max();
  for( const auto &entity : elements( gridView ) )
  {
    const auto geoEn = entity.geometry();
    auto v = velocity( geoEn.center() );
    dtMin = std::min( dtMin, geoEn.volume() / v.two_norm() );
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
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
  using Flags = Dune::VoF::Flags< GridView >;

  using DataWriter = Dune::VTKSequenceWriter< GridView >;

  using Polygon = Dune::VoF::Polygon< typename ReconstructionSet::Reconstruction::Coordinate >;

  // Testproblem
  using ProblemType = RotatingCircle< double, GridView::dimensionworld >;
  ProblemType circle;

  // calculate dt
  int level = parameters.get< int >( "grid.level" );
  double dt = parameters.get< double >( "scheme.cflFactor" )
    * initTimeStep( gridView, [ &circle ] ( const auto &x ) { DomainVector rot; circle.velocityField( x, 0.0, rot ); return rot; } );
  const double startTime = parameters.get< double >( "scheme.start", 0.0 );
  const double endTime = parameters.get< double >( "scheme.end", 10 );
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );
  const bool writeData = parameters.get< bool >( "io.writeData", 0 );

  int saveNumber = 1;
  const double saveInterval = std::max( parameters.get< double >( "io.saveInterval", 1 ), dt );
  double nextSaveTime = saveInterval;

  // build domain references for each cell
  Stencils stencils( gridView );

  // allocate and initialize objects for data representation
  ColorFunction colorFunction( gridView );
  ColorFunction update( gridView );
  ReconstructionSet reconstructionSet( gridView );
  Flags flags ( gridView );
  auto reconstruction = Dune::VoF::reconstruction( gridView, colorFunction, stencils );
  auto evolution = Dune::VoF::evolution(  reconstructionSet, colorFunction, eps );

  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( level );
  createDirectory( path.str() );

  DataWriter vtkwriter ( gridView, "vof", path.str(), "~/dune" );
  vtkwriter.addCellData ( colorFunction, "celldata" );

  std::vector< Polygon > recIO;
  VTUWriter< std::vector< Polygon > > vtuwriter( recIO );
  std::stringstream name;
  name.fill('0');
  name << "vof-rec-" << std::setw(5) << 0 << ".vtu";

  // Initial data
  Dune::VoF::circleInterpolation( circle.center( startTime ), circle.radius( startTime ), colorFunction );

  // Initial reconstruction
  flags.reflag( colorFunction, eps );
  reconstruction(  colorFunction, reconstructionSet, flags );

  if ( writeData )
  {
    filterReconstruction( gridView, reconstructionSet, recIO );
    vtkwriter.write( 0 );
    vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );
  }

  FixedStepTimeProvider tp ( dt, startTime );

  auto velocity = [ &tp, &circle ] ( const auto &x ) { DomainVector rot; circle.velocityField( x, tp.time(), rot ); return rot; };

  while ( tp.time() < endTime )
  {
    flags.reflag( colorFunction, eps );

    reconstruction( colorFunction, reconstructionSet, flags );
    evolution( colorFunction, reconstructionSet, velocity, tp, update, flags );
    colorFunction.axpy( 1.0, update );

    tp.next();

    if ( writeData && 2.0 * std::abs( tp.time() - nextSaveTime ) < tp.deltaT() )
    {
      vtkwriter.write( tp.time() );

      filterReconstruction( gridView, reconstructionSet, recIO );
      std::stringstream name_;
      name_.fill('0');
      name_ << "vof-rec-" << std::setw(5) << saveNumber << ".vtu" ;
      vtuwriter.write( Dune::concatPaths( path.str(), name_.str() ) );


      nextSaveTime += saveInterval;
      ++saveNumber;
    }

  }

  return Dune::VoF::exactL1Error( colorFunction, reconstructionSet, circle.center( tp.time() ), circle.radius( tp.time() ) );
}

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  double lastL1Error;

  //  create grid
  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;

  int level = parameters.get< int >( "grid.level" );
  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  grid.globalRefine( refineStepsForHalf * level );

  int numRuns = parameters.get<int>( "grid.runs", 1 );

  for ( int i = 0; i < numRuns; ++i )
  {
    // start time integration
    auto L1Error = algorithm( grid.leafGridView(), parameters );

      // print errors and eoc
    if ( i > 0 )
    {
      const double eoc = log( lastL1Error / L1Error ) / M_LN2;

      if( eoc < 1.5 )
        DUNE_THROW( Dune::InvalidStateException, "EOC check of 2d rotating circle problem failed.");

      std::cout << "EOC " << i << ": " << eoc << std::endl;
    }

    lastL1Error = L1Error;

    // refine
    parameters[ "grid.level" ] = std::to_string( ++level );
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
