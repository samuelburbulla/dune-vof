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
#include "exactevolution.hh"
#include "polygon.hh"
#include "io.hh"
#include "problems/rotatingcircle.hh"
#include "utility.hh"
#include "vtu.hh"

// filterReconstruction
// --------------------

template < class GridView, class ReconstructionSet, class Polygon >
void filterReconstruction( const GridView &gridView, const ReconstructionSet &reconstructionSet, std::vector< Polygon > &io )
{
  using Coord = typename Polygon::Position;

  io.clear();
  for ( const auto& entity : elements( gridView ) )
  {
    auto is = intersect( Dune::VoF::makePolytope( entity.geometry() ), reconstructionSet[ entity ].boundary() );
    auto intersection = static_cast< typename decltype( is )::Result > ( is );

    if ( intersection == typename decltype( is )::Result() )
      continue;

    std::vector< Coord > vertices;
    for ( std::size_t i = 0; i < intersection.size(); ++i )
      vertices.push_back( intersection.vertex( i ) );

    if ( vertices.size() > 0 )
      io.push_back( Polygon( vertices ), reconstructionSet[ entity ].innerNormal() );
  }
}


// TimeProvider
// ------------

class TimeProvider
{
 public:
  TimeProvider ( const double cfl, const double dt, const double start )
   : cfl_ ( cfl ), dt_( dt * cfl_ ), time_( start ), dtEst_ ( std::numeric_limits< double >::max() ) {}

  double time() const { return time_; }

  double deltaT() const { return dt_; }

  void next() {
    time_ += dt_;
    step_++;
    dt_ = dtEst_;
    dtEst_ = std::numeric_limits< double >::max();
  }

  void provideTimeStepEstimate( const double dtEst )
  {
    dtEst_ = std::min( dtEst * cfl_, dtEst_ );
  }

 private:
  double cfl_, dt_, time_, dtEst_;
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
double algorithm ( const GridView& gridView, const Dune::ParameterTree &parameters, const double cfl )
{
  using DomainVector = Dune::FieldVector< double, GridView::dimensionworld >;
  using ColorFunction = ColorFunction< GridView >;
  using Stencils = Dune::VoF::VertexNeighborsStencil< GridView >;
  using ReconstructionSet = Dune::VoF::ReconstructionSet< GridView >;
  using Flags = Dune::VoF::Flags< GridView >;

  using DataWriter = Dune::VTKSequenceWriter< GridView >;

  using Polygon = OutputPolygon< typename ReconstructionSet::DataType::Coordinate >;

  // Testproblem
  using ProblemType = RotatingCircle< double, GridView::dimensionworld >;
  ProblemType circle;

  // calculate dt
  int level = parameters.get< int >( "grid.level" );
  double dt = initTimeStep( gridView, [ &circle ] ( const auto &x ) { DomainVector rot; circle.velocityField( x, 0.0, rot ); return rot; } );
  const double eps = parameters.get< double >( "scheme.epsilon", 1e-6 );
  const bool writeData = parameters.get< bool >( "io.writeData", 0 );

  // build domain references for each cell
  Stencils stencils( gridView );

  // allocate and initialize objects for data representation
  ColorFunction colorFunction( gridView );
  ColorFunction errorFunction( gridView );
  ColorFunction uhExact( colorFunction );
  ColorFunction update( gridView );
  Flags flags ( gridView );
  auto exactEvolution = Dune::VoF::exactEvolution( colorFunction );

  std::stringstream path;
  path << "./" << parameters.get< std::string >( "io.folderPath" ) << "/vof-" << std::to_string( level );
  createDirectory( path.str() );

  DataWriter vtkwriter ( gridView, "vof", path.str(), "~/dune" );
  vtkwriter.addCellData ( colorFunction, "uh" );
  vtkwriter.addCellData ( uhExact, "u" );

  std::vector< Polygon > recIO;
  VTUWriter< std::vector< Polygon > > vtuwriter( recIO );

  // Initial data
  Dune::VoF::circleInterpolation( circle.center( 0.0 ), circle.radius( 0.0 ), colorFunction );
  Dune::VoF::circleInterpolation( circle.center( 0.0 ), circle.radius( 0.0 ), uhExact );

  // Initial flagging
  flags.reflag( colorFunction, eps );

  if ( writeData )
    vtkwriter.write( 0 );

  TimeProvider tp ( cfl, dt, 0.0 );

  auto velocity = [ &circle ] ( const auto &x, const auto &t ) { DomainVector rot; circle.velocityField( x, t, rot ); return rot; };

  exactEvolution( colorFunction, circle, velocity, tp, update, flags );
  colorFunction.axpy( 1.0, update );
  tp.next();

  Dune::VoF::circleInterpolation( circle.center( tp.time() ), circle.radius( tp.time() ), uhExact );

  if ( writeData )
    vtkwriter.write( 1 );

  ReconstructionSet reconstructionSet( gridView );
  auto reconstruction = Dune::VoF::reconstruction( gridView, colorFunction, stencils );
  flags.reflag( colorFunction, eps );
  reconstruction( colorFunction, reconstructionSet, flags );

  return Dune::VoF::exactL1Error( colorFunction, flags, reconstructionSet, circle.center( tp.time() ), circle.radius( tp.time() ) );
}

int main(int argc, char** argv)
try {
  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  for ( double cfl = 0.25; cfl >= 0.25; cfl *= 0.5 )
  {
    std::cout << std::endl << "CFL=" << cfl << std::endl;

    double lastError = 0.0;

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
      auto error = algorithm( grid.leafGridView(), parameters, cfl );


      // print errors and eoc
      if ( level > 0 )
      {
        const double eoc = ( log( lastError ) - log( error ) ) / M_LN2;

        std::cout << "  EOC " << level << ": " << eoc << std::endl;
      }

      std::cout << "Err " << level << ": " << error << std::endl;

      lastError = error;

      // refine
      parameters[ "grid.level" ] = std::to_string( level+1 );
      grid.globalRefine( refineStepsForHalf );

    }

    parameters[ "grid.level" ] = std::to_string( 0 );
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
