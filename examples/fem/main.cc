#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
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
#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/operator/projection/l2projection.hh>

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
#include "../problems/rotatingcircle.hh"


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

  std::string prefix () const
  {
    std::stringstream s;
    s << "vof-fem-" << level_ << "-";
    return s.str();
  }

private:
  int level_;
};


// Problem
// -------
template < class Base, class FunctionSpace >
struct Problem : public Base
{
  using FunctionSpaceType = FunctionSpace;
};


// algorithm
// ---------

template< class Grid >
std::tuple< double, double > algorithm ( Grid &grid, int level, double start, double end, double cfl, double eps )
{
  using GridType = Grid;
  using GridPartType =
    Dune::Fem::LeafGridPart< GridType >;

  using FunctionSpaceType =
    Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;

  using DiscreteFunctionSpaceType =
    Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType =
    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;

  // Testproblem
  using ProblemType =
    Problem< RotatingCircle< double, GridPartType::dimensionworld >, FunctionSpaceType >;

  using SolutionType =
    Dune::Fem::InstationaryFunction< ProblemType >;

  using GridSolutionType =
    Dune::Fem::GridFunctionAdapter< SolutionType, GridPartType >;

  using TimeProviderType =
    Dune::Fem::FixedStepTimeProvider< typename GridType::CollectiveCommunication >;

  using ReconstructionSet =
    Dune::VoF::ReconstructionSet< GridPartType >;
  using Polygon =
    Polygon< typename ReconstructionSet::Reconstruction::Coordinate >;

  using DataIOTupleType =
    std::tuple< DiscreteFunctionType * >;
  using DataOutputType =
    Dune::Fem::DataOutput< GridType, DataIOTupleType >;
  using RecOutputType =
    ReconstructionWriter< ReconstructionSet, Polygon >;

  // Stencil
  using Stencils =
    Dune::VoF::VertexNeighborsStencil< GridPartType >;

  using L1NormType =
    Dune::Fem::L1Norm< GridPartType >;
  using L2NormType =
    Dune::Fem::L2Norm< GridPartType >;

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType space( gridPart );

  Stencils stencils( gridPart );

  ProblemType problem;
  ReconstructionSet reconstructions( gridPart );

  L1NormType l1norm( gridPart );
  L2NormType l2norm( gridPart );

  DiscreteFunctionType uh( "uh", space );
  DiscreteFunctionType update( "update", space );

  auto cuh = Dune::VoF::discreteFunctionWrapper( uh );
  auto cupdate = Dune::VoF::discreteFunctionWrapper( update );

  auto evolution = Dune::VoF::evolution( reconstructions, cuh, eps );
  auto reconstruction = Dune::VoF::reconstruction( gridPart, cuh, stencils );
  auto flags = Dune::VoF::flags( gridPart );

  SolutionType solution( problem, start );
  GridSolutionType u( "solution", solution, gridPart, 9 );

  double timeStep = std::pow( 2, -(3 + level ) );
  timeStep *= cfl;
  timeStep /= ProblemType::maxVelocity();
  TimeProviderType timeProvider( start, timeStep, gridPart.comm() );

  DataIOTupleType dataIOTuple = std::make_tuple( &uh );
  DataOutputType dataOutput( grid, dataIOTuple, timeProvider, DataOutputParameters( level ) );
  RecOutputType recOutput ( level );

  Dune::Fem::L2Projection < GridSolutionType, DiscreteFunctionType > l2projection ( 15 );
  l2projection( u, uh );
  uh.communicate();

  flags.reflag( cuh, eps );
  reconstruction( cuh, reconstructions, flags );

  recOutput.write ( reconstructions );
  dataOutput.write( timeProvider );


  using DomainType = typename FunctionSpaceType::DomainType;
  auto velocity = [ &timeProvider, &problem ] ( const auto &x ) { DomainType u; problem.velocityField( x, timeProvider.time(), u ); return u; };

  std::size_t count = 0;
  for( ; timeProvider.time() <= end; )
  {
    if( Dune::Fem::Parameter::verbose() )
      std::cerr << "time step = " << timeProvider.timeStep() << ", "
                << "time = " << timeProvider.time() << ", "
                << "dt = " << timeProvider.deltaT() << std::endl;

    flags.reflag( cuh, eps );

    reconstruction( cuh, reconstructions, flags );
    evolution( cuh, reconstructions, velocity, timeProvider.deltaT(), cupdate, flags );
    uh.axpy( 1.0, update );

    timeProvider.next();

    solution.setTime( timeProvider.time() );

    if ( dataOutput.willWrite( timeProvider ) )
    {
      recOutput.write ( reconstructions );
      dataOutput.write( timeProvider );

      if( Dune::Fem::Parameter::verbose() )
          std::cout << "written reconstructions count=" << count << std::endl;
    }

    dataOutput.write( timeProvider );
  }

  DiscreteFunctionType uhEnd( "", space );
  l2projection( u, uhEnd );
  return std::make_tuple( l1norm.distance( uh, uhEnd ), l2norm.distance( uh, uhEnd ) );
}


int main(int argc, char** argv)
try {

  Dune::Fem::MPIManager::initialize( argc, argv );

   // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  using GridType = Dune::GridSelector::GridType;

  std::stringstream gridFile;
  gridFile << GridType::dimension << "dgrid.dgf";

   // create grid
  Dune::GridPtr< GridType > gridPtr( gridFile.str() );
  gridPtr->loadBalance();
  GridType& grid = *gridPtr;

  const int level = Dune::Fem::Parameter::getValue< int >( "level", 0 );
  const int repeats = Dune::Fem::Parameter::getValue< int >( "repeats", 2 );
  const int refineStepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();

  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  const double startTime = Dune::Fem::Parameter::getValue< double >( "start", 0.0 );
  const double endTime = Dune::Fem::Parameter::getValue< double >( "end", 2.5 );
  const double cfl = Dune::Fem::Parameter::getValue< double >( "cfl", 1.0 );
  const double eps = Dune::Fem::Parameter::getValue< double >( "eps", 1e-9 );

  auto oldErrors = algorithm( grid, level, startTime, endTime, cfl, eps );

  std::cout << "L1 Error: " << std::setw( 16 ) << std::get< 0 >( oldErrors ) << "   ";
  std::cout << "L2 Error: " << std::setw( 16 ) << std::get< 1 >( oldErrors ) << std::endl;

  for( int step = level+1; step <= level+repeats; ++step )
  {
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    auto newErrors = algorithm( grid, step, startTime, endTime, cfl, eps );

    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      const double l1eoc = log( std::get< 0 >( oldErrors )  / std::get< 0 >( newErrors ) ) / M_LN2;
      const double l2eoc = log( std::get< 1 >( oldErrors )  / std::get< 1 >( newErrors ) ) / M_LN2;

      std::cout << "L1 Error: " << std::setw( 16 ) << std::get< 0 >( newErrors ) << "   ";
      std::cout << "L2 Error: " << std::setw( 16 ) << std::get< 1 >( newErrors ) << std::endl;
      std::cout << "L1 EOC( " << std::setw( 2 ) << step << " ) = " << std::setw( 11 ) << l1eoc << "   ";
      std::cout << "L2 EOC( " << std::setw( 2 ) << step << " ) = " << std::setw( 11 ) << l2eoc << std::endl;
    }

    oldErrors = newErrors;
  }

  return 0;
}
catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
}
