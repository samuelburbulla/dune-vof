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

// dune-vof includes
#include <dune/vof/evolution.hh>
#include <dune/vof/vertexneighborsstencil.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/reconstruction.hh>

// local includes
#include "fvscheme.hh"
#include "model.hh"
#include "problems.hh"

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


// algorithm
// ---------

template< class Grid >
std::tuple< double, double > algorithm ( Grid &grid, int level, double start, double end, double cfl, double eps )
{
  using GridType = Grid;
  using GridPartType =
    Dune::Fem::AdaptiveLeafGridPart< GridType >;

  using FunctionSpaceType =
    Dune::Fem::FunctionSpace< double, double, GridPartType::dimensionworld, 1 >;

  using DiscreteFunctionSpaceType =
    Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >;
  using DiscreteFunctionType =
    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType>;

  using ProblemType =
    Problem< FunctionSpaceType >;
  using SolutionType =
    Dune::Fem::InstationaryFunction< ProblemType >;
  using ModelType =
    Model< ProblemType >;

  using GridSolutionType =
    Dune::Fem::GridFunctionAdapter< SolutionType, GridPartType >;

  using TimeProviderType =
    Dune::Fem::FixedStepTimeProvider< typename GridType::CollectiveCommunication >;

  using DataIOTupleType =
    std::tuple< DiscreteFunctionType *, GridSolutionType * >;
  using DataOutputType =
    Dune::Fem::DataOutput< GridType, DataIOTupleType >;



  using Stencils = 
    Dune::VoF::VertexNeighborsStencil< GridPartType >;
  using ReconstructionSet =
    Dune::VoF::ReconstructionSet< GridPartType >;
  using SchemeType =
    Dune::VoF::Evolution< ReconstructionSet, DiscreteFunctionType >;



  using L1NormType =
    Dune::Fem::L1Norm< GridPartType >;
  using L2NormType =
    Dune::Fem::L2Norm< GridPartType >;

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType space( gridPart );

  // build domain references for each cell
  Stencils stencils( gridPart );

  ModelType model;
  ReconstructionSet reconstructions( gridPart );
  auto reconstruction = Dune::VoF::reconstruction( grid, space, stencils, eps );

  SchemeType scheme( eps );

  L1NormType l1norm( gridPart );
  L2NormType l2norm( gridPart );

  DiscreteFunctionType uh( "uh", space );
  DiscreteFunctionType update( "update", space );
  
  SolutionType solution( model.problem(), start );
  GridSolutionType u( "solution", solution, gridPart, 5 );

  double timeStep = std::pow( 2, -(3 + level ) );
  timeStep *= cfl;
  TimeProviderType timeProvider( start, timeStep, gridPart.comm() );

  DataIOTupleType dataIOTuple = std::make_tuple( &uh, &u );
  DataOutputType dataOutput( grid, dataIOTuple, timeProvider, DataOutputParameters( level ) );

  
  Dune::Fem::interpolate( u, uh );
  uh.communicate();
  dataOutput.write( timeProvider );

  using RangeType = typename FunctionSpaceType::RangeType;
  auto velocity = [ &timeProvider, &model ] ( const auto &x ) { RangeType u; model.problem().evaluate( x, timeProvider.time(), u ); return u; };

  for( ; timeProvider.time() < end; timeProvider.next() )
  {
    if( Dune::Fem::Parameter::verbose() )
      std::cerr << "time step = " << timeProvider.timeStep() << ", "
                << "time = " << timeProvider.time() << ", "
                << "dt = " << timeProvider.deltaT() << std::endl;
    
    Dune::VoF::reconstruction( space, reconstructions );
    scheme( uh, reconstructions, velocity, timeProvider.deltaT(), update );
    uh.axpy( timeProvider.deltaT(), update );

    solution.setTime( timeProvider.time() );

    dataOutput.write( timeProvider );
  }

  return std::make_tuple( l1norm.distance( u, uh ), l2norm.distance( u, uh ) );
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

   //  create grid
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
