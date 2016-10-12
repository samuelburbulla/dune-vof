#include "config.h"

//- C++ includes
#include <cmath>
#include <vector>

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/parametertreeparser.hh>

//- dune-vof includes
#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>


int main(int argc, char** argv)
try {

  Dune::MPIHelper::instance( argc, argv );

  using GridType = Dune::GridSelector::GridType;
  using Coordinate = Dune::FieldVector< double, GridType::dimension >;
  using Polytope = typename std::conditional< Coordinate::dimension == 2, Dune::VoF::Polygon< Coordinate >, Dune::VoF::Polyhedron< Coordinate > >::type;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

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
    const auto& entity = grid.leafGridView().begin< 0 >();
    const auto& geoEn = entity.geometry();
    const auto& polytope = Dune::VoF::makePolytope( geoEn );

    Coordinate normal ( 0.0 );
    normal[ 0 ] = -1.0;

    Dune::VoF::HalfSpace< Coordinate > hs1 ( normal, geoEn.corner(0) );
    Polytope interface1 = Dune::VoF::intersect( polytope, hs1 );
    if ( interface1.volume() > std::numeric_limits< double >::epsilon() )
      DUNE_THROW( Dune::InvalidStateException, "Intersection to empty part failed.");

    Dune::VoF::HalfSpace< Coordinate > hs2 ( normal, geoEn.center() );
    Polytope interface2 = Dune::VoF::intersect( polytope, hs2 );
    if ( std::abs( interface2.volume() - 0.5 * geoEn.volume() ) > std::numeric_limits< double >::epsilon() )
      DUNE_THROW( Dune::InvalidStateException, "Intersection to half of cell failed.");

    Dune::VoF::HalfSpace< Coordinate > hs3 ( normal, geoEn.corner(1) );
    Polytope interface3 = Dune::VoF::intersect( polytope, hs3 );
    if ( std::abs( interface2.volume() - 1.0 * geoEn.volume() ) > std::numeric_limits< double >::epsilon() )
      DUNE_THROW( Dune::InvalidStateException, "Intersection to full cell failed.");

    Coordinate normal2 ( 1.0 );
    normal2 /= normal2.two_norm();

    Dune::VoF::HalfSpace< Coordinate > hs4 ( normal2, geoEn.corner(0) );
    Polytope interface4 = Dune::VoF::intersect( polytope, hs4 );
    if ( std::abs( interface2.volume() - 0.5 * geoEn.volume() ) > std::numeric_limits< double >::epsilon() )
      DUNE_THROW( Dune::InvalidStateException, "Intersection to diagonal half of cell failed.");

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
