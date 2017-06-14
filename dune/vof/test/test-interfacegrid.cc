#include <config.h>

#include <cmath>

#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/test/checkcommunicate.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkpartition.hh>
#include <dune/grid/test/gridcheck.hh>

#include <dune/vof/colorfunction.hh>
#include <dune/vof/interfacegrid/grid.hh>
#include <dune/vof/reconstruction.hh>
#include <dune/vof/reconstructionset.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>

#include "average.hh"
#include "problems/ellipse.hh"


template< class C, class R >
auto ifGrid ( const C& color, R rOp ) -> Dune::VoF::InterfaceGrid< R >
{
  return Dune::VoF::InterfaceGrid< R >( color, std::move( rOp ) );
}

// main
// ----

int main ( int argc, char** argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  typedef Dune::GridSelector::GridType Grid;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  //  create grid
  Dune::GridPtr< Grid > grid( std::to_string( Grid::dimension ) + "dgrid.dgf" );
  grid->loadBalance();

  int level = parameters.get< int >( "grid.level" );
  const int refineStepsForHalf = Dune::DGFGridInfo< Grid >::refineStepsForHalf();
  grid->globalRefine( refineStepsForHalf * level );

  typedef Grid::LeafGridView GridView;
  const GridView gridView = grid->leafGridView();

  typedef Dune::FieldVector< double, GridView::dimensionworld > Coordinate;

  // obtain initial color function
#if GRIDDIM == 2
  Coordinate axis( 1.0 / std::sqrt( static_cast< double > ( GridView::dimensionworld ) ) );
  Ellipse< double, GridView::dimensionworld > problem( { axis, Dune::VoF::generalizedCrossProduct( axis ) }, { 0.2, 0.5 } );
#elif GRIDDIM == 3
  Coordinate axis1 ( { 1.0, 0.0, 0.0 } );
  Coordinate axis2 ( { 0.0, 1.0, 0.0 } );
  Coordinate axis3 ( { 0.0, 0.0, 1.0 } );
  Ellipse< double, GridView::dimensionworld > problem( { axis1, axis2, axis3 }, { 0.2, 0.2, 0.4 } );
#endif

  Dune::VoF::ColorFunction< GridView > colorFunction( gridView );
  Dune::VoF::Average< Ellipse< double, GridView::dimensionworld > > average ( problem );
  average( colorFunction );

  typedef Dune::VoF::VertexNeighborsStencil< GridView > Stencils;
  Stencils stencils( gridView );

  // create interface grid
  // using Reconstruction = decltype( Dune::VoF::reconstruction( stencils ) );
  // typedef Dune::VoF::InterfaceGrid< Reconstruction > InterfaceGrid;

  // InterfaceGrid interfaceGrid( colorFunction, stencils );

  auto interfaceGrid = ifGrid( colorFunction, Dune::VoF::reconstruction( stencils ) );
  using InterfaceGrid = std::decay_t< decltype( interfaceGrid ) >;

  // perform check
  if( InterfaceGrid::dimension < 2 )
    gridcheck( interfaceGrid );
  checkIterators( interfaceGrid.leafGridView() );
  checkPartitionType( interfaceGrid.leafGridView() );
  if( InterfaceGrid::dimension < 2 )
    checkIntersectionIterator( interfaceGrid );
  checkCommunication( interfaceGrid, -1, std::cout );

  Dune::VoF::DataSet< GridView, std::size_t > indexSet( gridView );
  for ( const auto &entity : elements( gridView ) )
    indexSet[ entity ] = gridView.indexSet().index( entity );

  // perform VTK output
  Dune::VTKWriter< GridView > vtkWriter( gridView );
  vtkWriter.addCellData( colorFunction, "color" );
  vtkWriter.addCellData( indexSet, "index" );
  vtkWriter.write( "test-interfacegrid-color" );

  std::vector< double > normals( interfaceGrid.leafGridView().indexSet().size( 0 ) * GridView::dimensionworld );
  for ( const auto &entity : elements( interfaceGrid.leafGridView() ) )
    for ( std::size_t i = 0; i < GridView::dimensionworld; ++i )
      normals[ interfaceGrid.leafGridView().indexSet().index( entity ) * GridView::dimensionworld + i ] = interfaceGrid.dataSet().reconstructionSet()[ entity.impl().hostElement() ].innerNormal()[ i ];

  Dune::VTKWriter< InterfaceGrid::LeafGridView > interfaceVtkWriter( interfaceGrid.leafGridView() );
  interfaceVtkWriter.addCellData( normals, "normals", 3 );
  interfaceVtkWriter.write( "test-interfacegrid-reconstruction" );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
