#include <config.h>

#include <cmath>

#include <iostream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/test/checkcommunicate.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkpartition.hh>
#include <dune/grid/test/gridcheck.hh>

#include <dune/vof/interfacegrid/grid.hh>
#include <dune/vof/reconstruction/modifiedswartz.hh>
#include <dune/vof/reconstruction/modifiedyoungssecondorder.hh>
#include <dune/vof/reconstructionSet.hh>
#include <dune/vof/stencil/vertexneighborsstencil.hh>

#include "average.hh"
#include "colorfunction.hh"
#include "problems/ellipse.hh"


// main
// ----

int main ( int argc, char** argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  typedef Dune::GridSelector::GridType Grid;

  // set parameters
  Dune::ParameterTree parameters;
  Dune::ParameterTreeParser::readINITree( "parameter.ini", parameters );
  Dune::ParameterTreeParser::readOptions( argc, argv, parameters );

  //  create grid
  Dune::GridPtr< Grid > grid( std::to_string( Grid::dimension ) + "dgrid.dgf" );
  grid->loadBalance();

  int level = parameters.get< int >( "grid.level" );
  const int refineStepsForHalf = Dune::DGFGridInfo< Grid >::refineStepsForHalf();
  grid->globalRefine( refineStepsForHalf * level );

  typedef Grid::LeafGridView GridView;
  const GridView gridView = grid->leafGridView();

  // obtain initial color function
  Dune::FieldVector< double, GridView::dimensionworld > axis( { 1.0 / std::sqrt( 2 ), 1.0 / std::sqrt( 2 ) } );
  Ellipse< double, GridView::dimensionworld > problem( { axis, Dune::VoF::generalizedCrossProduct( axis ) }, { 0.2, 0.5 } );
  ColorFunction< GridView > colorFunction( gridView );
  Dune::VoF::average( colorFunction, problem );
  colorFunction.communicate();

  // create interface grid
  typedef Dune::VoF::ReconstructionSet< GridView > ReconstructionSet;
  typedef Dune::VoF::VertexNeighborsStencil< GridView > Stencils;
  typedef Dune::VoF::ModifiedSwartzReconstruction< ColorFunction< GridView >, ReconstructionSet, Stencils, Dune::VoF::ModifiedYoungsSecondOrderReconstruction< ColorFunction< GridView >, ReconstructionSet, Stencils > > Reconstruction;
  typedef Dune::VoF::InterfaceGrid< Reconstruction > InterfaceGrid;

  Stencils stencils( gridView );
  InterfaceGrid interfaceGrid( colorFunction, stencils );

  // perform check
  gridcheck( interfaceGrid );
  checkIterators( interfaceGrid.leafGridView() );
  checkPartitionType( interfaceGrid.leafGridView() );
  checkIntersectionIterator( interfaceGrid );
  checkCommunication( interfaceGrid, -1, std::cout );

  // perform VTK output
  Dune::VTKWriter< GridView > vtkWriter( gridView );
  vtkWriter.addCellData( colorFunction, "color" );
  vtkWriter.write( "test-interfacegrid-color" );

  Dune::VTKWriter< InterfaceGrid::LeafGridView > interfaceVtkWriter( interfaceGrid.leafGridView() );
  interfaceVtkWriter.write( "test-interfacegrid-reconstruction" );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
