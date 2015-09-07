#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <fstream>
#include <vector>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/yaspgrid.hh>
#include <cmath>
#include <functional>
#include <stdexcept>


#include <dune/vof/flagCells.hh>
#include <dune/vof/errors.hh>

#include <dune/vof/params.hh>
#include <dune/vof/io.hh>
#include <dune/vof/domain.hh>
#include <dune/vof/initialize.hh>


#include <dune/vof/evolve.hh>
#include <dune/vof/reconstruct.hh>


typedef Dune::FieldVector<double,2> fvector;

struct Parameters
{
  	// concetration eps for mixed cells
  	double eps = 1e-6;
  	// number of cells for the cartesian grid in one direction
  	int numberOfCells = 32;


  	// params for timeloop
	double dtAlpha = 0.1;
	double tEnd = 2.5;

	// params for saving data
	double saveInterval = 0.1;

	std::string folderPath = "results/";
};

void filterReconstruction( const std::vector< std::array<fvector,3 > > &rec, std::vector< std::vector<fvector> > &io )
{
	io.clear();
	for( auto && it: rec)
		if( it[ 2 ] == fvector( 0.0 ) )
		{
			std::vector< fvector > points;
			points.push_back( it[ 0 ] );
			points.push_back( it[ 1 ] );
			io.push_back( std::move( points ) );
		}
}



template < class Grid >
std::tuple<double, double> algorithm ( const Grid& grid, const Parameters &params )
{
		// build domain references for each cell
	Dune::VoF::DomainOfCells<Grid> domain ( grid );

	auto gridView = grid.leafGridView();

	int n = grid.leafGridView().size(0);



	// allocate and initialize vectors for data representation
	std::vector<double> concentration ( n );
	std::vector< std::array<fvector,3> > reconstruction( n );
	std::vector< std::vector<fvector> > recIO;
	std::vector<bool> cellIsMixed ( n );
	std::vector<bool> cellIsActive ( n );
	std::vector<fvector> velocityField ( n );


	Dune::VoF::initialize( grid, concentration, Dune::VoF::f0<fvector> );

	// calculate dt
	double dt = params.dtAlpha * ( 1.0 / params.numberOfCells ) / Dune::VoF::psiMax();


	double t = 0;
	int k = 0;


	Dune::VoF::flagCells( gridView, concentration, reconstruction, domain, cellIsMixed, cellIsActive, params.eps );

	Dune::VoF::clearReconstruction( reconstruction );
	Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, params.eps );
	filterReconstruction( reconstruction, recIO );

	int saveNumber = 1;
	double saveStep = params.saveInterval;


	// VTK Writer
	std::stringstream path;
	path << "./results/vof-" << std::to_string( params.numberOfCells );

	Dune::VoF::createDirectory( path.str() );

	Dune::VTKSequenceWriter<typename Grid::LeafGridView> vtkwriter ( gridView, "vof", path.str(), "~/dune" );

	vtkwriter.addCellData ( concentration, "celldata" );
	vtkwriter.addCellData ( cellIsMixed, "cellmixed" );
	vtkwriter.addCellData ( cellIsActive, "cellactive" );

	vtkwriter.write( 0 );

	std::cerr << std::endl << params.numberOfCells << " cells" << std::endl;

	std::vector<double> concentrationEnd( n, 0 );

	while ( t < params.tEnd )
	{
		++k;

		auto psit = std::bind( Dune::VoF::psi<fvector>, std::placeholders::_1, t);
		Dune::VoF::L1projection( grid, velocityField, psit );


		Dune::VoF::clearReconstruction( reconstruction );

		Dune::VoF::flagCells( grid.leafGridView(), concentration, reconstruction, domain, cellIsMixed, cellIsActive, params.eps );
		Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, params.eps );
		Dune::VoF::evolve( grid, concentration, reconstruction, domain, params.numberOfCells, t, dt, params.eps, cellIsMixed, cellIsActive, velocityField );
		filterReconstruction( reconstruction, recIO );

		t += dt;

		if ( t >= saveStep )
		{
		  vtkwriter.write( t );

		  saveStep += params.saveInterval;
		  ++saveNumber;
		}

		std::cerr << "\r" << "[" << (int)(( t / params.tEnd ) * 100) << "%]";
		//std::cerr << "s=" << grid.size(0) << " k=" << k << " t=" << t << " dt=" << dt << " saved=" << saveNumber-1 << std::endl;

	}

	auto ft = std::bind( Dune::VoF::f<fvector>, std::placeholders::_1, t);

	double L1 = Dune::VoF::l1error( grid, concentration, ft );
	double L2 = Dune::VoF::l2error( grid, concentration, ft );

	return std::tuple<double, double> ( L1, L2 );
}

int main(int argc, char** argv)
{

  	Dune::MPIHelper::instance( argc, argv );

  	try {


  		// type declarations
		const int dim = 2;
  		typedef Dune::YaspGrid<dim> GridType;
		typedef typename Dune::FieldVector<double,dim> fvector;

		// set parameters
  		Parameters params;

		Dune::VoF::handleInput ( argc, argv, params );


		std::tuple<double, double> lastErrorTuple;

		std::cout << "Cells \t\t L1 \t eoc \t\t L2 \t eoc" << std::endl << std::endl;


		for ( std::size_t i = 0; i < 7; ++i )
		{

			// build Grid
			fvector upper( 1.0 );
			Dune::array<int,dim> noc;
			std::fill( noc.begin(), noc.end(), params.numberOfCells );
			GridType grid( upper, noc );

			// start time integration
			auto errorTuple = algorithm( grid,  params );

		    // print errors and eoc
			if ( i > 0 )
			{
				const double eocL1 = log( std::get< 0 >( lastErrorTuple )  / std::get< 0 >( errorTuple ) ) / M_LN2;
				const double eocL2 = log( std::get< 1 >( lastErrorTuple )  / std::get< 1 >( errorTuple ) ) / M_LN2;

				std::cout << "\t\t\t\t\t " << eocL1 << " \t\t \t " << eocL2 << std::endl;

			}

			std::cout << params.numberOfCells << "\t\t\t" << std::get<0> ( errorTuple ) << " \t\t \t " << std::get<1> ( errorTuple ) << std::endl;

			lastErrorTuple = errorTuple;

			// refine
			params.numberOfCells *= 2;
		}



		return 0;

	}

	catch ( Dune::Exception &e )
	{
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch( int e )
	{
		if ( e == 10 )
			std::cerr << "Error: No intersection in cell with his reconstruction." << std::endl;
		else
			throw;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}

}
