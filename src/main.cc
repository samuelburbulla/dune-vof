#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <fstream>
#include <vector>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <cmath>
#include <functional>
#include <stdexcept>


#include <dune/vof/flagCells.hh>
#include <dune/vof/errors.hh>

#include <dune/vof/main.hh>
#include <dune/vof/vtkout.hh>
#include <dune/vof/domain.hh>
#include <dune/vof/initialize.hh>


#include <dune/vof/evolve.hh>
#include <dune/vof/reconstruct.hh>


class Parameters 
{
public:
  	Parameters () {
  		folderName = folderPath + nameOfSeries + "/";
  	};

  	// concetration eps for mixed cells
  	double eps = 1e-3;

  	// number of cells for the cartesian grid in one direction
  	int numberOfCells = 32;

  	// params for timeloop
	double dtAlpha = 0.1;
	double tEnd = 10.0;

	// params for saving data
	double saveInterval = 0.01;

	std::string folderPath = "results/";
	std::string nameOfSeries = "default";
	std::string folderName;
	
};


template < class Grid >
void algorithm ( const Grid& grid, const Parameters &params )
{
	typedef Dune::FieldVector<double,2> fvector;

	// build domain references for each cell
	Dune::VoF::DomainOfCells<Grid> domain ( grid );

	auto gridView = grid.leafGridView();

	int n = grid.leafGridView().size(0);

	// allocate and initialize vectors for data representation
	std::vector<double> concentration( n );
	std::vector< std::array<fvector,3> > reconstruction( n );
	std::vector<bool> cellIsMixed ( n );
	std::vector<bool> cellIsActive ( n );


	Dune::VoF::initialize( grid, concentration, Dune::VoF::f0<fvector> );
	 



	// calculate dt
	double dt = params.dtAlpha * ( 1.0 / params.numberOfCells ) / Dune::VoF::psiMax();


	double t = 0;
	int k = 0;

	  
	Dune::VoF::flagCells( gridView, concentration, reconstruction, domain, params.numberOfCells, cellIsMixed, cellIsActive, params.eps );
	Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, params.eps ); 
	Dune::VoF::vtkout( grid, concentration, "concentration", params.folderName, params.numberOfCells, reconstruction, 0, cellIsMixed, cellIsActive, t );
	  
	int saveNumber = 1;  
	double saveStep = 0.01;
	  
	while ( t < params.tEnd )
	{
		++k;   


		Dune::VoF::clearReconstruction( reconstruction );

		Dune::VoF::flagCells( grid.leafGridView(), concentration, reconstruction, domain, params.numberOfCells, cellIsMixed, cellIsActive, params.eps );
		Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, params.eps ); 
		Dune::VoF::evolve( grid, concentration, reconstruction, domain, params.numberOfCells, t, dt, params.eps, cellIsMixed, cellIsActive );

		t += dt;

		if ( t >= saveStep )
		{
		  Dune::VoF::vtkout( grid, concentration, "concentration", params.folderName, params.numberOfCells, reconstruction, saveNumber, cellIsMixed, cellIsActive, t );
		  saveStep += params.saveInterval;
		  ++saveNumber;
		}


		



		std::cerr << "s=" << grid.size(0) << " k=" << k << " t=" << t << " dt=" << dt << " saved=" << saveNumber-1 << std::endl;
	}

	std::vector<double> concentrationEnd( n, 0 );

	auto ft = std::bind( Dune::VoF::f<fvector>, std::placeholders::_1, t);
	Dune::VoF::initialize( grid, concentrationEnd, ft );

	double L1 = Dune::VoF::l1error( grid, concentration, concentrationEnd );    
	double L2 = Dune::VoF::l2error( grid, concentration, concentrationEnd );   
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
		  



		// build Grid
		fvector upper( 1.0 );
		Dune::array<int,dim> noc;
		std::fill( noc.begin(), noc.end(), params.numberOfCells );
		GridType grid( upper, noc );

	  

		    
		// start time integration
		algorithm( grid,  params ); 

		    

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
