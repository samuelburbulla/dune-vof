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

#include <dune/vof/flagCells.hh>
#include <dune/vof/errors.hh>

#include <dune/vof/vof.hh>
#include <dune/vof/vtkout.hh>
#include <dune/vof/domain.hh>
#include <dune/vof/initialize.hh>


#include <dune/vof/evolve.hh>
#include <dune/vof/reconstruct.hh>


    double eps = 1e-3;


    template < class G, class C, class R, class D >
    void timeloop ( const G& grid, C& concentration, R& reconstruction, D& domain, 
		    std::vector<bool>& cellIsMixed, std::vector<bool>& cellIsActive, double tEnd, double dt, 
		    const int numberOfCells, const double eps, const char* folderName )
    {
      double t = 0;
      int k = 0;
      const double saveInterval = 0.01;
      double saveStep = 0.01;
      int saveNumber = 1;
      
      
      Dune::VoF::flagCells( grid.leafGridView(), concentration, reconstruction, domain, numberOfCells, cellIsMixed, cellIsActive, eps );
      Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, eps ); 
      Dune::VoF::vtkout( grid, concentration, "concentration", folderName, numberOfCells, reconstruction, 0, cellIsMixed, cellIsActive, t );
      
            
      
      while ( t < tEnd )
      {
	++k;   
	
	Dune::VoF::flagCells( grid.leafGridView(), concentration, reconstruction, domain, numberOfCells, cellIsMixed, cellIsActive, eps );
	
	Dune::VoF::reconstruct( grid, concentration, reconstruction, cellIsMixed, domain, eps ); 
	  
	Dune::VoF::evolve( grid, concentration, reconstruction, domain, numberOfCells, t, dt, eps, cellIsMixed, cellIsActive );
	
	t += dt;
	
	if ( t >= saveStep )
	{
	  Dune::VoF::vtkout( grid, concentration, "concentration", folderName, numberOfCells, reconstruction, saveNumber, cellIsMixed, cellIsActive, t );
	  saveStep += saveInterval;
	  ++saveNumber;
	}
	
	std::cout << "s=" << grid.size(0) << " k=" << k << " t=" << t << " dt=" << dt << " saved=" << saveNumber-1 << std::endl;
      }
    }




    int main(int argc, char** argv)
    {
      Dune::MPIHelper::instance( argc, argv );
	
      try {
	
	//mesh size
	int numberOfCells = 32;
	// folder to store the plots
	char folderName[128] = "default";
	
	Dune::VoF::handleInput ( argc, argv, numberOfCells, folderName );
      
	Dune::VoF::createFolders( folderName );
	
	
	// type declarations
	const int dim = 2;
	typedef Dune::YaspGrid<dim> GridType;
	typedef typename Dune::FieldVector<double,dim> fvector;
	
	
	// build Grid
	fvector upper( 1.0 );
	Dune::array<int,dim> noc;
	std::fill( noc.begin(), noc.end(), numberOfCells );
	GridType grid( upper, noc );

	
	// build domain references for each cell
	Dune::VoF::DomainOfCells<GridType> domain ( grid );
	
	// build an PNG-File from all edges to plot the reconstruction
	Dune::VoF::createMeshforPlot( numberOfCells, folderName );
	  
	
	
	int n = grid.leafGridView().size(0);
	
	// allocate and initialize a vector for the concentration
	std::vector<double> concentration( n );
	Dune::VoF::initialize( grid, concentration );
      
	
	// allocate and initialize a vector for the reconstruction
	std::vector< std::array<fvector,3> > reconstruction( n );
	  
    
	
	// allocate vectors for flagging cells
	std::vector<bool> cellIsMixed ( n );
	std::vector<bool> cellIsActive ( n );
	  
	
	
	// calculate dt
	double dt = 0.01 * ( 1.0 / numberOfCells ) / Dune::VoF::psiMax();
	
	
	
	
	
	
	// concentration on start to check the error
	std::vector<double> concentrationStart ( concentration );
	
	    
	// start time integration
	timeloop( grid, concentration, reconstruction, domain, cellIsMixed, cellIsActive, 2 * M_PI, dt, numberOfCells, eps, folderName ); 
	
	
	//concentration on the end to check the error
	std::vector<double> concentrationEnd ( concentration );
	
	
	
	
	// calculate the L1-error
	double L1 = Dune::VoF::l1error( grid, concentrationStart, concentrationEnd );    
	std::cout << "L1-Fehler: " << L1 << std::endl; 
	
	// calculate the L2-error
	double L2 = Dune::VoF::l2error( grid, concentrationStart, concentrationEnd );    
	std::cout << "L2-Fehler: " << L2 << std::endl;   
	    
	
	return 0;
      }
      catch (Dune::Exception &e){
	std::cerr << "Dune reported error: " << e << std::endl;
      }
      catch (...){
	std::cerr << "Unknown exception thrown!" << std::endl;
      }
    }
