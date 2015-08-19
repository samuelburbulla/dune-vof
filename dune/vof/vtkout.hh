# ifndef __DUNE_REC_VOL_TRACK_VTKOUT_HH__
# define __DUNE_REC_VOL_TRACK_VTKOUT_HH__

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>

namespace Dune
{
  namespace VoF
  { 


    void createMeshforPlot( const int numberOfCells, const std::string &folderName )
    {
      std::ofstream myfile;
      std::string name ( folderName + "values/mesh" );
      myfile.open ( name );
      
      double dx = 1.0 / numberOfCells;
      
      // vertical lines
      for ( int i = 1; i < numberOfCells; i++ )
      {
	myfile << i * dx << " " << 0 << std::endl;
	myfile << i * dx << " " << 1 << std::endl;
	
	myfile << std::endl;
      }
      
      // horizontal lines
      for ( int i = 1; i < numberOfCells; i++ )
      {
	myfile << 0 << " " << i * dx << std::endl;
	myfile << 1 << " " << i * dx << std::endl;
	
	myfile << std::endl;
      }
      
      myfile.close();
    }

    template < class R, class G >
    void writeReconstructionToGnuplotFiles( const std::string s, const R& reconstruction, 
				       const std::vector<bool>& cellIsMixed, 
				       const std::vector<bool>& cellIsActive,
				       const G& grid,
				       const int numberOfCells, const std::string &folderName )
    {  
      
      std::ofstream valFile;
      valFile.open ( folderName + "values/val_" + s );
      
      std::ofstream mixedFile;
      mixedFile.open ( folderName + "values/mixed_" + s );
      
      std::ofstream activeFile;
      activeFile.open ( folderName + "values/active_" + s );
    
      
      for ( auto entity : elements( grid.leafGridView() ) )
      {
	int i = grid.leafGridView().indexSet().index( entity );
	
	valFile << reconstruction[i][0][0] << " " << reconstruction[i][0][1] << std::endl;
	valFile << reconstruction[i][1][0] << " " << reconstruction[i][1][1] << std::endl;
	valFile << std::endl;
	
	auto geo = entity.geometry();
	
	if ( cellIsMixed[i] )
	{
	  mixedFile << geo.center()[0] << " " << geo.center()[1] << std::endl;
	}
	if ( cellIsActive[i] )
	{
	  activeFile << geo.center()[0] << " " << geo.center()[1] << std::endl;
	}
	
	
      }
      valFile.close();
      mixedFile.close();
      activeFile.close();      


    }
    
    template < class C >
    void writeSumOfMassToFile( const C& c, const double time, int numberOfCells, const std::string &folderName )
    {
      std::ofstream massenFile;
      massenFile.open ( folderName + "massenerhaltung", std::fstream::app );
      
      double sum = 0;
      for ( std::size_t i = 0; i < c.size(); ++i )
	sum += c[i] / ( numberOfCells * numberOfCells );
      
      double referenceMass = M_PI * 0.15 * 0.15;
      
      massenFile << time << " " << std::abs( sum - referenceMass ) << std::endl;
      
      massenFile.close();
    }

    
       
    
    template < class G, class V, class R >
    void vtkout( const G& grid, const V& c, const std::string &name, const std::string &folderName, 
		 const int numberOfCells, const R& reconstruction, int k,
		 const std::vector<bool>& cellIsMixed, 
		 const std::vector<bool>& cellIsActive,
		 double time = 0.0 )
    {
      
      std::stringstream ss;
      ss << std::setw(6) << std::setfill('0') << k;
      std::string s = ss.str();
      
      // write VTK file
      Dune::VTKWriter <typename G::LeafGridView> vtkwriter( grid.leafGridView() );
      std::string fname ( folderName + "vtk/concentration-" + s );
      std::string sername ( "plot.series" );
      vtkwriter.addCellData ( c, "celldata" );
      vtkwriter.addCellData ( cellIsMixed, "cellmixed" );
      vtkwriter.addCellData ( cellIsActive, "cellactive" );
      
      vtkwriter.write( fname, Dune::VTK::ascii );
      
      
      
      writeReconstructionToGnuplotFiles( s, reconstruction, cellIsMixed, cellIsActive, grid, numberOfCells, folderName );
      
      writeSumOfMassToFile( c, time, numberOfCells, folderName );

      
    }

  }
}

# endif
