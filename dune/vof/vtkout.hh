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
    
    int createFolders( const char* folderName ) {
      
      char path[128];
      sprintf( path, "../Animations/%s", folderName );
      mkdir( path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
      
      char folder[128];
      sprintf( folder, "%s/recValues", path );
      mkdir( folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
      sprintf( folder, "%s/recImages", path );
      mkdir( folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
      sprintf( folder, "%s/vtk", path );
      mkdir( folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

      std::ofstream massenFile;
      char name2[128];
      sprintf( name2, "../Animations/%s/massenerhaltung", folderName );
      massenFile.open ( name2 );
      massenFile << "";
      massenFile.close();
      
      return 0;
    }


    void createMeshforPlot( const int numberOfCells, const char* folderName )
    {
      std::ofstream myfile;
      char name[128];
      sprintf( name, "../Animations/%s/recValues/mesh", folderName );
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

    template < class R >
    void writeReconstructionToGnuplot( const int k, const R& reconstruction, const char* folderName )
    {  
      std::ofstream myfile;
      char name[128];
      sprintf( name, "../Animations/%s/recValues/recVal_%d", folderName, k );

      myfile.open ( name );
      
      for ( std::size_t i = 0; i < reconstruction.size(); i++ )
      {
	myfile << reconstruction[i][0][0] << " " << reconstruction[i][0][1] << std::endl;
	myfile << reconstruction[i][1][0] << " " << reconstruction[i][1][1] << std::endl;
	
	myfile << std::endl;
      }
      myfile.close();
	
	
      
      FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
      
      fprintf( gnuplotPipe, "set terminal png size 2000,2000 enhanced font \"Helvetica,30\"\n" );
      fprintf( gnuplotPipe, "set output \"../Animations/%s/recImages/recImage_%d.png\"\n", folderName, k );  
      fprintf( gnuplotPipe, "set xrange[0:1]\n" );
      fprintf( gnuplotPipe, "set yrange[0:1]\n" );
      fprintf( gnuplotPipe, "plot \"%s\" with lines ls 1 linecolor 1 linewidth 2, \"../Animations/%s/recValues/mesh\" with lines ls 1 lt rgb \"#CCCCCC\" \n", name, folderName );
      
      fclose(gnuplotPipe);
      


    }

    template < class G, class V, class R >
    void vtkout( const G& grid, const V& c, const char * name, const char* folderName, const int numberOfCells, const R& reconstruction, int k, double time = 0.0 )
    {
      
      
      // write VTK file
      Dune::VTKWriter <typename G::LeafGridView> vtkwriter( grid.leafGridView() );
      char fname[128];
      char sername[128];
      sprintf( fname, "../Animations/%s/vtk/%s-%05d", folderName, name, k );


      sprintf( sername, "%s.series", name );
      vtkwriter.addCellData ( c, "celldata" );
      
      vtkwriter.write( fname, Dune::VTK::ascii );
      
      
      writeReconstructionToGnuplot( k, reconstruction, folderName );
      
      //Massenerhaltung
      std::ofstream massenFile;
      char name2[128];
      sprintf( name2, "../Animations/%s/massenerhaltung", folderName );
      massenFile.open ( name2, std::fstream::app );
      
      double sum = 0;
      for ( std::size_t i = 0; i < c.size(); ++i )
	sum += c[i] / ( numberOfCells * numberOfCells );
      
      massenFile << time << " " << std::abs( sum - M_PI * 0.15 * 0.15 ) << std::endl;
      
      massenFile.close();
      
    }

  }
}

# endif
