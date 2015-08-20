# ifndef __DUNE_REC_VOL_TRACK_VTKOUT_HH__
# define __DUNE_REC_VOL_TRACK_VTKOUT_HH__

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

namespace Dune
{
  namespace VoF
  {


    void createMeshforPlot ( const int numberOfCells, const std::string &folderName )
    {
      std::ofstream myfile;
      std::string name( folderName + "values/mesh" );
      myfile.open( name );

      double dx = 1.0 / numberOfCells;

      // vertical lines
      for( int i = 1; i < numberOfCells; i++ )
      {
        myfile << i * dx << " " << 0 << std::endl;
        myfile << i * dx << " " << 1 << std::endl;

        myfile << std::endl;
      }

      // horizontal lines
      for( int i = 1; i < numberOfCells; i++ )
      {
        myfile << 0 << " " << i * dx << std::endl;
        myfile << 1 << " " << i * dx << std::endl;

        myfile << std::endl;
      }

      myfile.close();
    }

    template< class R, class G >
    void writeReconstructionToGnuplotFiles ( const std::string s, const R &reconstruction,
                                             const std::vector< bool > &cellIsMixed,
                                             const std::vector< bool > &cellIsActive,
                                             const G &grid,
                                             const int numberOfCells, const std::string &folderName )
    {

      std::ofstream valFile;
      valFile.open( folderName + "values/val_" + s );

      std::ofstream mixedFile;
      mixedFile.open( folderName + "values/mixed_" + s );

      std::ofstream activeFile;
      activeFile.open( folderName + "values/active_" + s );


      for( auto entity : elements( grid.leafGridView() ) )
      {
        int i = grid.leafGridView().indexSet().index( entity );

        valFile << reconstruction[ i ][ 0 ][ 0 ] << " " << reconstruction[ i ][ 0 ][ 1 ] << std::endl;
        valFile << reconstruction[ i ][ 1 ][ 0 ] << " " << reconstruction[ i ][ 1 ][ 1 ] << std::endl;
        valFile << std::endl;

        auto geo = entity.geometry();

        if( cellIsMixed[ i ] )
          mixedFile << geo.center()[ 0 ] << " " << geo.center()[ 1 ] << std::endl;
        if( cellIsActive[ i ] )
          activeFile << geo.center()[ 0 ] << " " << geo.center()[ 1 ] << std::endl;


      }
      valFile.close();
      mixedFile.close();
      activeFile.close();


    }

  } // end of namespace VoF
} // end of namespace Dune

# endif
