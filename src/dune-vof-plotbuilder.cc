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



int plotbuilder( int max_k, const std::string &folderPath )
{
  
  for ( int k = 0; k < max_k; ++k )  
  {
  
    std::string meshData = folderPath + "values/mesh";
    std::string reconstructionData = folderPath + "values/val_" + std::to_string( k );
    std::string mixedCellsData = folderPath + "values/mixed_" + std::to_string(  k );
    std::string activeCellsData = folderPath + "values/active_" + std::to_string( k );
    
    std::string imagePath = folderPath + "images/plot_" + std::to_string( k ) + ".png";
    
    
    // plot with gnuplot
    std::ofstream gnuplotPipe;
    gnuplotPipe.open( "gnuplot -persistent" );
    
    gnuplotPipe << "set terminal png size 2000,2000 enhanced font \"Helvetica,30\"\n";
    gnuplotPipe << "set output \"" << imagePath << "\"\n";  
    gnuplotPipe << "set xrange[0:1]\n";
    gnuplotPipe << "set yrange[0:1]\n";
    gnuplotPipe << "plot \"" << reconstructionData << "\" with lines ls 1 linecolor 1 linewidth 2,"
		<< "\"" << mixedCellsData << "\" w p ls 1 pt 5,"
		<< "\"" << activeCellsData << "\" w p ls 2 pt 5," 
		<< "\"" << meshData << "\" with lines ls 1 lt rgb \"#CCCCCC\" \n";

    gnuplotPipe.close();
  
  }
  return 0;

}
