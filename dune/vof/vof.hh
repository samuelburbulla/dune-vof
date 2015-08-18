#ifndef DUNE_VOF_HH
#define DUNE_VOF_HH

namespace Dune
{
  namespace VoF
  {
	
    // handle commandline input
    void handleInput ( int argc, char** argv, int numberOfCells, std::string &nameOfSeries )
    {
      // given mesh size
      if (argc > 1)
      {
	numberOfCells = atoi( argv[1] );
     
	nameOfSeries = numberOfCells + "x" + numberOfCells;
      }
	
    }

  }
}

#endif // DUNE_VOF_HH
