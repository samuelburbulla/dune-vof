#ifndef DUNE_MAIN_HH
#define DUNE_MAIN_HH

namespace Dune
{
  namespace VoF
  {

    // handle commandline input
    template< class Parameters >
    void handleInput ( int argc, char **argv, Parameters &params )
    {
      // given mesh size
      if( argc > 1 )
        params.numberOfCells = atoi( argv[ 1 ] );

      //params.nameOfSeries = params.numberOfCells + "x" + params.numberOfCells;

      //params.folderName = params.folderPath + params.nameOfSeries + "/";

    }


  } // end of namespace VoF
} // end of namespace Dune

#endif // DUNE_VOF_HH
