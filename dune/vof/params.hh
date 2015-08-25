#ifndef DUNE_VOF_PARAMS_HH
#define DUNE_VOF_PARAMS_HH

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
    }


  } // end of namespace VoF
} // end of namespace Dune

#endif // DUNE_VOF_PARAMS_HH
