#ifndef DUNE_VOF_HH
#define DUNE_VOF_HH

namespace Dune
{
  namespace VoF
  {
	
    // handle commandline input
    void handleInput ( int argc, char** argv, int numberOfCells, char* folderName )
    {
      // given mesh size
      if (argc > 1)
	numberOfCells = atoi( argv[1] );
	

	// arguments
	if (argc > 2)
	{
	  for( int i = 2;  i < argc; ++i)
	  {
	  
	    if ( strcmp( argv[i], "folder") == 0 )
	    {
	      if ( i+1 < argc )
		sprintf( folderName , "%s", argv[i+1] );
	    }
	  }
	}
    }

  }
}

#endif // DUNE_VOF_HH
