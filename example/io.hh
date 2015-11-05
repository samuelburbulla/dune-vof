# ifndef IO_HH
# define IO_HH

// C++ inludes
#include <fstream>
#include <iostream>
#include <string>

// C includes
#include <cstdio>
#include <sys/stat.h>
#include <dirent.h>


    // Copied from Dune-Fem (dune/fem/io/io.cc)
    // ----------------------------------------

    bool createDirectory ( const std::string &inName )
    {
      std::string name = inName;

      // strip of last character if it is a '/'
      if( name[ name.size() - 1 ] == '/' )
        name = name.substr( 0, name.size() -1 );

      // try to open directory (returns null pointer, if path does not exist)
      DIR *dir = opendir( name.c_str() );
      if( dir != 0 )
      {
        if( closedir( dir ) < 0 )
          std::cerr << "Error: Could not close directory." << std::endl;
        return true;
      }

      // try to create the father directory
      size_t pos = name.rfind( '/' );
      if( pos != std::string::npos )
      {
        const std::string father = name.substr( 0, pos );
        if( !createDirectory( father ) )
          return false;
      }

      // try to create the new child directory
      mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
      return (mkdir( name.c_str(), mode ) >= 0);
    }

# endif // #ifndef IO_HH
