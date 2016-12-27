#ifndef BINARYWRITER_HH
#define BINARYWRITER_HH

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>

// dune-common includes
#include <dune/common/path.hh>

// dune-fem includes
#include <dune/fem/io/io.hh>
#include <dune/fem/io/streams/binarystreams.hh>


// BinaryDataWriter
// ================
class BinaryWriter
{
public:
  using BinaryStream = Dune::Fem::BinaryFileOutStream;

  BinaryWriter ( const Dune::ParameterTree parameters, int level, double time ) : level_( level )
  {
    saveTime_ = time;
    saveStep_ = parameters.get< double >( "io.savestep", 0.1 );
    path_ = parameters.get< std::string >( "io.path", "data" );
    prefix_ = parameters.get< std::string >( "io.prefix", "vof" );
    Dune::Fem::createDirectory ( path_ );
  }

  const bool willWrite ( double time )
  {
    return time - saveTime_ >= -0.5 * saveStep_;
  }

  template < class Grid, class DF >
  const void write ( const Grid& grid, const DF &uh, double time, const bool forced = false )
  {

    if ( willWrite( time ) || forced )
    {
      std::stringstream name;
      name.fill('0');
      name << "s" << std::setw(4) << grid.comm().size() << "-p" << std::setw(4) << grid.comm().rank()
        << "-" << prefix_ << "-" << std::to_string( level_ ) << "-" << std::setw(5) << std::to_string( writeStep_ );

      std::stringstream dfname;
      dfname << name.str() << ".bin";
      BinaryStream binaryStream ( Dune::concatPaths( path_, dfname.str() ) );
      binaryStream << time;
      uh.write( binaryStream );

      saveTime_ += saveStep_;
      writeStep_++;
    }
  }

private:
  std::string path_, prefix_;
  const std::size_t level_;
  double saveStep_, saveTime_ ;
  std::size_t writeStep_ = 0;
};

#endif
