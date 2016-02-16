#ifndef BINARYWRITER_HH
#define BINARYWRITER_HH

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>

// dune includes
#include <dune/fem/io/io.hh>


// BinaryDataWriter
// ================
class BinaryWriter
{
public:
  using BinaryStream = Dune::Fem::BinaryFileOutStream;

  template < class TP >
  BinaryWriter ( const std::size_t level, const TP& timeProvider ) : level_( level )
  {
    saveTime_ = timeProvider.time();
    path_ = Dune::Fem::Parameter::getValue< typename std::string >( "fem.io.path", "./data/" );
    saveStep_ = Dune::Fem::Parameter::getValue< double >( "fem.io.savestep", 0.1 );
    Dune::Fem::createDirectory ( path_ );
  }

  template < class TimeProviderType >
  const bool willWrite ( const TimeProviderType &timeProvider )
  {
    return timeProvider.time() - saveTime_ >= -0.5 * saveStep_;
  }

  template < class DF >
  const void write ( const DF &uh )
  {
    std::stringstream name;
    name.fill('0');
    name << "vof-fem-" << std::to_string( level_ ) << "-" << std::setw(5) << std::to_string( writeStep_ ) << ".bin";
    BinaryStream binaryStream ( Dune::concatPaths( path_, name.str() ) );
    uh.write( binaryStream );

    saveTime_ += saveStep_;
    writeStep_++;
  }

private:
  std::string path_;
  const std::size_t level_;
  double saveStep_, saveTime_ ;
  std::size_t writeStep_ = 0;
};

#endif
