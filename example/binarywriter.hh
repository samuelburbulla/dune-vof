#ifndef BINARYWRITER_HH
#define BINARYWRITER_HH

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>

// dune-common includes
#include <dune/common/path.hh>
#include <dune/common/parametertree.hh>

// dune-fem includes
#include <dune/fem/io/io.hh>
#include <dune/fem/io/streams/binarystreams.hh>


// BinaryDataWriter
// ================
template< class GridView, class DF >
class BinaryWriter
{
public:
  using BinaryStream = Dune::Fem::BinaryFileOutStream;

  BinaryWriter ( const GridView& gridView, const DF& uh, const Dune::ParameterTree parameters, int level ) : gridView_( gridView ), uh_( uh ), level_( level )
  {
    saveTime_ = std::numeric_limits< double >::min();
    saveStep_ = parameters.get< double >( "io.savestep", 0.1 );
    path_ = parameters.get< std::string >( "io.path", "data" );
    prefix_ = parameters.get< std::string >( "io.prefix", "vof" );
    writeData_ = parameters.get< bool >( "io.writeData", true );
    Dune::Fem::createDirectory ( path_ );
  }

  const bool willWrite ( double time )
  {
    return writeData_ && time - saveTime_ >= -0.5 * saveStep_;
  }

  const void write ( double time, const bool forced = false )
  {

    if ( willWrite( time ) || forced )
    {
      std::stringstream name;
      name.fill('0');
      name << "s" << std::setw(4) << gridView_.comm().size() << "-p" << std::setw(4) << gridView_.comm().rank()
        << "-" << prefix_ << "-" << std::to_string( level_ ) << "-" << std::setw(5) << std::to_string( writeStep_ );

      std::stringstream dfname;
      dfname << name.str() << ".bin";
      BinaryStream binaryStream ( Dune::concatPaths( path_, dfname.str() ) );
      binaryStream << time;
      uh_.write( binaryStream );

      saveTime_ += saveStep_;
      writeStep_++;
    }
  }

private:
  const GridView& gridView_;
  const DF &uh_;
  std::string path_, prefix_;
  const std::size_t level_;
  double saveStep_, saveTime_ ;
  std::size_t writeStep_ = 0;
  bool writeData_;
};

#endif
