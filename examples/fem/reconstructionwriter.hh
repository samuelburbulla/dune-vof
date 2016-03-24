#ifndef RECONSTRUCTIONWRITER_HH
#define RECONSTRUCTIONWRITER_HH

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/vof/geometry/intersection.hh>
#include <dune/vof/geometry/polygon.hh>
#include <dune/vof/geometry/polytope.hh>

// local includes
#include "vtu.hh"



// ReconstructionWriter
// --------------------

template < class GridView, class ReconstructionSet, class Polygon >
struct ReconstructionWriter
{
  ReconstructionWriter ( const GridView& gridView, const std::size_t level ) : gridView_( gridView ), level_ ( level ) {};

  void write( const ReconstructionSet &reconstructionSet )
  {
    using Dune::VoF::make_polygon;
    using Dune::VoF::intersect;

    std::vector< Polygon > io;
    for ( const auto& entity : Dune::elements( gridView_ ) )
    {
      Dune::VoF::Line< typename Polygon::Position > intersection = intersect( make_polygon( entity.geometry() ), reconstructionSet[ entity ].boundary() );

      std::vector< typename Polygon::Position > is;
      for( int i = 0; i < intersection.size(); ++i )
        is.push_back( intersection.vertex( i ) );

      if ( !is.empty() )
        io.push_back( Polygon ( is, reconstructionSet[ entity ].innerNormal() ) );
    }

    VTUWriter< std::vector< Polygon > > vtuwriter( io );

    std::stringstream path;
    path << "./data/";

    std::stringstream name;
    name.fill('0');
    name << "s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-p" << std::setw(4) <<  Dune::Fem::MPIManager::rank()
      << "-vof-rec-" << std::to_string( level_ ) << "-" << std::setw(5) << count_ << ".vtu";

    vtuwriter.write( Dune::concatPaths( path.str(), name.str() ) );

    if ( Dune::Fem::MPIManager::rank() == 0 )
      writeRecPVTUFile( path );

    count_++;
  }

  std::size_t& count () { return count_; }

private:

  void writeRecPVTUFile( const std::stringstream& path ) const
  {
    const std::size_t size = Dune::Fem::MPIManager::size();

    std::stringstream name;
    name.fill('0');
    name << "s" << std::setw(4) << size << "-vof-rec-" << std::to_string( level_ ) << "-" << std::setw(6) << count_ << ".pvtu";

    std::stringstream content;
    content.fill('0');
    content <<
  R"(<?xml version="1.0"?>
  <VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">
    <PUnstructuredGrid GhostLevel="0">
      <PCellData Normals="normals">
        <PDataArray type="Float64" NumberOfComponents="3" Name="normals" format="ascii"/>
      </PCellData>
      <PPoints>
        <PDataArray type="Float64" NumberOfComponents="3" Name="Coordinates" format="ascii"/>
      </PPoints>
      <PCells>
        <PDataArray type="Int32" NumberOfComponents="1" Name="connectivity" format="ascii"/>
        <PDataArray type="Int32" NumberOfComponents="1" Name="offsets" format="ascii"/>
        <PDataArray type="UInt8" NumberOfComponents="1" Name="types" format="ascii"/>
      </PCells>
)";

    for ( std::size_t i = 0; i < size; ++i )
      content << "      <Piece  Source=\"s" << std::setw(4) << Dune::Fem::MPIManager::size() << "-p" << std::setw(4) << i << "-vof-rec-" << level_ << "-"
        << std::setw(5) << count_ << ".vtu\"/>" << std::endl;

    content <<
  R"(  </PUnstructuredGrid>
  </VTKFile>
    )";

    std::fstream f;
    f.open( Dune::concatPaths( path.str(), name.str() ), std::ios::out );
    f << content.str();
    f.close();
  }

  const GridView gridView_;
  const std::size_t level_;
  std::size_t count_ = 0;
  mutable std::vector< Polygon > io;

};

#endif
