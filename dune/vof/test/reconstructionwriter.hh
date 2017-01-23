#ifndef RECONSTRUCTIONWRITER_HH
#define RECONSTRUCTIONWRITER_HH

// C++ includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/path.hh>

#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

// local includes
#include "vtu.hh"



// ReconstructionWriter
// --------------------

template < class GridView, class ReconstructionSet, class Flags, class DataParameters, class OutputPolygon >
struct ReconstructionWriter
{
  ReconstructionWriter ( const GridView& gridView, const ReconstructionSet &reconstructionSet, const Flags &flags, const DataParameters &dataParameters )
   : gridView_( gridView ), reconstructionSet_( reconstructionSet ), flags_( flags ), dataParameters_( dataParameters ), count_( dataParameters.startcounter() ) {};

  void write()
  {
    using Dune::VoF::intersect;

    std::vector< OutputPolygon > io;
    for ( const auto& entity : Dune::elements( gridView_ ) )
    {
      if ( !flags_.isMixed( entity ) )
        continue;

      auto polytope = Dune::VoF::makePolytope( entity.geometry() );
      auto it = intersect( polytope, reconstructionSet_[ entity ].boundary() );
      auto intersection = static_cast< typename decltype( it )::Result > ( it );

      std::vector< typename OutputPolygon::Position > is;
      for( std::size_t i = 0; i < intersection.size(); ++i )
        is.push_back( intersection.vertex( i ) );

      if ( !is.empty() )
        io.push_back( OutputPolygon ( is, reconstructionSet_[ entity ].innerNormal() ) );
    }

    VTUWriter< std::vector< OutputPolygon > > vtuwriter( io );

    std::stringstream name;
    name.fill('0');
    name << "s" << std::setw(4) << gridView_.grid().comm().size() << "-p" << std::setw(4) <<  gridView_.grid().comm().rank()
      << "-" << dataParameters_.prefix() << std::setw(5) << count_ << ".vtu";

    vtuwriter.write( Dune::concatPaths( dataParameters_.path(), name.str() ) );

    if ( gridView_.grid().comm().rank() == 0 )
      writeRecPVTUFile();

    count_++;
  }

private:
  void writeRecPVTUFile() const
  {
    const std::size_t size = gridView_.grid().comm().size();

    std::stringstream name;
    name.fill('0');
    name << "s" << std::setw(4) << size << "-" << dataParameters_.prefix() << std::setw(5) << count_ << ".pvtu";

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
      content << "      <Piece  Source=\"s" << std::setw(4) << size << "-p" << std::setw(4) << i << "-" << dataParameters_.prefix()
        << std::setw(5) << count_ << ".vtu\"/>" << std::endl;

    content <<
  R"(  </PUnstructuredGrid>
  </VTKFile>
    )";

    std::fstream f;
    f.open( Dune::concatPaths( dataParameters_.path(), name.str() ), std::ios::out );
    f << content.str();
    f.close();
  }

  const GridView gridView_;
  const ReconstructionSet &reconstructionSet_;
  const Flags &flags_;
  const DataParameters dataParameters_;

  std::size_t count_ = 0;
  //mutable std::vector< OutputPolygon > io;

};

#endif