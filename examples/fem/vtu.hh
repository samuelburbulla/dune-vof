#ifndef VTU_HH
#define VTU_HH

#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

#include <dune/common/fvector.hh>


template< class Data >
struct VTUDataType;

template<>
struct VTUDataType< float >
{
  static std::string type () { return "Float32"; }
  static int components () { return 1; }

  static std::string toString ( const float &value )
  {
    std::ostringstream ret;
    ret << std::scientific << std::setprecision( 8 ) << value;
    return ret.str();
  }
};

template<>
struct VTUDataType< double >
{
  static std::string type () { return "Float64"; }
  static int components () { return 1; }

  static std::string toString ( const double &value )
  {
    std::ostringstream ret;
    ret << std::scientific << std::setprecision( 16 ) << value;
    return ret.str();
  }
};

template<>
struct VTUDataType< std::int32_t >
{
  static std::string type () { return "Int32"; }
  static int components () { return 1; }
  static std::string toString ( const std::int32_t &value ) { return std::to_string( value ); }
};

template<>
struct VTUDataType< std::uint8_t >
{
  static std::string type () { return "UInt8"; }
  static int components () { return 1; }
  static std::string toString ( const std::uint8_t &value ) { return std::to_string( value ); }
};

template< >
struct VTUDataType< Dune::FieldVector< double, 2 > >
{
  typedef VTUDataType< double > VTUDouble;

  static std::string type () { return VTUDouble::type(); }
  static int components () { return 3; }

  static std::string toString ( const Dune::FieldVector< double, 2 > &value )
  {
    return VTUDouble::toString( value[ 0 ] ) + " " + VTUDouble::toString( value[ 1 ] ) + " " + VTUDouble::toString( double( 0 ) );
  }
};



// VTUWriter
// ---------

template< class PolygonVector >
class VTUWriter
{
  typedef typename PolygonVector::value_type Polygon;
  typedef typename Polygon::Position Position;

public:
  VTUWriter ( const PolygonVector &v ) : v_( v )
  {
    resize();
  }

  void write ( const std::string &name ) const
  {
    resize();
    std::ofstream vtu( name );
    vtu << "<?xml version=\"1.0\"?>" << std::endl;
    vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << std::endl;
    vtu << "  <UnstructuredGrid>" << std::endl;
    vtu << "    <Piece NumberOfPoints=\"" << overallSize_ << "\" "
        << "NumberOfCells=\"" << v_.size() << "\">" << std::endl;

    vtu << "      <CellData Normals=\"normals\">" << std::endl;
    writeNormals( vtu );
    vtu << "      </CellData>" << std::endl;

    vtu << "      <Points>" << std::endl;
    writeCoordinates( vtu );
    vtu << "      </Points>" << std::endl;

    vtu << "      <Cells>" << std::endl;
    writeConnectivity( vtu );
    writeOffsets( vtu );
    writeTypes( vtu );
    vtu << "      </Cells>" << std::endl;

    vtu << "    </Piece>" << std::endl;
    vtu << "  </UnstructuredGrid>" << std::endl;
    vtu << "</VTKFile>" << std::endl;
  }

private:

  void resize () const
  {
    overallSize_ = 0;
    for( const auto &entry : v_ )
      overallSize_ += entry.size();
  }

  void writeCoordinates ( std::ostream &vtu ) const
  {
    std::vector< Position > points( overallSize_ );
    std::size_t index = 0;

    for( const Polygon &pol : v_ )
      for( const Position &pos : pol )
        points[ index++ ] = pos;

    writeDataArray( vtu, "Coordinates", points );
  }

  void writeConnectivity ( std::ostream &vtu ) const
  {
    std::vector< int > connectivity( overallSize_ );
    for( std::size_t i = 0; i < overallSize_; ++i )
      connectivity[ i ] = i;

    writeDataArray( vtu, "connectivity", connectivity );
  }

  void writeOffsets ( std::ostream &vtu ) const
  {
    const std::size_t size = v_.size();
    std::vector< std::int32_t > offsets( size );
    offsets[ 0 ] = v_[ 0 ].size();
    for( std::size_t i = 1; i < size; ++i )
      offsets[ i ] = offsets[ i - 1 ] + v_[ i ].size();
    writeDataArray( vtu, "offsets", offsets );
  }

  void writeNormals( std::ostream &vtu ) const
  {
    std::vector< Position > normals( overallSize_ );
    for( std::size_t i = 0; i < overallSize_; ++i  )
      normals[ i ] = v_[ i ].normal();

    writeDataArray( vtu, "normals", normals );
  }

  void writeTypes ( std::ostream &vtu ) const
  {
    std::vector< std::uint8_t > types( v_.size(), 4 );
    writeDataArray( vtu, "types", types );
  }

  template< class Data >
  static void writeDataArray ( std::ostream &vtu, const std::string &name, const std::vector< Data > &data )
  {
    vtu << "        <DataArray type=\"" << VTUDataType< Data >::type() << "\" "
        << "NumberOfComponents=\"" << VTUDataType< Data >::components() << "\" "
        << "Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    std::size_t linelength = 0;
    for( const Data &value : data )
    {
      const std::string s = VTUDataType< Data >::toString( value );
      const std::size_t sz = s.size();
      if( linelength + sz > 79 )
      {
        vtu << std::endl;
        linelength = 0;
      }
      std::string indent( linelength < 10 ? "          " : " " );
      vtu << indent << s;
      linelength += (sz + indent.size());
    }
    if( linelength > 0 )
      vtu << std::endl;
    vtu << "        </DataArray>" << std::endl;
  }

  const PolygonVector &v_;
  mutable std::size_t overallSize_;
};

#endif // #ifndef VTU_HH
