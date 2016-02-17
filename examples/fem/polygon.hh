#ifndef POLYGON_HH
#define POLYGON_HH

#include <initializer_list>

template< class PointType >
class Polygon
{
  typedef std::vector< PointType > Container;
public:
  typedef PointType Position;
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;

  Polygon () {}

  Polygon ( std::vector< PointType > const &l, const PointType &normal ) : data_( l.begin(), l.end() ), normal_( normal ) {}

  Position &operator[] ( std::size_t i ) { return data_[ i ]; }
  const Position &operator[] ( std::size_t i ) const { return data_[ i ]; }

  Position& normal () { return normal_; }
  const Position& normal () const { return normal_; }

  iterator begin () { return data_.begin(); }
  const_iterator begin () const { return data_.begin(); }

  iterator end () { return data_.end(); }
  const_iterator end () const { return data_.end(); }

  void resize ( std::size_t newSize ) { data_.resize( newSize ); }

  std::size_t size () const { return data_.size(); }

  bool operator== ( const Polygon &other ) { return data_ == other.data_; }
  bool operator!= ( const Polygon &other ) { return data_ != other.data_; }

private:
  Container data_;
  Position normal_;
};


#endif
