#ifndef __DUNE_GRID_REC_VOL_TRACKEVOLVE_HH__
#define __DUNE_GRID_REC_VOL_TRACKEVOLVE_HH__

#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace Dune
{
  namespace VoF
  {
   
    template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
    }

    template < class V >
    void rotate90degreesCounterClockwise( V& v )
    {
	std::swap( v[0], v[1] );
	v[0] = -v[0];
    }

    // Gerade der Form 0 = n*x + p
    template < class V >
    class Line2D
    {

    public:
	Line2D<V> () 
	{
	  n = 0;
	  p = 0;
	}
	Line2D<V> ( const V normal, const V point ) 
	{
	  n = normal;
	  p = n * point;
	  p = -p;
	}
	Line2D<V> ( const Line2D<V>& l ) 
	{
	  n( l.n );
	  p( l.p );
	}  
	Line2D<V> ( const V n2, const double p2 )
	{
	  n = n2;
	  p = p2;
	}
	
	Line2D<V> ( const std::pair<V,V> pair )
	{
	  n = pair.first;
	  n -= pair.second;
	  rotate90degreesCounterClockwise( n );
	  p = n * pair.first;
	  p = -p;
	}

      V n;
      double p;
	
    };

    template < class V >
    const V lineIntersection ( const Line2D<V>& g, const Line2D<V>& l  )
    {
      
      Dune::DynamicMatrix<double> A ( 2, 2 );
      A[0][0] = g.n[0];
      A[0][1] = g.n[1];
      A[1][0] = l.n[0];
      A[1][1] = l.n[1];
      
      Dune::DynamicVector<double> b ( 2 );
      b[0] = -g.p;
      b[1] = -l.p;
      
      
      Dune::DynamicVector<double> x ( 2 );
      
      A.solve( x, b );
      
      
      return x;  
    }


    // konvexes 2D-Polygon
    template < class V >
    class Polygon2D
    {

    public:
	Polygon2D<V> () 
	{}
	
	Polygon2D<V> ( const std::vector<V> vertices )
	{
	  points = vertices;
	}
	
	const V operator[]( const int i ) const
	{
	  return points[i % points.size()]; 
	}

	void addVertex ( const V &vertex, const double TOL = 2e-14 )
	{
	  std::size_t n = points.size();
	  
	  // check, if point is already here
	  for ( std::size_t i = 0; i < n; i++ )
	    if ( (vertex - points[i]).two_norm() < TOL ) return; 

	  if ( n == 0 || n == 1)
	  {
	    points.push_back( vertex );
	    return;
	  }
	  else 
	  {
	    
	    for ( std::size_t i = 0; i < points.size(); i++ )
	    {
	      V normal = points[(i+1)%n] - points[i];
	      rotate90degreesCounterClockwise( normal );
	      
	      if ( normal * ( vertex - points[i] ) < 0 )   
	      {
		points.insert( points.begin() + i + 1, vertex );
		return;
	      }
	    }
	  }
	  
	}
	
	const std::size_t corners() const
	{
	  return points.size();
	}    
	  
	double volume() const
	{
	  double sum = 0;
	  int n = points.size();
	  for ( int i = 0; i < n; i++ )
	  {
	    sum += ( points[i][1] + points[(i+1)%n][1] ) * ( points[i][0] - points[(i+1)%n][0] );
	  }
	  return sum / 2.0;
	}
	
	const bool pointInBorders( const V vertex ) const
	{
	  int n = points.size();
	  
	  if ( n == 0 ) return false;
	  
	  
	  bool inside = true;
	  
	  for ( int i = 0; i < n; ++i )
	    inside = inside && this->SkalarProdTest( vertex, points[i], points[(i+1)%n] );
	  
	  return inside;
	}
	
	std::vector<V> points;
	double vol = 0;
	
    private:
    
	const int SkalarProdTest ( V vertex, V p1, V p2, const double TOL = 1e-14 ) const
	{
	  V n = p2 - p1;
	  rotate90degreesCounterClockwise( n );
	  
	  double skalar = n * ( vertex - p1);
	  
	  return skalar >= 0 || std::abs( skalar ) < TOL;  
	  
	}
	
    };

    template < class V >
    double polygonIntersectionVolume( const Polygon2D<V> polygon1, const Polygon2D<V> polygon2 )
    {
      if ( polygon1.corners() <= 2 || polygon2.corners() <= 2 ) return 0;
      
      
      Polygon2D<V> intersectionPolygon;
      bool lastInside = false, thisInside = false;
      std::vector<std::pair<V,V>> cuttingLines;
      
      
      // Schleife einmal fuer jedes Polygon
      Polygon2D<V> p1 = polygon1;
      Polygon2D<V> p2 = polygon2;
      
      for ( std::size_t p = 0; p <= 1; p++ )
      {
	lastInside = false;
	cuttingLines.clear();
      
	// Ecken in anderem Polygon
	for ( std::size_t i = 0;  i < p1.corners() + 2; ++i )
	{
	  if ( p2.pointInBorders( p1[i] ) )
	  {
	    intersectionPolygon.addVertex( p1[i] );
	    
	    thisInside = true;
	  }
	  
	  if ( lastInside != thisInside && i > 0 )
	  {
	    cuttingLines.push_back( std::pair<V,V>( p1[i-1], p1[i] ) );
	  }
	  
	  lastInside = thisInside;  
	  thisInside = false;
	}
	
	//Schnittpunkte
	V is;
	for ( std::size_t i = 0;  i < p2.corners(); ++i )
	{
	  for ( std::size_t j = 0; j < cuttingLines.size(); j++ )
	  {
	    Line2D<V> cut ( cuttingLines[j] );
	    
	    is = lineIntersection( Line2D<V>( std::pair<V,V>( p2[i], p2[i+1] ) ), cut );
	    
	    if ( (cuttingLines[j].first - is) * (cuttingLines[j].second - is) < 0 
	      && (p2[i] - is) * (p2[i+1] - is) < 0 )
	      intersectionPolygon.addVertex( is );
	  }
	}    
	    
	p1 = polygon2;
	p2 = polygon1;
      }
            
      return intersectionPolygon.volume();
      
    }


    template < class V >
    void insertElementIfNotExists( const V& v, std::vector<V>& list, const double TOL = 1e-14 )    
    {
      bool found = false;
      
      for ( std::size_t i = 0; i < list.size(); ++i )
	if ( ( list[i] - v).two_norm() < TOL ) 
	  found = true;
      
      if ( !found )
	list.push_back( v );
      
    }


    template < class GV, class E, class Geo, class V >
    std::vector<V> lineIntersectionPoints( const GV& gridView, const E& entity, const Geo& geo, const Line2D<V>& g, const double TOL = 1e-14 )
    {
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename GV::Intersection Intersection;
	
      std::vector<V> intersectionPoints;
      
      for ( IntersectionIterator is = gridView.ibegin( entity ); is != gridView.iend( entity ); ++is)
      {
	  const Intersection &edge  = *is;
	  
	  const V c0 = edge.geometry().corner(0);
	  const V c1 = edge.geometry().corner(1);
	  
	  // reconstruction line intersects with edge
	  if ( 
		( isInner( c0, g, TOL ) && isOuter( c1, g, TOL ) )
	     || ( isInner( c1, g, TOL ) && isOuter( c0, g, TOL ) )   )
	  {
	  
	    // build line through edge for intersection
	    const Line2D<V> lineThroughEdge ( edge.centerUnitOuterNormal(), c0 );
	    
	    // add intersection point
	    intersectionPoints.push_back( lineIntersection( g, lineThroughEdge ) );
	  }
	  
	  // reconstruction line intersects with c0
	  else if ( isOnRecLine( c0, g, TOL ) )
	  { 
	    insertElementIfNotExists( c0, intersectionPoints );
	  }
	  //reconstruction line intersects with c1
	  else if ( isOnRecLine( c1, g, TOL ) )
	  {
	    insertElementIfNotExists( c1, intersectionPoints );
	  }
	  
      }
     
      return intersectionPoints;
    }

    
    

    template < class V >
    bool isInner ( const V& vertex, const Line2D<V>& g, const double TOL )
    {
      return vertex * g.n + g.p >= TOL; 
    }
    
    template < class V >
    bool isOuter ( const V& vertex, const Line2D<V>& g, const double TOL )
    {
      return vertex * g.n + g.p <= -TOL; 
    }
    
    template < class V >
    bool isOnRecLine ( const V& vertex, const Line2D<V>& g, const double TOL )
    {
      return std::abs( vertex * g.n + g.p ) < TOL; 
    }
    
    
    
    
    template < class Geo, class V >
    std::vector<V> getInnerVertices( const Geo& geo, const Line2D<V>& g, const double TOL = 1e-8 )
    {
      std::vector<V> innerVertices;
    
      for ( int i = 0; i < geo.corners(); ++i )
	
	if ( isInner( geo.corner(i), g, TOL ) )
	  innerVertices.push_back( geo.corner(i) );
	
	
      return innerVertices;
    }  

    
    
    
    template < class Geo, class V >
    std::vector<V> getOuterVertices( const Geo& geo, const Line2D<V>& g, const double TOL = 1e-8 )
    {
      std::vector<V> outerVertices;
    
      for ( int i = 0; i < geo.corners(); ++i )
	
	if ( isOuter( geo.corner(i), g, TOL ) )
	  outerVertices.push_back( geo.corner(i) );
	
	
      return outerVertices;
    }   


    template < class GV, class E, class Geo, class V >
    double getVolumeFraction ( const GV& gridView, const E& entity, const Geo& geo, const Line2D<V>& g )
    {
      Polygon2D<V> polygonVertices;
      
      for ( auto v : lineIntersectionPoints( gridView, entity, geo, g ) )
	  polygonVertices.addVertex( v );
      
      for ( auto v : getInnerVertices( geo, g ) )
	polygonVertices.addVertex( v );
      
      
      assert( polygonVertices.corners() <= 5 );
      
      return polygonVertices.volume() / geo.volume();
      
    }
 
 
  }
}

#endif


