#ifndef __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__
#define __DUNE_GRID_REC_VOL_RECONSTRUCT_HH__

#include <stdio.h>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include "geometricaltoolbox.hh"


namespace Dune
{
  namespace VoF
  {

    template < class G, class C, class E, class V, class D >
    void guessNormal( const G &grid, const C &c, const E &entity, const D& domain, Line2D<V>& g );

    template < class GV, class E, class V >
    double computeInterfaceLinePosition( const GV& gridView, const E& entity, const V& n, double concentration );

    template < class GV, class E, class V >
    V interfaceLineCentroid( const GV& gridView, const E& entity, const V& normal, const double pEntity );

    double brentsMethod ( double a, double b, double (*f) ( double p ) );


      
    template < class G, class C, class R, class F, class D >
    void reconstruct ( const G& grid, C& concentration, R& reconstruction, const F& cellIsMixed, const D& domain, const double eps )
    {
      const int dimworld = G::dimensionworld;
      typedef typename Dune::FieldVector<double, dimworld> fvector;
      typedef typename G::LeafGridView GridView;
      typedef typename GridView::template Codim<0>::Iterator LeafIterator;
      typedef typename LeafIterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;  

      // clear reconstruction
      fvector null; null = 0;
      for ( std::size_t i = 0; i < reconstruction.size(); ++i )
	reconstruction[i] = std::array<fvector,3> ( { null, null, null } );
      
      
      
      int count = 0;
	
      const GridView gridView = grid.leafGridView();

      fvector normalOld;
      Line2D<fvector> g;
      
      for ( LeafIterator leaf = gridView.template begin<0>(); leaf != gridView.template end<0>(); ++leaf )
      {
	count = 0;
	    
	const Entity &entity = *leaf;  
	int entityIndex = gridView.indexSet().index( entity );
	
	if ( cellIsMixed[ entityIndex ] && concentration[ entityIndex ] <= 1 - eps ) 
	{
	  const Geometry geoEntity = entity.geometry();

	
	  guessNormal( grid, concentration, entity, entityIndex, domain, g );
	  
	  
	  double sumCount;
	  fvector centroidLineNormal, sumNormals;
	    
	  do 
	  {
	    sumCount = 0;
	    sumNormals = 0;

	    
	    auto entityIntersections = 
	      computeInterfaceLinePosition( gridView, entity, geoEntity, concentration[ entityIndex ], g );
	    
	    fvector centroid1 = entityIntersections[0] + entityIntersections[1];
	    centroid1 *= 0.5;
	    
	    for( int neighborIndex : domain.cellsInDomain[ entityIndex ] )
	    {
		      
		if ( cellIsMixed[ neighborIndex ] && concentration[ neighborIndex ] <= 1 - eps ) 
		{
		  
		  const Entity &neighbor = grid.entity( domain.seeds[ neighborIndex ] );
			  
		  Line2D<fvector> h;
		  
		  // nimm diese Zelle nur, wenn die geratene Normal schon in eine aehnliche Richtung zeigt
		  guessNormal( grid, concentration, neighbor, neighborIndex, domain, h );
		  if ( g.n * h.n < 0 )
		    continue;
		  
		  
		  
		  const Geometry geoNeighbor = neighbor.geometry();

		  h.n = g.n;     
			  
		  auto neighborIntersections = 
		    computeInterfaceLinePosition( gridView, neighbor, geoNeighbor, concentration[ neighborIndex ], h );  	      
			
		  
		  assert( entityIntersections.size() == 2 );
		  assert( neighborIntersections.size() == 2 );
		  
		  fvector centroid2 = neighborIntersections[0] + neighborIntersections[1];
		  centroid2 *= 0.5;
				  
		  centroidLineNormal = centroid2 - centroid1;    
		  rotate90degreesCounterClockwise( centroidLineNormal );
		  
		  
		  assert( centroidLineNormal.two_norm() > 0 );
		  
		  
		  
		  // errechnete Normale muss in eine aehnliche Richtung wie geschaetze Normale auf jeder Zelle zeigen
		  if ( centroidLineNormal * g.n < 0 )
		    centroidLineNormal *= -1.0;      


		  sumNormals.axpy( 1.0 / (geoEntity.center() - centroid2).two_norm(), centroidLineNormal );
		  sumCount++;
		}
	      
	    }
	    

	    normalOld = g.n;
	    
	    if( sumCount > 0 )
	    {
	      g.n = sumNormals;
	      g.n *= 1.0 / sumCount;
	    }
	    
	    assert( g.n.two_norm() > 1e-14 );
	      
	    count++;
	    
	  } while ( (g.n - normalOld).two_norm() > 1e-8 && count < 30 );   // limit number of loops
	
	  
	  auto entityIntersections = 
	      computeInterfaceLinePosition( gridView, entity, geoEntity, concentration[ entityIndex ], g );
	      
	  reconstruction[entityIndex] = std::array<fvector,3> ( { entityIntersections[0], entityIntersections[1], g.n } );
	  
	}
      }
      
      return;
    }



    template < class G, class C, class E, class V, class D >
    void guessNormal( const G &grid, const C &c, const E &entity, const int entityIndex, const D& domain, Line2D<V>& g )
    {

      const int n = domain.cellsInDomain[ entityIndex ].size();
      
      Dune::DynamicMatrix<double> A ( n, 2 );
      Dune::DynamicVector<double> b ( n );
      
      V xi = entity.geometry().center();
      
      for ( int i = 0; i < n; ++i )
      {
	int neighborIndex = domain.cellsInDomain[ entityIndex ][i];
	const E &neighbor = grid.entity( domain.seeds[ neighborIndex ] );
	
	
	V xk = neighbor.geometry().center();
	double wk = 1.0 / ( sqr ( (xk - xi).one_norm() ) );

	A[i][0] = wk * ( xk[0] - xi[0] );
	A[i][1] = wk * ( xk[1] - xi[1] );
	
	b[i] = wk * ( c[ neighborIndex ] - c[ entityIndex ] );
      }
      
      Dune::DynamicMatrix<double> AtA ( 2, 2 );
      
      for( int i = 0; i < 2; ++i )
	for ( int j = 0; j < 2; ++j )
	{
	  AtA[i][j] = 0;
	  
	  for ( int k = 0; k < n; ++k )
	    AtA[i][j] += A[k][i] * A[k][j];
	}
      
      Dune::DynamicVector<double> x ( 2 );
      Dune::DynamicVector<double> Atb ( 2 );
      
      A.mtv( b, Atb );
      
      AtA.solve( x, Atb ); 
      
      assert ( x.two_norm() > 0 );
      x *= 1.0 / x.two_norm();
      
      g.n = x;
      
    }



    template < class GV, class E, class Geo, class V >
    std::vector<V> computeInterfaceLinePosition( const GV& gridView, const E& entity, const Geo& geo, double concentration, Line2D<V>& g )
    {  
      double pMin = 0, pMax = 0, volume;
      //use bigger range than [0,1] initially
      double volMin = -1; 
      double volMax = 2;
	
      // Initial guess for p
      for ( int i = 0; i < geo.corners(); ++i )
      {
	g.p = geo.corner(i) * g.n;
	g.p = -g.p;
	  
	volume = getVolumeFraction( gridView, entity, geo, g );
	
	
	if ( volume <= volMax && volume >= concentration )
	{
	    pMax = g.p;
	    volMax = volume;
	}
	
	if ( volume >= volMin && volume <= concentration )
	{
	    pMin = g.p;
	    volMin = volume;
	}
      }
	
      
      g.p = brentsMethod( pMin, pMax, 
		  [&gridView, &entity, &geo, &concentration, &g] ( double p ) -> double
		  { Line2D<V> h (g.n, p); return getVolumeFraction ( gridView, entity, geo, h) - concentration; } );
      
      return lineIntersectionPoints( gridView, entity, geo, g );
    }



    //Quelle: Wikipedia (Brent-Verfahren)
    template < class F >
    double brentsMethod ( double a, double b, F f )
    {
      const double TOL = 1e-14;
      
      double fa,fb,fc,c,d,e,p,q,m,s,tol,r;
      
      
      fa = f(a);
      fb = f(b);
      
      if ( std::abs( fa ) < TOL ) return a;
      if ( std::abs( fb ) < TOL ) return b;
      
      assert( fa * fb < 0 );
      
      c = a; fc = fa;
      d = b - a; e = d;
    
      int iter = 0;
      int maxiter = 1000;
    
      while ( iter < maxiter )
      {
	
	iter++;
    
	if ( fb * fc > 0 )
	{
	  c = a; 
	  fc = fa; 
	  d = b - a;
	  e = d;
	}
    
	if ( std::abs( fc ) < std::abs( fb ) )
	{
	  a = b; 
	  b = c; 
	  c = a;
	  fa = fb; 
	  fb = fc; 
	  fc = fa;
	}
	
	tol = 2.0 * 1e-8 * std::abs(b) + TOL;
	m = ( c - b ) / 2.0; 
    
	if ( std::abs( m ) > tol && std::abs( fb ) > 0 )
	{
	  
	  if ( std::abs( e ) < tol || std::abs( fa ) <= std::abs( fb ) )
	  {
	    d = m;
	    e = m;
	  }
	  else
	  {
	    s = fb / fa;
	    
	    if ( a == c )
	    {
	      p = 2.0 * m * s; 
	      q = 1 - s;
	    }
	    else
	    {
	      q = fa / fc; 
	      r = fb / fc;
	      p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
	      q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
	    }
	    
	    if ( p > 0 )
	      q = -q;
	    else
	      p = -p;
	    
	    s = e;
	    e = d;
	    
	    if ( 2.0 * p < 3.0 * m * q - std::abs( tol * q ) && p < std::abs( s * q / 2.0 ) )
	    { 
	      d = p / q;
	    }
	    else
	    {
	      d = m; 
	      e = m;
	    }
	  }
	  
	  a = b;
	  fa = fb;
	
	  if ( std::abs( d ) > tol )
	  {
	    b = b + d;
	  }
	  else
	  {
	    if ( m > 0 )
	    {
	      b = b + tol;
	    }
	    else
	    {
	      b = b - tol;
	    }
	  }
	}
	else
	{
	  break;
	}
    
	fb = f(b);
      }

      return b;

    }


  }
}

#endif

