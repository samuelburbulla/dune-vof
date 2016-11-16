#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cmath>

template < class ctype, int dim >
struct RotatingCircle {};


template < class ctype >
struct RotatingCircle < ctype, 2 >
{
  using DomainType = Dune::FieldVector< ctype, 2 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  DomainType center( const double t ) const
  {
    DomainType center = rotationCenter();
    DomainType offset { 0, 0 };
    DomainType normal { 0, 0 };

    center.axpy( std::cos( ( 2 * M_PI / 10 ) * t ), offset );
    center.axpy( std::sin( ( 2 * M_PI / 10 ) * t ), normal );

    return center;
  }

  ctype radius( const double t ) const
  {
    return 0.4;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &r ) const
  {
    DomainType center = rotationCenter();

    r[ 0 ] = x[ 1 ] - center[ 1 ];
    r[ 1 ] = center[ 0 ] - x[ 0 ];

    r *= 2 * M_PI / 10;
  }

private:
  DomainType rotationCenter () const
  {
    return DomainType { 0.5, 0.5 };
  }

};

template < class ctype >
struct RotatingCircle < ctype, 3 >
{
  using DomainType = Dune::FieldVector< ctype, 3 >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = 3 };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    DomainType center{ 0.5, 0.5, 0.5 };
    center.axpy( std::cos( (2 * M_PI / 10)*t ), DomainType{ 0.25, 0.0, 0.0 } );
    center.axpy( std::sin( (2 * M_PI / 10)*t ), DomainType{ 0.0,  0.25, 0.0 } );
    center.axpy( - std::sin( (2 * M_PI / 10)*t ), DomainType{ 0.0, 0.0, 0.25 } );

    double dist = ( x - center ).two_norm();

    u = ( dist < 0.15 ) ? RangeType( 1.0 ) : RangeType( 0.0 );
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    DomainType c = x;
    c -= DomainType{ 0.5, 0.5, 0.5 };

    rot[ 0 ] = - c[ 1 ] + c[ 2 ];
    rot[ 1 ] = c[ 0 ];
    rot[ 2 ] = - c[ 0 ];
    rot *= 2 * M_PI / 10 / std::sqrt(2);
  }

};


template < class ctype, int dim >
struct LinearWall
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    double pos = 0.4;
    double width = 0.3;

    u = ( std::abs( x[0] - t * speed_ - pos ) <= width ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot[ 0 ] = speed_;
    rot[ 1 ] = 0.0;
  }

private:
  double speed_ = 0.02;

};


template < class ctype, int dim >
struct Diagonal
{
  using DomainType = Dune::FieldVector< ctype, dim >;
  using RangeType = Dune::FieldVector< ctype, 1 >;

  enum { dimDomain = dim };
  enum { dimRange = 1 };

  void evaluate ( const DomainType &x, double t, RangeType &u ) const
  {
    double pos = 0.3;

    u = ( x[0] + x[1] < pos + t * speed_ ) ? 1.0 : 0.0;
  }

  void velocityField ( const DomainType &x, const double t, DomainType &rot ) const
  {
    rot[ 0 ] = speed_;
    rot[ 1 ] = speed_;
    rot /= std::sqrt( 2 );
  }

private:
  double speed_ = 0.05;

};

#endif //#ifndef PROBLEM_HH
