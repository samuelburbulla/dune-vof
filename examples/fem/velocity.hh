#ifndef VELOCITY_HH
#define VELOCITY_HH

template< class Problem >
class Velocity
{
  using DomainType = typename Problem::DomainType;
  using RangeType = typename Problem::RangeType;

public:
  Velocity ( const Problem &problem ) : problem_ ( problem ) {}

  DomainType operator() ( const DomainType &x, const double t ) const
  {
    DomainType v ( 0.0 );
    problem_.velocityField( x, t, v );
    return v;
  }

  template < class Intersection >
  void bind ( const Intersection &intersection ) {}

private:
  const Problem &problem_;
};


#endif
