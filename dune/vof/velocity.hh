#ifndef VELOCITY_HH
#define VELOCITY_HH

template< class Problem, class GridView >
class Velocity
{
  using DomainType = typename Problem::DomainType;
  using RangeType = typename Problem::RangeType;
  using Intersection = typename GridView::template Codim< 0 >::Entity::Geometry;
  using LocalCoordinate = typename Intersection::LocalCoordinate;

public:
  Velocity ( const Problem &problem, double time ) : problem_ ( &problem ), time_( time ) {}

  DomainType operator() ( const LocalCoordinate &local ) const
  {
    DomainType v ( 0.0 );
    problem().velocityField( intersection().global( local ), time_, v );
    return v;
  }

  void bind ( const Intersection &intersection ) { intersection_ = &intersection; }

  void unbind () { intersection_ = nullptr; }

private:
  const Problem& problem() const { assert( problem_ ); return *problem_; }
  const Intersection& intersection() const { assert( intersection_ ); return *intersection_; }

  const Intersection* intersection_ = nullptr;
  const Problem* problem_;
  double time_;
};


#endif
