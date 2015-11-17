#ifndef FVSCHEME_HH
#define FVSCHEME_HH

// dune-grid includes
#include <dune/grid/common/rangegenerators.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

template< class DiscreteFunction, class Model >
struct FiniteVolumeScheme
{
  using DiscreteFunctionType = DiscreteFunction;
  using DiscreteFunctionSpaceType =
    typename DiscreteFunctionType::DiscreteFunctionSpaceType;

  using ModelType = Model;

private:
  using GridPartType =
    typename DiscreteFunctionSpaceType::GridPartType;
  using IndexSetType =
    typename GridPartType::IndexSetType;

  using DomainType =
    typename ModelType::DomainType;
  using RangeType =
    typename ModelType::RangeType;
  using ctype =
    typename ModelType::ctype;

public:
  FiniteVolumeScheme ( const DiscreteFunctionSpaceType &space,
                       const ModelType &model )
  : space_( space ),
    model_( model )
  {}

  void operator() ( double time, const DiscreteFunctionType &uh, DiscreteFunctionType &upd ) const
  {
    upd.clear();

    for( auto &&entity : space() )
    {
      auto indexEn = indexSet().index( entity );
      
      const auto uEn = uh.localFunction( entity );
      auto updEn = upd.localFunction( entity );

      double volEn = entity.geometry().volume();

      auto iend = gridPart().iend( entity );
      for( auto iit = gridPart().ibegin( entity ); iit  != iend ; ++iit )
      {
        const auto &intersection = *iit;

        if( intersection.neighbor() )
        {
          const auto &neighbor = intersection.outside();
          auto indexNb  = indexSet().index( neighbor );

          if( indexEn < indexNb )
          {
            const auto uNb = uh.localFunction( neighbor );
            auto updNb = upd.localFunction( neighbor );
            double volNb = neighbor.geometry().volume();

            handleSurface( intersection, time, volEn, volNb, uEn, uNb, updEn, updNb );
          }
        }
        else if( intersection.boundary() )
        {
          handleBoundary( intersection, time, volEn, uEn, updEn );
        }
      }
    }
    
    upd.communicate();
  }

private:
  template< class Intersection, class LocalFunction >
  void handleSurface( const Intersection &intersection, double time,
                      double volEn, double volNb,
                      const LocalFunction &uEn, const LocalFunction &uNb,
                      LocalFunction &updEn, LocalFunction &updNb ) const 
  {
    using ReferenceElementsType
      = Dune::ReferenceElements< typename Intersection::ctype, Intersection::mydimension >;
    const auto &x = ReferenceElementsType::general( intersection.type() ).position( 0, 0 );
    const auto &xEn = intersection.geometryInInside().center();
    const auto &xNb = intersection.geometryInOutside().center();

    RangeType valueEn, valueNb, flux;
    uEn.evaluate( xEn, valueEn );
    uNb.evaluate( xNb, valueNb );

    model().numericalFlux( intersection.integrationOuterNormal( x ),
                           time, intersection.geometry().center(),
                           valueEn, valueNb, flux );

    RangeType fluxEn = flux;
    fluxEn *= -1.0 / volEn;
    updEn.axpy( xEn, fluxEn );

    RangeType fluxNb = flux;
    fluxNb *= 1.0 / volNb;
    updNb.axpy( xNb, fluxNb );
  }

  template< class Intersection, class LocalFunction >
  void handleBoundary( const Intersection &intersection, double time, double volEn,
                       const LocalFunction &uEn, LocalFunction &updEn ) const
  {
    using ReferenceElementsType
      = Dune::ReferenceElements< typename Intersection::ctype, Intersection::mydimension >;
    const auto &x = ReferenceElementsType::general( intersection.type() ).position( 0, 0 );
    const auto &xEn = intersection.geometryInInside().center();

    RangeType valueEn, flux;
    uEn.evaluate( xEn, valueEn );

    model().boundaryFlux( intersection.integrationOuterNormal( x ),
                          time, intersection.geometry().center(),
                          valueEn, flux );

    flux *= -1.0 / volEn;
    updEn.axpy( xEn, flux );
  }

  const DiscreteFunctionSpaceType &space () const { return space_; }
  const GridPartType &gridPart () const { return space().gridPart(); }
  const ModelType &model () const { return model_; }
  const IndexSetType &indexSet() const { return gridPart().indexSet(); }

  const DiscreteFunctionSpaceType &space_;
  const ModelType &model_;
};

#endif // FVSCHEME_HH
