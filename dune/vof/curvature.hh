#ifndef DUNE_VOF_CURVATURE_HH
#define DUNE_VOF_CURVATURE_HH

#include <algorithm>
#include <cmath>
#include <utility>

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

    // Calculate the curvature of the inferface in each mixed cell.

    /**
     * \ingroup Method
     * \brief set of curvatures
     *
     * \tparam  GV  grid view
     * \tparam  ST  stencils
     * \tparam  RS  reconstruction set
     * \tparam  FL  flags
     */
    template< class GV, class ST, class DF, class RS, class FL >
    struct Curvature
    {
      using GridView = GV;
      using Stencils = ST;
      using DiscreteFunction = DF;
      using ReconstructionSet = RS;
      using Flags = FL;
      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      static constexpr int dim = Coordinate::dimension;
      static constexpr int derivatives = dim; //dim + dim * ( dim + 1 ) / 2;
      using ctype = typename Coordinate::value_type;
      using Matrix = FieldMatrix< ctype, derivatives, derivatives >;
      using Vector = FieldVector< ctype, derivatives >;

    private:
      using IndexSet = decltype( std::declval< GridView >().indexSet() );
      using Index = decltype( std::declval< IndexSet >().index( std::declval< Entity >() ) );

    public:
      explicit Curvature ( const GridView &gridView, const Stencils &stencils )
       : gridView_( gridView ), stencils_( stencils ), curvature_( indexSet().size( 0 ) )
      {}

      void operator() ( const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags )
      {
        for ( const auto& entity : elements( gridView() ) )
        {
          curvature_[ index( entity ) ] = 0.0;

          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applyLocal( entity, uh, reconstructions, flags );
        }

        auto tmpCurvature1 ( curvature_ );
        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applySmoothing1st( entity, uh, tmpCurvature1 );
        }

        auto tmpCurvature2 ( curvature_ );
        for ( const auto& entity : elements( gridView() ) )
        {
          if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
            continue;

          applySmoothing2nd( entity, uh, reconstructions, tmpCurvature2 );
        }

      }

      double operator[] ( const Entity& entity ) const
      {
        return curvature_[ index( entity ) ];
      }

      double operator[] ( const int index ) const
      {
        return curvature_[ index ];
      }

      std::size_t size () const { return curvature_.size(); }

      const GridView &gridView () const { return gridView_; }

    private:
      double smoothingWeight( const double color )
      {
        return std::pow( 1.0 - 2.0 * std::abs( 0.5 - color ), 8.0 );
      }

      void applySmoothing1st ( const Entity &entity, const DiscreteFunction &uh, const std::vector< double > tmpCurvature )
      {
        double kappa = tmpCurvature[ index( entity ) ] * smoothingWeight( uh[ entity ] );
        double sumWeights = smoothingWeight( uh[ entity ] );

        for ( const auto &neighbor : stencils_[ entity ] )
        {
          kappa += tmpCurvature[ index( neighbor ) ] * smoothingWeight( uh[ neighbor ] );
          sumWeights += smoothingWeight( uh[ neighbor ] );
        }

        kappa /= sumWeights;
        curvature_[ index( entity ) ] = kappa;
      }

      double additionalSmoothingWeight( const Coordinate &normal, const Coordinate &delta )
      {
        return std::pow( std::abs( normal * delta ), 8.0 );
      }

      void applySmoothing2nd ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const std::vector< double > tmpCurvature )
      {
        Coordinate normal = reconstructions[ entity ].innerNormal();

        double kappa = tmpCurvature[ index( entity ) ] * smoothingWeight( uh[ entity ] );
        double sumWeights = smoothingWeight( uh[ entity ] );

        for ( const auto &neighbor : stencils_[ entity ] )
        {
          Coordinate delta = neighbor.geometry().center() - entity.geometry().center();
          delta /= delta.two_norm();

          kappa += tmpCurvature[ index( neighbor ) ] * smoothingWeight( uh[ neighbor ] ) * additionalSmoothingWeight( normal, delta );
          sumWeights += smoothingWeight( uh[ neighbor ] ) * additionalSmoothingWeight( normal, delta );
        }

        kappa /= sumWeights;
        curvature_[ index( entity ) ] = kappa;
      }

      void applyLocal ( const Entity &entity, const DiscreteFunction &uh, const ReconstructionSet &reconstructions, const Flags &flags )
      {
        /*
        double tol = std::numeric_limits< double >::epsilon();
        // Height function
        Coordinate normal = reconstructions[ entity ].innerNormal();
        Coordinate center = entity.geometry().center();

        Coordinate direction ( 0.0 );
        Coordinate directionOrth ( 0.0 );
        if ( std::abs( normal[ 0 ] ) > std::abs( normal[ 1 ] ) )
        {
         direction[ 0 ] = 1.0;
         directionOrth[ 1 ] = 1.0;
        }
        else
        {
         direction[ 1 ] = 1.0;
         directionOrth[ 0 ] = 1.0;
        }

        double left = 0.0, right = 0.0;
        double middle = uh[ entity ];
        double dx = std::numeric_limits< double >::max();
        for( const auto& neighbor : stencils_[ entity ] )
        {
          Coordinate d = center - neighbor.geometry().center();

          if ( std::abs( d * directionOrth ) < tol )
            middle += uh[ neighbor ];
          else if ( d * directionOrth > 0 )
            left += uh[ neighbor ];
          else if ( d * directionOrth < 0 )
            right += uh[ neighbor ];

          dx = std::min( dx, d.two_norm() );
        }

        right *= dx;
        left *= dx;
        middle *= dx;

        if ( dx < middle && middle < 2 * dx )
        {
          // Workaround for boundary! (linear interpolation)
          if ( right < tol )
            right = middle + ( middle - left );
          if ( left < tol )
            left = middle + ( middle - right );

          double Hx = ( right - left ) / ( 2 * dx );
          double Hxx = ( right - 2 * middle + left ) / ( dx * dx );

          curvature_[ index( entity ) ] = - Hxx / std::pow( 1.0 + Hx * Hx, 3.0 / 2.0 );
        }
        */

        // Least squares for gradients
        //auto interfaceEn = interface( entity, reconstructions );
        //Coordinate center = interfaceEn.centroid();
        Coordinate center = entity.geometry().center();

        Matrix nablaN( 0.0 );

        for ( std::size_t k = 0; k < dim; ++k )
        {
          Matrix AtA( 0.0 );
          Coordinate Atb( 0.0 );

          for( const auto &neighbor : stencils_[ entity ] )
          {
            //if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
              //continue;
            //auto interfaceNb = interface( neighbor, reconstructions );
            //Coordinate centerNb = interfaceNb.centroid();

            Coordinate centerNb = neighbor.geometry().center();

            Coordinate d = centerNb - center;
            const ctype weight = 1.0 / d.two_norm2();
            d *= weight;
            AtA += outerProduct( d, d );
            Atb.axpy( weight * ( reconstructions[ neighbor ].innerNormal()[ k ] - reconstructions[ entity ].innerNormal()[ k ] ), d );
          }

          Coordinate dNk ( 0.0 );
          AtA.solve( dNk, Atb );
          curvature_[ index( entity ) ] -= dNk[ k ];
          //nablaN[ k ] = dNk;
        }

        //Coordinate n = reconstructions[ entity ].innerNormal();
        //Coordinate nablaNn;
        //nablaN.mv( n, nablaNn );
        //curvature_[ index( entity ) ] -= n * nablaNn;




        /*
        // Second order Taylor series expansion with least squares
        if ( stencils_[ entity ].size() < 8 )
          return;

        const Coordinate center = entity.geometry().center();
        const double colorEn = uh[ entity ];

        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );
        for ( const auto &neighbor : stencils_[ entity ] )
        {
          Coordinate d = neighbor.geometry().center() - center;

          Vector v; // v = ( dx,  dy, dxx, dyy, dxy )
          for ( std::size_t i = 0; i < dim; ++i )
            v[ i ] = d[ i ];

          for ( std::size_t i = 0; i < dim; ++i )
            v[ i + dim ] = 0.5 * d[ i ] * d[ i ];

          int n = 0;
          for ( std::size_t i = 0; i < dim; ++i )
            for ( std::size_t j = i+1; j < dim; ++j )
              if ( i == j )
                continue;
              else
              {
                v[ n + 2 * dim ] = d[ i ] * d[ j ];
                n++;
              }

          const ctype weight = 1.0 / d.two_norm2();
          v *= weight;
          AtA += outerProduct( v, v );
          Atb.axpy( weight * ( uh[ neighbor ] - colorEn ), v );
        }

        Vector solution( 0.0 );
        AtA.solve( solution, Atb );

        double dxu = solution[ 0 ];
        double dyu = solution[ 1 ];
        double dxxu = solution[ 2 ];
        double dyyu = solution[ 3 ];
        double dxyu = solution[ 4 ];

        if ( dxu * dxu + dyu * dyu < std::numeric_limits< double >::epsilon() )
          curvature_[ index( entity ) ] = 0.0;
        else
          curvature_[ index( entity ) ] = - ( dxxu * dyu * dyu + dyyu * dxu * dxu - 2.0 * dxu * dyu * dxyu ) / std::pow( dxu * dxu + dyu * dyu, 3.0 / 2.0 );
        */

        /*
        // New finite differences
        double AtA = 0.0;
        double Atb = 0.0;
        double BtB = 0.0;
        double Btb = 0.0;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();
        Coordinate normalOrth = generalizedCrossProduct( normalEn );

        double xe = normalOrth * centroidEn;
        double ue = normalEn * centroidEn;

        for( const auto& neighbor1 : stencils_[ entity ] )
        {
          if ( ( !flags.isMixed( neighbor1 ) && !flags.isFullAndMixed( neighbor1 ) ) )
            continue;

          auto interfaceNb1 = interface( neighbor1, reconstructions );
          Coordinate centroidNb1 = interfaceNb1.centroid();

          double xn1 = normalOrth * centroidNb1;
          double un1 = normalEn * centroidNb1;

          double dx = ( xn1 - xe );
          BtB += dx * dx;
          Btb += dx * ( un1 - ue );

          for( const auto& neighbor2 : stencils_[ entity ] )
          {
            if ( ( !flags.isMixed( neighbor2 ) && !flags.isFullAndMixed( neighbor2 ) ) )
              continue;

            if ( neighbor1 == neighbor2 )
              continue;

            auto interfaceNb2 = interface( neighbor2, reconstructions );
            Coordinate centroidNb2 = interfaceNb2.centroid();

            double xn2 = normalOrth * centroidNb2;
            double un2 = normalEn * centroidNb2;

            double tmp = 0.5 * ( xn2 - xn1 ) * ( xn2 - xe ) * ( xn1 - xe );
            AtA += tmp * tmp;
            Atb += tmp * ( ( un2 - ue ) * ( xn1 - xe ) - ( un1 - ue ) * ( xn2 - xe ) );
          }
        }
        double dxu = Btb / BtB;
        double dx2u = Atb / AtA;

        curvature_[ index( entity ) ] = dx2u / std::pow( 1.0 + dxu * dxu, 3.0 / 2.0 );
        */
        /*
        // Finite differences
        double divN ( 0.0 );
        double dN ( 0.0 );
        int n = 0;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffC = centroidEn - centroidNb;

          double dx = generalizedCrossProduct( normalEn ) * diffC;
          double dy = normalEn * diffC;
          divN += dy / ( dx * dx );
          dN += dy / dx;
          n++;
        }
        divN /= n / 2.0;
        dN /= n / 2.0;

        curvature_[ index( entity ) ] = - divN / std::pow( 1 + dN * dN, 3.0 / 2.0 );
        */
        /*
        // Interpolate with circle through points
        Matrix AtA( 0.0 );
        Coordinate Atb( 0.0 );
        std::vector< Coordinate > points;

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          points.push_back( centroidNb );
        }

        for ( const auto &point : points )
        {
          Coordinate p = point;
          p -= centroidEn;
          const ctype weight = p.two_norm2();
          p *= 2.0 * weight;
          AtA += outerProduct( p, p );
          Atb.axpy( weight * ( point.two_norm2() - centroidEn.two_norm2() ), p );
        }
        Coordinate centerOfCircle;
        AtA.solve( centerOfCircle, Atb );

        double radius = ( centerOfCircle - centroidEn ).two_norm();
        double sign = - 1 + 2 * ( 0 < ( centerOfCircle - centroidEn ) * reconstructions[ entity ].innerNormal() );
        curvature_[ index( entity ) ] = sign / radius;
        */
        /*
        // Interpolate with parabole through points
        Matrix AtA( 0.0 );
        Vector Atb( 0.0 );

        auto interfaceEn = interface( entity, reconstructions );
        Coordinate centroidEn = interfaceEn.centroid();
        Coordinate normalEn = reconstructions[ entity ].innerNormal();

        Vector p ( { 0, 0, 1 } );
        const ctype weight = 1.0;
        p *= weight;
        AtA += outerProduct( p, p );
        Atb.axpy( 0.0, p );

        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( !flags.isMixed( neighbor ) && !flags.isFullAndMixed( neighbor ) )
            continue;

          auto interfaceNb = interface( neighbor, reconstructions );
          Coordinate centroidNb = interfaceNb.centroid();

          Coordinate diffC = centroidNb - centroidEn;

          double dx = generalizedCrossProduct( normalEn ) * diffC;
          double dy = normalEn * diffC;


          Vector p ( { dx * dx, dx, 1.0 } );
          const ctype weight = 1.0;
          p *= weight;
          AtA += outerProduct( p, p );
          Atb.axpy( weight * dy, p );

          // with normal
          Coordinate n = reconstructions[ neighbor ].innerNormal();
          Vector p2 ( { 2.0 * dx, 1.0, 0.0 } );
          const ctype weight2 = 1.0;
          p2 *= weight2;
          AtA += outerProduct( p2, p2 );
          Atb.axpy( weight2 * ( - ( n * generalizedCrossProduct( normalEn ) ) / ( n * normalEn ) ), p2 );

        }
        Vector abc;
        AtA.solve( abc, Atb );

        curvature_[ index( entity ) ] = 2.0 * abc[0] / std::pow( 1.0 + abc[1] * abc[1], 3.0 / 2.0 );
        */
      }

      void averageCurvature( const Entity &entity )
      {
        int n = 0;
        for( const auto& neighbor : stencils_[ entity ] )
        {
          if ( curvature_[ index( neighbor ) ] == 0.0 )
            continue;

          curvature_[ index( entity ) ] += curvature_[ index( neighbor ) ];
          n++;
        }
        curvature_[ index( entity ) ] /= n;
      }

      auto interface( const Entity &entity, const ReconstructionSet &reconstructions ) const
      {
        auto polygon = makePolytope( entity.geometry() );
        auto it = intersect( polygon, reconstructions[ entity ].boundary() );
        auto interface = static_cast< typename decltype( it )::Result > ( it );
        return interface;
      }

      Matrix outerProduct ( const Vector &a, const Vector &b ) const
      {
        Matrix m( 0.0 );
        for ( std::size_t i = 0; i < derivatives; ++i )
          m[ i ].axpy( a[ i ], b );

        return m;
      }

      const IndexSet& indexSet () const { return gridView().indexSet(); }
      Index index ( const Entity &entity ) const { return indexSet().index( entity ); }

      GridView gridView_;
      const Stencils &stencils_;
      std::vector< double > curvature_;
    };

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_CURVATURE_HH
