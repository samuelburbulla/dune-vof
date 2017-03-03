#ifndef DUNE_VOF_RECONSTRUCTION_HEIGHTFUNCTION_HH
#define DUNE_VOF_RECONSTRUCTION_HEIGHTFUNCTION_HH

#include <cmath>

#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dune/vof/geometry/algorithm.hh>
#include <dune/vof/geometry/utility.hh>

#include <dune/vof/stencil/heightfunctionstencil.hh>


namespace Dune
{
  namespace VoF
  {

    /**
     * \ingroup Reconstruction
     * \brief   height function reconstruction operator
     *
     * \tparam DF   discrete function type
     * \tparam RS   reconstruction set type
     * \tparam StS  stencils type
     * \tparam IR   initial reconstruction type
     */
    template< class DF, class RS, class StS, class IR >
    struct HeightFunctionReconstruction
    {
      using ColorFunction = DF;
      using ReconstructionSet = RS;
      using StencilSet = StS;
      using InitialReconstruction = IR;

      using GridView = typename ColorFunction::GridView;

    private:
      using Reconstruction = typename ReconstructionSet::DataType;
      using Stencil = typename StencilSet::Stencil;

      using Entity = typename ColorFunction::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      using HeightFunctionStencils = Dune::VoF::HeightFunctionStencils< GridView >;

      static constexpr int dim = Coordinate::dimension;

      using Heights = Dune::FieldVector< double, HeightFunctionStencils::Stencil::noc >;

    public:
      explicit HeightFunctionReconstruction ( const StencilSet &stencils )
       : stencils_( stencils ), heightFunctionStencils_( stencils.gridView() ), initializer_( stencils )
      {}

      /**
       * \brief   (global) operator application
       *
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags, bool communicate = false ) const
      {
        initializer_( color, reconstructions, flags );

        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, flags, color, reconstructions[ entity ] );
        }

        if ( communicate )
          reconstructions.communicate();
      }

      /**
       * \brief   (local) operator application
       *
       * \tparam  Flags
       * \param   entity          current element
       * \param   flags           set of flags
       * \param   color           color functions
       * \param   reconstructions set of reconstruction
       */
      template< class Flags >
      void applyLocal ( const Entity &entity, const Flags &flags, const ColorFunction &color, Reconstruction &reconstruction ) const
      {
        const Coordinate &normal = reconstruction.innerNormal();
        auto stencil = heightFunctionStencils_( normal, entity );

        Heights heights ( 0.0 );

        for( std::size_t i = 0; i < stencil.columns(); ++i )
          for( int t = stencil.tdown(); t <= stencil.tup(); ++t )
          {
            if ( !stencil.valid( i, t ) )
              continue;

            double u = color[ stencil( i, t ) ];

            // local monotonic variation
            if ( t < 0 )
            {
              if ( u < color[ stencil( i, t+1 ) ] - 1e-8  )
                u = 1.0;
            }
            else if ( t > 0 )
            {
              if ( u > color[ stencil( i, t-1 ) ] + 1e-8 )
                u = 0.0;
            }

            heights[ i ] += u;
          }

        // Check constraint
        double uMid = heights[ ( heights.size() - 1 ) / 2 ];
        int effTdown = stencil.effectiveTdown();
        if ( uMid < effTdown || uMid > effTdown + 1 )
          return;

        Coordinate newNormal = computeNormal( heights, getOrientation( normal ) );

        reconstruction = locateHalfSpace( makePolytope( entity.geometry() ), newNormal, color[ entity ] );
      }

    private:
      template< class Coordinate >
      auto computeNormal ( const Heights &heights, const Coordinate &orientation ) const
      -> typename std::enable_if< Coordinate::dimension == 2, Coordinate >::type
      {
        double Hx;
        if ( heights[ 0 ] == 0.0 )
          Hx = ( heights[ 2 ] - heights[ 1 ] );
        else if ( heights[ 2 ] == 0.0 )
          Hx = ( heights[ 1 ] - heights[ 0 ] );
        else
          Hx = ( heights[ 2 ] - heights[ 0 ] ) / 2.0;

        Coordinate newNormal ( { Hx, -1.0 } );
        normalize( newNormal );

        if ( orientation == Coordinate( { 0.0, -1.0 } ) || orientation == Coordinate( { 1.0, 0.0 } ) )
          newNormal *= -1.0;
        if ( orientation == Coordinate( { 1.0, 0.0 } ) || orientation == Coordinate( { -1.0, 0.0 } ) )
          newNormal = generalizedCrossProduct( newNormal );

        return newNormal;
      }

      template< class Coordinate >
      auto computeNormal ( const Heights &heights, const Coordinate &orientation ) const
      -> typename std::enable_if< Coordinate::dimension == 3, Coordinate >::type
      {
        double Hx;
        if ( heights[ 5 ] == 0.0 )
          Hx = ( heights[ 4 ] - heights[ 3 ] );
        else if ( heights[ 3 ] == 0.0 )
          Hx = ( heights[ 5 ] - heights[ 4 ] );
        else
          Hx = ( heights[ 5 ] - heights[ 3 ] ) / 2.0;

        double Hy;
        if ( heights[ 7 ] == 0.0 )
          Hy = ( heights[ 5 ] - heights[ 1 ] );
        else if ( heights[ 1 ] == 0.0 )
          Hy = ( heights[ 7 ] - heights[ 5 ] );
        else
          Hy = ( heights[ 7 ] - heights[ 1 ] ) / 2.0;

        Coordinate newNormal;

        if ( orientation == Coordinate( { -1.0, 0.0, 0.0 } ) )
          newNormal = Coordinate( { 1.0, -Hx, Hy } );

        if ( orientation == Coordinate( { 0.0, -1.0, 0.0 } ) )
          newNormal = Coordinate( { -Hx, 1.0, Hy } );

        if ( orientation == Coordinate( { 0.0, 0.0, -1.0 } ) )
          newNormal = Coordinate( { -Hx, Hy, 1.0 } );

        if ( orientation == Coordinate( { 1.0, 0.0, 0.0 } ) )
          newNormal = Coordinate( { -1.0, -Hx, Hy } );

        if ( orientation == Coordinate( { 0.0, 1.0, 0.0 } ) )
          newNormal = Coordinate( { Hx, -1.0, Hy } );

        if ( orientation == Coordinate( { 0.0, 0.0, 1.0 } ) )
          newNormal = Coordinate( { Hx, -Hy, -1.0 } );

        normalize( newNormal );
        return newNormal;
      }

      static inline Coordinate getOrientation( const Coordinate &normal )
      {
        std::size_t dir = 0;
        double max = std::numeric_limits< double >::min();
        for ( std::size_t i = 0; i < dim; ++i )
        if ( std::abs( normal[ i ] ) > max )
        {
          dir = i;
          max = std::abs( normal[ i ] );
        }
        Coordinate orientation;
        orientation[ dir ] = ( normal[ dir ] > 0 ) ? -1.0 : 1.0;

        return orientation;
      }

      const Stencil &stencil ( const Entity &entity ) const { return stencils_[ entity ]; }

      const StencilSet &stencils_;
      const HeightFunctionStencils heightFunctionStencils_;
      InitialReconstruction initializer_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
