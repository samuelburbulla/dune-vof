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
      using VertexStencilSet = StS;
      using InitialReconstruction = IR;

      using GridView = typename ColorFunction::GridView;

    private:
      using Reconstruction = typename ReconstructionSet::DataType;
      using VertexStencil = typename VertexStencilSet::Stencil;

      using Entity = typename ColorFunction::Entity;
      using Coordinate = typename Entity::Geometry::GlobalCoordinate;

      using Stencil = HeightFunctionStencil< GridView >;

      static constexpr int dim = Coordinate::dimension;

      using Heights = Dune::FieldVector< double, Stencil::noc >;
      using Orientation = std::tuple< int, int >;

    public:
      explicit HeightFunctionReconstruction ( const VertexStencilSet &vertexStencilSet )
       : vertexStencilSet_( vertexStencilSet ), initializer_( vertexStencilSet )
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
        const Orientation orientation = getOrientation( normal );

        const auto entityInfo = GridView::Grid::getRealImplementation( entity ).entityInfo();
        const auto stencil = Stencil( color.gridView(), entityInfo, orientation );

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
        //double uMid = heights[ ( heights.size() - 1 ) / 2 ];
        //int effTdown = stencil.effectiveTdown();
        //if ( uMid < effTdown || uMid > effTdown + 1 )
          //return;

        Coordinate newNormal = computeNormal( heights, orientation );

        reconstruction = locateHalfSpace( makePolytope( entity.geometry() ), newNormal, color[ entity ] );
      }

    private:
    #if GRIDDIM == 2
      Coordinate computeNormal ( const Heights &heights, const Orientation &orientation ) const
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

        int i = std::get< 0 >( orientation );
        int j = std::get< 1 >( orientation );

        if ( ( i == 0 && j == 1 ) || ( i == 1 && j == -1 ) )
          newNormal *= -1.0;
        if ( ( i == 0 && j == 1 ) || ( i == 0 && j == -1 ) )
          newNormal = generalizedCrossProduct( newNormal );

        return newNormal;
      }

    #elif GRIDDIM == 3
      Coordinate computeNormal ( const Heights &heights, const Orientation &orientation ) const
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


        int i = std::get< 0 >( orientation );
        int j = std::get< 1 >( orientation );

        Coordinate newNormal;
        newNormal[ i ] = -j;
        newNormal[ (i+1)%3 ] = Hx * -j;
        newNormal[ (i+2)%3 ] = Hy;

        normalize( newNormal );
        return newNormal;
      }
    #endif

      static inline Orientation getOrientation( const Coordinate &normal )
      {
        std::size_t dir = 0;
        double max = std::numeric_limits< double >::min();
        for ( std::size_t i = 0; i < dim; ++i )
        if ( std::abs( normal[ i ] ) > max )
        {
          dir = i;
          max = std::abs( normal[ i ] );
        }
        int sign = ( normal[ dir ] > 0 ) ? -1.0 : 1.0;

        return std::make_tuple( dir, sign );
      }

      const VertexStencil &vertexStencil ( const Entity &entity ) const { return vertexStencilSet_[ entity ]; }

      const VertexStencilSet &vertexStencilSet_;
      InitialReconstruction initializer_;
    };

  }       // end of namespace VoF
} // end of namespace Dune

#endif // #ifndef DUNE_VOF_RECONSTRUCTION_MODIFIEDYOUNGS_HH
