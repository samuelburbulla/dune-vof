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
    template< class GV, class StS, class IR >
    struct HeightFunctionReconstruction
    {
      using GridView = GV;
      using VertexStencilSet = StS;
      using InitialReconstruction = IR;

    private:
      using VertexStencil = typename VertexStencilSet::Stencil;

      using Entity = typename decltype(std::declval< GridView >().template begin< 0 >())::Entity;
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
       * \tparam  ColorFunction
       * \tparam  ReconstructionSet
       * \tparam  Flags
       * \param   color           color function
       * \param   reconstructions set of interface
       * \param   flags           set of flags
       */
      template< class ColorFunction, class ReconstructionSet, class Flags >
      void operator() ( const ColorFunction &color, ReconstructionSet &reconstructions, const Flags &flags, bool communicate = false ) const
      {
        initializer_( color, reconstructions, flags );

        for ( const auto &entity : elements( color.gridView(), Partitions::interiorBorder ) )
        {
          if ( !flags.isMixed( entity ) )
            continue;

          applyLocal( entity, color, flags, reconstructions[ entity ] );
        }

        if ( communicate )
          reconstructions.communicate();
      }

      /**
       * \brief   (local) operator application
       *
       * \tparam  ColorFunction
       * \tparam  Flags
       * \tparam  Reconstruction
       * \param   entity          current element
       * \param   color           color functions
       * \param   flags           set of flags
       * \param   reconstruction  single reconstruction
       */
      template< class ColorFunction, class Flags, class Reconstruction >
      void applyLocal ( const Entity &entity, const ColorFunction &color, const Flags &flags, Reconstruction &reconstruction ) const
      {
        const Coordinate &normal = reconstruction.innerNormal();
        const Orientation orientation = getOrientation( normal );

        const auto entityInfo = GridView::Grid::getRealImplementation( entity ).entityInfo();
        const auto stencil = Stencil( color.gridView(), entityInfo, orientation );

        Heights heights ( 0.0 );
        double tol = 1e-12;

        for( std::size_t i = 0; i < stencil.columns(); ++i )
        {
          if ( !stencil.valid( i, 0 ) )
            continue;

          double u0 = color[ stencil( i, 0 ) ];
          heights[ i ] += u0;

          // upwards
          double lastU = u0;

          for( int t = 1; t <= stencil.tup(); ++t )
          {
            if ( !stencil.valid( i, t ) )
              break;

            double u = color[ stencil( i, t ) ];

            if ( u > lastU - tol )
              break;

            heights[ i ] += u;
            lastU = u;
          }

          lastU = u0;

          // downwards
          for( int t = -1; t >= stencil.tdown(); --t )
          {
            if ( !stencil.valid( i, t ) )
              break;

            double u = color[ stencil( i, t ) ];

            if ( u < lastU + tol )
              u = 1.0;

            heights[ i ] += u;
            lastU = u;
          }
        }

        // Check constraint
        double uMid = heights[ ( heights.size() - 1 ) / 2 ];
        int effTdown = stencil.effectiveTdown();
        if ( uMid < effTdown || uMid > effTdown + 1 )
          return;

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
