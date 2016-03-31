#ifndef DUNE_VOF_FEMDFWRAPPER_HH
#define DUNE_VOF_FEMDFWRAPPER_HH

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>


namespace Dune
{

  namespace VoF
  {


    // FemDiscreteFunctionWrapper
    // --------------------------

    /**
     * \ingroup Other
     * \brief wrapper for a dune-fem discrete function
     *
     * \tparam  DF  discrete function type
     */
    template< class DF >
    struct FemDiscreteFunctionWrapper
    {
      typedef DF DiscreteFunction;
      typedef typename DiscreteFunction::GridPartType GridView;
      typedef typename GridView::template Codim< 0 >::EntityType Entity;
      typedef typename DiscreteFunction::RangeFieldType ctype;

    public:
      FemDiscreteFunctionWrapper ( DiscreteFunction &discreteFunction )
      : discreteFunction_( discreteFunction )
      {}

      ctype& operator[] ( const Entity& entity ) { return discreteFunction().localFunction( entity )[0]; }
      const ctype& operator[] ( const Entity& entity ) const { return discreteFunction().localFunction( entity )[0]; }

      const GridView &gridView () const { return discreteFunction().gridPart(); }

      void axpy ( const ctype a, DiscreteFunction &x )
      {
        discreteFunction().axpy( a, x );
      }

      void clear () { discreteFunction().clear(); }

      DiscreteFunction& discreteFunction () { return discreteFunction_; }
      const DiscreteFunction& discreteFunction () const { return discreteFunction_; }

    private:
      DiscreteFunction &discreteFunction_;
    };


    template< class DiscreteFunction >
    inline static FemDiscreteFunctionWrapper< DiscreteFunction > discreteFunctionWrapper ( DiscreteFunction &discreteFunction )
    {
      return FemDiscreteFunctionWrapper< DiscreteFunction > ( discreteFunction );
    }


  } // namespace Dune

} // namespace VoF

#endif // #ifndef DUNE_VOF_FEMDFWRAPPER_HH
