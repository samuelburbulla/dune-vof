#ifndef DFWRAPPER_HH
#define DFWRAPPER_HH

//- dune-grid includes
#include <dune/vof/mcmgmapper.hh>


// DFWRAPPER
// -------------

namespace Dune {
  namespace VoF {

    template< class DF >
    struct ColorFunction
    {
      typedef DF DiscreteFunction;
      typedef typename DiscreteFunction::GridPartType GridView;
      typedef typename GridView::template Codim< 0 >::EntityType Entity;
      typedef typename DiscreteFunction::RangeFieldType ctype;

    public:
      ColorFunction ( DiscreteFunction &discreteFunction )
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
    inline static ColorFunction< DiscreteFunction > colorFunction ( DiscreteFunction &discreteFunction ) 
    {
      return ColorFunction< DiscreteFunction > ( discreteFunction );
    }


  }
}

#endif // #ifndef DFWRAPPER_HH


