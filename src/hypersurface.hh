#ifndef DUNE_VOF_HYPERSURFACE_HH
#define DUNE_VOF_HYPERSURFACE_HH


namespace Dune
{
  namespace VoF
  {    

    // Representation for a hypersurface n * x + p = 0
    template < class fvector >
    class HyperSurface
    {      
      public:
        typedef typename fvector::value_type ctype;

        HyperSurface () 
         : _n (0), _p (0) 
        {};

        HyperSurface ( const HyperSurface &h ) 
         : _n ( h.normal() ), _p ( h.p() ) 
        {};

        HyperSurface ( const fvector& normal, const ctype skalar ) 
         : _n ( normal ), _p( skalar ) 
        {}

        HyperSurface ( const fvector& normal, const fvector& point ) 
         : _n ( normal )
        {
            _p = normal * point;
            _p *= -1.0;
        }
     

        const fvector& normal() const
        {
          return _n;
        }

        const ctype& p() const
        {
          return _p;
        }


        fvector& normal() 
        {
          return _n;
        }

        ctype& p() 
        {
          return _p;
        }

     private:
        fvector _n; 
        ctype _p;
    };

/*
    template < int dim >
    const fvector intersection (Dune::VoF::HyperSurface< dim > h, HyperSurface g ) {}

    template<>
    const fvector intersection ( Dune::VoF::Hypersurface< 2 >,  ... )



    template< class A, class B >
    struct C;

    template< class B >
    sruct C< double, B >
    {

    };
*/


  } // end of namespace VoF
} // end of namespace Dune

#endif