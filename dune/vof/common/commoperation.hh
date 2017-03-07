#ifndef DUNE_VOF_COMMON_COMMOPERATION_HH
#define DUNE_VOF_COMMON_COMMOPERATION_HH

namespace Dune
{
  namespace VoF
  {

    namespace CommOperation
    {
      struct Add
      {
        template< class T >
        T operator()( T a, T b ) const { return a + b; }
      };

      struct Copy
      {
        template< class T >
        T operator()( T a, T b ) const { return a; }
      };
    }

  }
}

#endif // #ifndef DUNE_VOF_COMMON_COMMOPERATION_HH
