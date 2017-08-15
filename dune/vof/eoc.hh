#ifndef DUNE_VOF_EOC_HH
#define DUNE_VOF_EOC_HH

#include <iostream>
#include <fstream>

namespace Dune
{

  struct EocOutput
  {
  public:
    EocOutput ( const std::string &filename, const double dx0, const int level0 = 0, const double refinementFactor = 2.0 )
     : refinementFactor_( refinementFactor ), lnH( std::log( refinementFactor ) ), dx( dx0 ), level0( level0 ), level( level0 )
    {
      errorsfile.open ( filename );
      errorsfile << "#    dx   \t error " << std::endl;
    };

    ~EocOutput ()
    {
      errorsfile.close();
    }

    void add ( const double error )
    {
      const double eoc = std::log( oldError / error ) / lnH;

      std::cout << "Error =" << std::setw( 16 ) << error << std::endl;

      errorsfile << std::setw( 10 ) << dx << "\t" << error << std::endl;

      if ( level > level0 )
      {
        std::cout << "   EOC( " << std::setw( 2 ) << level << " ) = " << std::setw( 11 ) << eoc << std::endl;
      }

      ++level;
      dx /= refinementFactor_;
      oldError = error;
    }

  private:
    std::ofstream errorsfile;
    const double refinementFactor_;
    const double lnH;
    double dx;
    const int level0;
    int level;
    double oldError = -1.0;
  };

}  // namespace Dune

#endif // #ifndef DUNE_VOF_EOC_HH
