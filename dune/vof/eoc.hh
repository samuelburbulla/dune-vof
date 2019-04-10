#ifndef DUNE_VOF_EOC_HH
#define DUNE_VOF_EOC_HH

#include <cmath>

#include <iomanip>
#include <iostream>
#include <fstream>

namespace Dune
{

  struct EocOutput
  {
  public:
    EocOutput ( const std::string &filename, const double dx0, const int level0, const double refinementFactor, const int dim )
     : refinementFactor_( refinementFactor ), lnH( std::log( refinementFactor ) ), dx( dx0 ), level0( level0 ), level( level0 ), dim_( dim )
    {
      errorsfile.open ( filename );
      errorsfile << "#   dx   \t error " << std::endl;

      errortable.open ( filename + "_latextable" );
      errortable << "Gitter & $E_{L^1}$ & EOC   \\\\" << std::endl << "\\hline \\\\" << std::endl;
    };

    ~EocOutput ()
    {
      errorsfile.close();
      errortable.close();
    }

    void add ( const double error )
    {
      const double eoc = std::log( oldError / error ) / lnH;

      std::cout << "Error =" << std::setw( 16 ) << error << std::endl;

      errorsfile << std::setw( 10 ) << dx << "\t" << error << std::endl;

      int cells = 8 * std::pow( 2, level );

      errortable << "$" << cells << "^" << dim_ << "$     \t & " << parseExponent( error ) << " & \t ";

      if ( level == level0 )
        errortable << "\\ \\ \\ - \\\\" << std::endl;
      else
      {
        errortable << std::fixed << std::setprecision( 2 ) << eoc << " \\" << std::endl;

        std::cout << "   EOC( " << std::setw( 2 ) << level << " ) = " << std::setw( 11 ) << eoc << std::endl;
      }

      ++level;
      dx /= refinementFactor_;
      oldError = error;
    }

  private:

    std::string parseExponent ( const double number )
    {
      std::stringstream parsed;
      parsed << std::scientific << std::uppercase << std::setprecision( 2 ) << number << std::nouppercase;
      std::string parsedstr = parsed.str();
      auto pos = parsedstr.find('E') + 2;
      while ( parsedstr[ pos ] == '0' )
        parsedstr.erase( pos, 1 );

      return parsedstr;
    }

    std::ofstream errorsfile, errortable;
    const double refinementFactor_;
    const double lnH;
    double dx;
    const int level0;
    int level;
    const int dim_;
    double oldError = -1.0;
  };

}  // namespace Dune

#endif // #ifndef DUNE_VOF_EOC_HH
