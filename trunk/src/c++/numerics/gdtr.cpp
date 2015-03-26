/*							gdtr.c
 *
 *	Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, gdtr();
 *
 * y = gdtr( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from zero to x of the gamma probability
 * density function:
 *
 *
 *                x
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               0
 *
 *  The incomplete gamma integral is used, according to the
 * relation
 *
 * y = igam( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See igam().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtr domain         x < 0            0.0
 *
 */
/*							gdtrc.c
 *
 *	Complemented gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, gdtrc();
 *
 * y = gdtrc( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from x to infinity of the gamma
 * probability density function:
 *
 *
 *               inf.
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               x
 *
 *  The incomplete gamma integral is used, according to the
 * relation
 *
 * y = igamc( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See igamc().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtrc domain         x < 0            0.0
 *
 */

/*							gdtr()  */


/*
Cephes Math Library Release 2.3:  March,1995
Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

#include "numerics/cephes.h"

namespace SAGE
{

double gdtr( double a, double b, double x )
{

  if( x < 0.0 )
    return std::numeric_limits<double>::quiet_NaN();

  return(  igam( b, a * x )  );
}



double gdtrc( double a, double b, double x )
{
  if( x < 0.0 )
    return std::numeric_limits<double>::quiet_NaN();

  return(  igamc( b, a * x )  );
}

}
