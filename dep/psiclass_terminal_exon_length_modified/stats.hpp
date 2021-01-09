#ifndef _LSONG_STATS
#define _LSONG_STATS

#include <stdio.h>
#include <math.h>

#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

/** The digamma function in long double precision.
* @param x the real value of the argument
* @return the value of the digamma (psi) function at that point
* @author Richard J. Mathar
* @since 2005-11-24
*/
long double digammal(long double x) ;

//https://people.sc.fsu.edu/~jburkardt/cpp_src/asa121/asa121.hpp
double trigamma ( double x, int *ifault );


double LogGammaDensity( double x, double k, double theta ) ;
double MixtureGammaAssignment( double x, double pi, double* k, double *theta ) ;


// http://people.sc.fsu.edu/~jburkardt/c_src/asa091/asa091.hpp
double alnorm ( double x, bool upper );
double gammad ( double x, double p, int *ifault );
double r8_min ( double x, double y );

double chicdf( double x, double df ) ;

#endif
