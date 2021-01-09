#include "stats.hpp"

/** The digamma function in long double precision.
* @param x the real value of the argument
* @return the value of the digamma (psi) function at that point
* @author Richard J. Mathar
* @since 2005-11-24
*/
long double digammal(long double x) 
{
	/* force into the interval 1..3 */
	if( x < 0.0L )
		return digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
	else if( x < 1.0L )
		return digammal(1.0L+x)-1.0L/x ;
	else if ( x == 1.0L)
		return -M_GAMMAl ;
	else if ( x == 2.0L)
		return 1.0L-M_GAMMAl ;
	else if ( x == 3.0L)
		return 1.5L-M_GAMMAl ;
	else if ( x > 3.0L)
		/* duplication formula */
		return 0.5L*(digammal(x/2.0L)+digammal((x+1.0L)/2.0L))+M_LN2l ;
	else
	{
		/* Just for your information, the following lines contain
		* the Maple source code to re-generate the table that is
		* eventually becoming the Kncoe[] array below
		* interface(prettyprint=0) :
		* Digits := 63 :
		* r := 0 :
		* 
		* for l from 1 to 60 do
		* 	d := binomial(-1/2,l) :
		* 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
		* 	evalf(r) ;
		* 	print(%,evalf(1+Psi(1)-r)) ;
		*o d :
		* 
		* for N from 1 to 28 do
		* 	r := 0 :
		* 	n := N-1 :
		*
 		*	for l from iquo(n+3,2) to 70 do
		*		d := 0 :
 		*		for s from 0 to n+1 do
 		*		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
 		*		od :
 		*		if 2*l-n > 1 then
 		*		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
 		*		fi :
 		*	od :
 		*	print(evalf((-1)^n*2*r)) ;
 		*od :
 		*quit :
		*/
		static long double Kncoe[] = { .30459198558715155634315638246624251L,
		.72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
		.27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
		.17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
		.11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
		.83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
		.59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
		.42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
		.304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
		.21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
		.15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
		.11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
		.80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
		.58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
		.41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;

		register long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
		register long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
		register long double resul = Kncoe[0] + Kncoe[1]*Tn ;

		x -= 2.0L ;
		int n ;

		for( n = 2 ; n < int( sizeof(Kncoe)/sizeof(long double) ) ; n++)
		{
			const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
			resul += Kncoe[n]*Tn1 ;
			Tn_1 = Tn ;
			Tn = Tn1 ;
		}
		return resul ;
	}
}



double trigamma ( double x, int *ifault )

//****************************************************************************
//  purpose:
//
//    trigamma calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
//
//  licensing:
//
//    this code is distributed under the gnu lgpl license. 
//
//  modified:
//
//    19 january 2008
//
//  author:
//
//    original fortran77 version by be schneider.
//    c++ version by john burkardt.
//
//  reference:
//
//    be schneider,
//    algorithm as 121:
//    trigamma function,
//    applied statistics,
//    volume 27, number 1, pages 97-99, 1978.
//
//  parameters:
//
//    input, double x, the argument of the trigamma function.
//    0 < x.
//
//    output, int *ifault, error flag.
//    0, no error.
//    1, x <= 0.
//
//    output, double trigamma, the value of the trigamma function at x.
//
{
	double a = 0.0001;
	double b = 5.0;
	double b2 =  0.1666666667;
	double b4 = -0.03333333333;
	double b6 =  0.02380952381;
	double b8 = -0.03333333333;
	double value;
	double y;
	double z;
	//
	//  check the input.
	//
	if ( x <= 0.0 )
	{
		*ifault = 1;
		value = 0.0;
		return value;
	}

	*ifault = 0;
	z = x;
	//
	//  use small value approximation if x <= a.
	//
	if ( x <= a )
	{
		value = 1.0 / x / x;
		return value;
	}
	//
	//  increase argument to ( x + i ) >= b.
	//
	value = 0.0;

	while ( z < b )
	{
		value = value + 1.0 / z / z;
		z = z + 1.0;
	}
	//
	//  apply asymptotic formula if argument is b or greater.
	//
	y = 1.0 / z / z;

	value = value + 0.5 *
		y + ( 1.0 + y * ( b2+ y * ( b4 + y * ( b6+ y * b8 )))) / z;

	return value;
}


double LogGammaDensity( double x, double k, double theta )
{
	return -k * log( theta ) + ( k - 1 ) * log( x ) - x / theta - lgamma( k ) ;
}

double MixtureGammaAssignment( double x, double pi, double* k, double *theta )
{
	if ( pi == 1 )
		return 1 ;
	else if ( pi == 0 )
		return 0 ;

	double lf0 = LogGammaDensity( x, k[0], theta[0] ) ;
	double lf1 = LogGammaDensity( x, k[1], theta[1] ) ; 
	return (double)1.0 / ( 1.0 + exp( lf1 + log( 1 - pi ) - lf0 - log( pi ) ) ) ;
}

//****************************************************************************80

double alnorm ( double x, bool upper )

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
	double a1 = 5.75885480458;
	double a2 = 2.62433121679;
	double a3 = 5.92885724438;
	double b1 = -29.8213557807;
	double b2 = 48.6959930692;
	double c1 = -0.000000038052;
	double c2 = 0.000398064794;
	double c3 = -0.151679116635;
	double c4 = 4.8385912808;
	double c5 = 0.742380924027;
	double c6 = 3.99019417011;
	double con = 1.28;
	double d1 = 1.00000615302;
	double d2 = 1.98615381364;
	double d3 = 5.29330324926;
	double d4 = -15.1508972451;
	double d5 = 30.789933034;
	double ltone = 7.0;
	double p = 0.398942280444;
	double q = 0.39990348504;
	double r = 0.398942280385;
	bool up;
	double utzero = 18.66;
	double value;
	double y;
	double z;

	up = upper;
	z = x;

	if ( z < 0.0 )
	{
		up = !up;
		z = - z;
	}

	if ( ltone < z && ( ( !up ) || utzero < z ) )
	{
		if ( up )
		{
			value = 0.0;
		}
		else
		{
			value = 1.0;
		}
		return value;
	}

	y = 0.5 * z * z;

	if ( z <= con )
	{
		value = 0.5 - z * ( p - q * y 
				/ ( y + a1 + b1 
				/ ( y + a2 + b2 
				/ ( y + a3 ))));
	}
	else
	{
		value = r * exp ( - y ) 
			/ ( z + c1 + d1 
			/ ( z + c2 + d2 
			/ ( z + c3 + d3 
			/ ( z + c4 + d4 
			/ ( z + c5 + d5 
			/ ( z + c6 ))))));
	}

	if ( !up )
	{
		value = 1.0 - value;
	}

	return value;
}

//****************************************************************************80

double gammad ( double x, double p, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMAD computes the Incomplete Gamma Integral
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by B Shea.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    B Shea,
//    Algorithm AS 239:
//    Chi-squared and Incomplete Gamma Integral,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, double X, P, the parameters of the incomplete 
//    gamma ratio.  0 <= X, and 0 < P.
//
//    Output, int IFAULT, error flag.
//    0, no error.
//    1, X < 0 or P <= 0.
//
//    Output, double GAMMAD, the value of the incomplete 
//    Gamma integral.
//
{
	double a;
	double an;
	double arg;
	double b;
	double c;
	double elimit = - 88.0;
	double oflo = 1.0E+37;
	double plimit = 1000.0;
	double pn1;
	double pn2;
	double pn3;
	double pn4;
	double pn5;
	double pn6;
	double rn;
	double tol = 1.0E-14;
	bool upper;
	double value;
	double xbig = 1.0E+08;

	value = 0.0;
	//
	//  Check the input.
	//
	if ( x < 0.0 )
	{
		*ifault = 1;
		return value;
	}

	if ( p <= 0.0 )
	{
		*ifault = 1;
		return value;
	}

	*ifault = 0;

	if ( x == 0.0 )
	{
		value = 0.0;
		return value;
	}
	//
	//  If P is large, use a normal approximation.
	//
	if ( plimit < p )
	{
		pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 ) 
				+ 1.0 / ( 9.0 * p ) - 1.0 );

		upper = false;
		value = alnorm ( pn1, upper );
		return value;
	}
	//
	//  If X is large set value = 1.
	//
	if ( xbig < x )
	{
		value = 1.0;
		return value;
	}
	//
	//  Use Pearson's series expansion.
	//  (Note that P is not large enough to force overflow in ALOGAM).
	//  No need to test IFAULT on exit since P > 0.
	//
	if ( x <= 1.0 || x < p )
	{
		arg = p * log ( x ) - x - lgamma ( p + 1.0 );
		c = 1.0;
		value = 1.0;
		a = p;

		for ( ; ; )
		{
			a = a + 1.0;
			c = c * x / a;
			value = value + c;

			if ( c <= tol )
			{
				break;
			}
		}

		arg = arg + log ( value );

		if ( elimit <= arg )
		{
			value = exp ( arg );
		}
		else
		{
			value = 0.0;
		}
	}
	//
	//  Use a continued fraction expansion.
	//
	else 
	{
		arg = p * log ( x ) - x - lgamma ( p );
		a = 1.0 - p;
		b = a + x + 1.0;
		c = 0.0;
		pn1 = 1.0;
		pn2 = x;
		pn3 = x + 1.0;
		pn4 = x * b;
		value = pn3 / pn4;

		for ( ; ; )
		{
			a = a + 1.0;
			b = b + 2.0;
			c = c + 1.0;
			an = a * c;
			pn5 = b * pn3 - an * pn1;
			pn6 = b * pn4 - an * pn2;

			if ( pn6 != 0.0 )
			{
				rn = pn5 / pn6;

				if ( fabs ( value - rn ) <= r8_min ( tol, tol * rn ) )
				{
					break;
				}
				value = rn;
			}

			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
			//
			//  Re-scale terms in continued fraction if terms are large.
			//
			if ( oflo <= fabs ( pn5 ) )
			{
				pn1 = pn1 / oflo;
				pn2 = pn2 / oflo;
				pn3 = pn3 / oflo;
				pn4 = pn4 / oflo;
			}
		}

		arg = arg + log ( value );

		if ( elimit <= arg )
		{
			value = 1.0 - exp ( arg );
		}
		else
		{
			value = 1.0;
		}
	}

	return value;
}

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}

// This function is implemented by Li Song. 
// If Y follows (theta, k)=(1/alpha, k)-gamma distribution, then P_Y(Y<t)=gammad(alpha*t, k).
// Furthurmore, chi-square with d.f. k is gamma distribution (2, k/2)
double chicdf( double x, double df ) 
{
	int ifault ;
	return gammad( x / 2.0, df / 2.0, &ifault ) ;
}

