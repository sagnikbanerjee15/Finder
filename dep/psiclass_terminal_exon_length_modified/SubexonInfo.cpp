// report the total depth for the segment
// Format: ./a.out alignments.bam splice_site
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <algorithm>
#include <vector>
#include <math.h>

#include "alignments.hpp"
#include "blocks.hpp"
#include "stats.hpp"

#define ABS(x) ((x)<0?-(x):(x))

char usage[] = "./subexon-info alignment.bam intron.splice [options]\n"
		"options:\n"
		"\t--minDepth INT: the minimum coverage depth considered as part of a subexon (default: 2)\n"
		"\t--noStats: do not compute the statistical scores (default: not used)\n" ;
char buffer[4096] ;

int gMinDepth ;

bool CompSplitSite( struct _splitSite a, struct _splitSite b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else if ( a.type != b.type )
		return a.type < b.type ; // We want the start of exons comes first, since we are scanning the genome from left to right,
					 // so we can terminate the extension early and create a single-base subexon later.
	else
		return a.oppositePos < b.oppositePos ;
}

bool CompBlocksByAvgDepth( struct _block a, struct _block b )
{
	double avgA = a.depthSum / (double)( a.end - a.start + 1 ) ;
	double avgB = b.depthSum / (double)( b.end - b.start + 1 ) ;
	return avgA < avgB ;
}

bool CompBlocksByRatio( struct _block a, struct _block b )
{
	return a.ratio < b.ratio ;	
}

int CompDouble( const void *p1, const void *p2 )
{
	double a = *(double *)p1 ;
	double b = *(double *)p2 ;
	if ( a > b )
		return 1 ;
	else if ( a < b )
		return -1 ;
	else
		return 0 ;
}

// Clean up the split sites;
void FilterAndSortSplitSites( std::vector<struct _splitSite> &sites )
{
	std::sort( sites.begin(), sites.end(), CompSplitSite ) ;	
	int i, j, k, l ;
	int size = sites.size() ;

	for ( i = 0 ; i < size ; )
	{
		for ( j = i + 1 ; j < size ; ++j )
			if ( sites[j].chrId != sites[i].chrId || sites[j].pos != sites[i].pos || sites[j].type != sites[i].type )
				break ;
		int maxSupport = 0 ;
		for ( k = i ; k < j ; ++k )
			if ( sites[k].support > maxSupport )
				maxSupport = sites[k].support ;
		
		int strandCnt[2] = {0, 0} ;
		char strand = sites[i].strand ;
		for ( k = i ; k < j ; ++k )
		{
			if ( sites[k].strand == '-' )
				strandCnt[0] += sites[k].support ; 
			else if ( sites[k].strand == '+' )
				strandCnt[1] += sites[k].support ;
			
		}
		if ( strandCnt[0] > strandCnt[1] )
			strand = '-' ;
		else if ( strandCnt[1] > strandCnt[0] )
			strand = '+' ;

		bool allOneExceptMax = false ;
		if ( maxSupport >= 20 )
		{
			allOneExceptMax = true ;
			for ( k = i ; k < j ; ++k )	
			{
				if ( sites[k].support == maxSupport )	
					continue ;
				if ( sites[k].support > 1 )
				{
					allOneExceptMax = false ;
					break ;
				}
			}
		}
		
		for ( k = i ; k < j ; ++k )
		{
			if ( ( sites[k].support < 0.01 * maxSupport && sites[k].support <= 3 ) 
				|| sites[k].support < 0.001 * maxSupport || sites[k].strand != strand // The introns from the different strand are filtered.
				|| ( sites[k].support < 0.02 * maxSupport && sites[k].mismatchSum >= 2 * sites[k].support )
				|| ( allOneExceptMax && sites[k].support == 1 ) 
				|| ( sites[k].support <= 2 && sites[k].mismatchSum >= 2 * sites[k].support )
				|| ( maxSupport >= 2 && sites[k].support == 1 && ( ABS( sites[k].oppositePos - sites[k].pos - 1 ) >= 10000 || sites[k].mismatchSum != 0 ) ) ) 
			{
				for ( l = i - 1 ; l >= 0 && sites[l].chrId == sites[i].chrId && sites[l].pos >= sites[k].oppositePos ; --l )
				{
					if ( sites[l].pos == sites[k].oppositePos && sites[l].oppositePos == sites[k].pos )
					{
						sites[l].support = -1 ;
						sites[k].support = -1 ;
						break ;
					}
				}

				for ( l = j ; l < size && sites[l].chrId == sites[i].chrId && sites[l].pos <= sites[k].oppositePos ; ++l )
				{
					if ( sites[l].pos == sites[k].oppositePos && sites[l].oppositePos == sites[k].pos )
					{
						sites[l].support = -1 ;
						sites[k].support = -1 ;
						break ;
					}
				}
			}
			/*else if ( sites[k].support <= 1 && sites[k].oppositePos - sites[k].pos + 1 >= 30000 )
			{
				for ( l = j ; l < size && sites[l].chrId == sites[i].chrId && sites[l].pos <= sites[k].oppositePos ; ++l )
				{
					if ( sites[l].pos == sites[k].oppositePos && sites[l].oppositePos == sites[k].pos )
					{
						if ( l - j >= 20 )
						{
							sites[l].support = -1 ;
							sites[k].support = -1 ;
							break ;
						}
					}
				}
				
			}*/
		}
		i = j ;
	}

	k = 0 ;
	for ( i = 0 ; i < size ; ++i )
	{
		if ( sites[i].support > 0 )
		{
			sites[k] = sites[i] ;
			if ( sites[k].strand == '?' )
				sites[k].strand = '.' ;
			++k ;
		}
	}
	sites.resize( k ) ;
}

// Remove the same sites.
void KeepUniqSplitSites( std::vector< struct _splitSite> &sites )
{
	int i, j ;
	int size = sites.size() ;
	int k = 0 ;
	for ( i = 0 ; i < size ; )
	{
		for ( j = i + 1 ; j < size ; ++j )
			if ( sites[j].chrId != sites[i].chrId || sites[j].pos != sites[i].pos || sites[j].type != sites[i].type )
				break ;
		sites[k] = sites[i] ;
		++k ;
		/*else
		{
			if ( sites[i].type != sites[k-1].type )
			{
				printf( "%d\n", sites[i].pos ) ;
			}
		}*/
		i = j ;
	}
	sites.resize( k ) ;

	// For the sites that corresponds to the start of an exon, we remove the adjust to it and 
	/*for ( i = 1 ; i < size ; ++i )
	{
		if ( sites[i].pos == sites[i - 1].pos && sites[i - 1].type == 2 )
		{
			if ( sites[i].type == 1 )	
				++sites[i].pos ;
		}
	}*/
}

// Filter split sites that are extremely close to each other but on different strand.
void FilterNearSplitSites( std::vector< struct _splitSite> &sites )
{
	int i, j ;
	int size = sites.size() ;
	int k = 0 ; 
	for ( i = 0 ; i < size - 1 ; ++i )
	{
		if ( sites[i].support < 0 || sites[i].type != sites[i + 1].type || sites[i].chrId != sites[i + 1].chrId )	
			continue ;
		if ( sites[i + 1].pos - sites[i].pos <= 7 && 
			( sites[i + 1].strand != sites[i].strand || 
				sites[i].strand == '?' ) ) 	
		{
			int tag = i ;
			if ( sites[i + 1].support < sites[i].support )
				tag = i + 1 ;
			sites[tag].support = -1 ;
			int direction ;
			if ( sites[tag].oppositePos < sites[tag].pos )
				direction = -1 ;
			else
				direction = 1 ;

			for ( j = tag ; j >= 0 && j < size ; j += direction )
				if ( sites[j].pos == sites[tag].oppositePos && sites[j].oppositePos == sites[tag].pos )
				{
					sites[j].support = -1 ;
					break ;
				}
		}
	}

	for ( i = 0 ; i < size ; ++i )
	{
		if ( sites[i].support > 0 )
		{
			sites[k] = sites[i] ;
			++k ;
		}
	}
	sites.resize( k ) ;
}

void FilterRepeatSplitSites( std::vector<struct _splitSite> &sites )
{
	int i, j ;
	int size = sites.size() ;
	int k = 0 ;
	for ( i = 0 ; i < size ; )
	{
		for ( j = i + 1 ; j < size ; ++j )
		{
			if ( sites[j].pos != sites[i].pos || sites[j].type != sites[i].type || sites[i].chrId != sites[j].chrId )
				break ;	
		}
		int max = -1 ;
		int maxtag = 0 ; 
		for ( k = i ; k < j ; ++k )
		{
			if ( sites[k].uniqSupport > max )
			{
				max = sites[k].uniqSupport ;
				maxtag = k ;
			}
			else if ( sites[k].uniqSupport == max && sites[k].support > sites[maxtag].support )
			{
				maxtag = k ;
			}
		}

		if ( max > -1 )
		{
			if ( sites[maxtag].uniqSupport > sites[maxtag].support * 0.1 )
			{
				for ( k = i ; k < j ; ++k )
					if ( sites[k].uniqSupport < 0.05 * sites[k].support )
					{
						sites[k].support = -1 ;

						int direction ;
						if ( sites[k].oppositePos < sites[k].pos )
							direction = -1 ;
						else
							direction = 1 ;
						int l ;
						for ( l = k ; l >= 0 && l < size ; l += direction )
							if ( sites[l].pos == sites[k].oppositePos && sites[l].oppositePos == sites[k].pos )
							{
								sites[l].support = -1 ;
								break ;
							}
					}
			}
			else 
			{
				for ( k = i ; k < j ; ++k )
				{
					if ( sites[k].support <= 10 )
					{
						sites[k].support = -1 ;

						int direction ;
						if ( sites[k].oppositePos < sites[k].pos )
							direction = -1 ;
						else
							direction = 1 ;
						int l ;
						for ( l = k ; l >= 0 && l < size ; l += direction )
							if ( sites[l].pos == sites[k].oppositePos && sites[l].oppositePos == sites[k].pos )
							{
								sites[l].support = -1 ;
								break ;
							}
					}
				}
			}
		}

		i = j ;
	}

	k = 0 ;
	for ( i = 0 ; i < size ; ++i )
	{
		if ( sites[i].support > 0 )
		{
			sites[k] = sites[i] ;
			++k ;
		}
	}
	sites.resize( k ) ;

}


// When k=1, the gamma distribution becomes exponential distribution, and can be optimized analytically..
// Maximize: sum_i z_i( log( 1/theta e^{-x_i / theta} )
/*double ThetaOfExponentialDistribution( double *x, double *z, int n )
{
	double sumZ = 0;
	double sumZX = 0 ;
	for ( i = 0 ; i < n ; ++i )
	{
		sumZ += z[i] ;
		sumZX += z[i] * x[i] ;
	}
}*/

// for boundK, if it is positive, it represent the upper bound. If it is negative, -boundK will be the lower bound for k.
// if boundK==0, there is no extra bound.
// The same logic for boundProduct, which bounds k*theta
void GradientDescentGammaDistribution( double &k, double &theta, double initK, double initTheta, double lowerBoundK, double upperBoundK, 
	double lowerBoundMean, double upperBoundMean, double *x, double *z, int n ) 
{
	int i ;
	k = initK ;
	theta = initTheta ;
	double c = 0.5 ;
	int iterCnt = 1 ;

	double sumZ = 0 ;
	double sumZX = 0 ;
	double sumZLogX = 0 ;
	double Hessian[2][2] ; // 0 for k, 1 for theta
	double inverseHessian[2][2] ;
	int tmp ;

	double positiveKRecord = -5, positiveThetaRecord = -5 ; // record the value of k, theta when theta, k becomes non-positive.

	for ( i = 0 ; i < n ; ++i )	
	{
		sumZ += z[i] ;
		sumZX += z[i] * x[i]  ;
		sumZLogX += z[i] * log( x[i] ) ;
	}

	while ( 1 )
	{
		double gradK = 0 ;
		double gradTheta = 0 ;

		double prevK = k ;
		double prevTheta = theta ;
		double digammaK = digammal( k ) ;

		gradK = sumZ * ( -log( theta ) - digammaK ) + sumZLogX ;
		gradTheta = -sumZ * ( k / theta ) + sumZX / ( theta * theta ) ;

		Hessian[0][0] = -sumZ * trigamma( k, &tmp ) ;
		Hessian[0][1] = -sumZ / theta ; // \partial l / ( \partial k \partial theta)
		Hessian[1][0] = -sumZ / theta ;
		Hessian[1][1] = sumZ * k / ( theta * theta ) - 2 * sumZX / ( theta * theta * theta ) ;

		double det = Hessian[0][0] * Hessian[1][1] - Hessian[0][1] * Hessian[1][0] ;
		/*printf( "%s iter %d:\n", __func__, iterCnt ) ;
		printf( "%lf %lf %lf %lf\n", sumZ, k, theta, sumZX ) ;	
		printf( "%lf %lf %lf\n", gradK, gradTheta, det ) ;	
		printf( "%lf %lf %lf %lf\n", Hessian[0][0], Hessian[0][1], Hessian[1][0], Hessian[1][1] ) ;*/
		if ( det <= 1e-4 && det >=-1e-4 )
		{
			k = k + c / iterCnt * gradK ;
			theta = theta + c / iterCnt * gradTheta ;
		}
		else
		{
			inverseHessian[0][0] = Hessian[1][1] / det ;
			inverseHessian[0][1] = -Hessian[0][1] / det ;
			inverseHessian[1][0] = -Hessian[1][0] / det ;
			inverseHessian[1][1] = Hessian[0][0] / det ;
			//printf( "%lf %lf %lf %lf: %lf\n=====\n", inverseHessian[0][0], inverseHessian[0][1], inverseHessian[1][0], inverseHessian[1][1],
			//	Hessian[1][0] * inverseHessian[0][1] + Hessian[1][1] * inverseHessian[1][1] ) ;	
			double step = 0.5 ;
			k = k - step * ( inverseHessian[0][0] * gradK + inverseHessian[0][1] * gradTheta ) ;
			theta = theta - step * ( inverseHessian[1][0] * gradK + inverseHessian[1][1] * gradTheta ) ;
			
			bool flag = false ;
			if ( k <= 1e-6 )
			{
				step = ( prevK - 1e-6 ) / ( inverseHessian[0][0] * gradK + inverseHessian[0][1] * gradTheta ) ;

				if ( ABS( theta - positiveThetaRecord ) < 1e-5 )
					flag = true ;
				positiveThetaRecord = theta ;
			}
			if ( theta <= 1e-6 )
			{
				double tmp = ( prevTheta - 1e-6 ) / ( inverseHessian[1][0] * gradK + inverseHessian[1][1] * gradTheta ) ;
				if ( tmp < step )
					step = tmp ;

				if ( ABS( k - positiveKRecord ) < 1e-5 )
					flag = true ;
				positiveKRecord = k ;
			}

			if ( step != 0.5 )
			{
				k = prevK - step * ( inverseHessian[0][0] * gradK + inverseHessian[0][1] * gradTheta ) ;
				theta = prevTheta - step * ( inverseHessian[1][0] * gradK + inverseHessian[1][1] * gradTheta ) ;
			}

			/*if ( flag ) 
			{
				k = prevK ;
				theta = prevTheta ;
				break ;
			}*/
		}

		if ( upperBoundK > 0 && k > upperBoundK )
			k = upperBoundK ;
		else if ( lowerBoundK > 0 && k < lowerBoundK )
			k = lowerBoundK ;

		if ( upperBoundMean > 0 && k * theta > upperBoundMean )
		{
			theta = upperBoundMean / k ;
		}
		else if ( lowerBoundMean > 0 && k * theta < lowerBoundMean )
		{
			theta = lowerBoundMean / k ;
		}

		if ( k <= 1e-6 )
		{
			k = 1e-6 ;
		}

		if ( theta <= 1e-6 )
		{
			theta = 1e-6 ;
		}

		double diff = ABS( prevK - k ) + ABS( prevTheta - theta ) ;
		if ( diff < 1e-5 ) 
			break ;
		
		++iterCnt ;
		//if ( det <= 1e-4 && det >=-1e-4 && k >= 5000 ) //&& diff < 1 )
		if ( k >= 10000 ) //&& diff < 1 )
		{
			k = prevK ;
			theta = prevTheta ;
			break ;
		}
		if ( iterCnt == 1000 )
			break ;
	}
}

double MixtureGammaEM( double *x, int n, double &pi, double *k, double *theta, int tries, double meanBound[2], int iter = 1000 )
{
	int i ;
	if ( n <= 0 )	
		return 0 ;
	double *z = new double[n] ; // the expectation that it assigned to model 0.
	double *oneMinusZ = new double[n] ;
	int t = 0 ;
	double history[5] = {-1, -1, -1, -1, -1} ;
	double maxX = -1 ;
	double sumX = 0 ;

	for ( i = 0 ; i < n ; ++i )
	{
		sumX += x[i] ;
		if ( x[i] > maxX )
			maxX = x[i] ;
	}
	if ( maxX > meanBound[1] && meanBound[1] >= 0 )
		maxX = meanBound[1] ;

	/*if ( meanBound[1] == -1 )
	{
		// The EM for coverage
		maxX = 10.0 ;
	}*/
		
	while ( 1 )
	{
		double npi, nk[2], ntheta[2] ;
		double sum = 0 ;
		for ( i = 0 ; i < n ; ++i )	
		{
			//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ); 
			//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ); 
			//z[i] = exp( lf0 + log( pi ) ) / ( exp( lf0 + log( pi ) ) + exp( lf1 + log( 1 - pi ) ) ) ;
			if ( pi != 0 )
				z[i] = MixtureGammaAssignment( x[i], pi, k, theta ) ;
			else
				z[i] = 0 ;
			/*if ( isnan( z[i] ) )
			{
				printf( "nan: %lf %lf %lf %lf\n", x[i], pi, k, theta ) ;
			}*/
			oneMinusZ[i] = 1 - z[i] ;
			sum += z[i] ;
		}

		// compute new pi.
		npi = sum / n ;

		// Use gradient descent to compute new k and theta.
		if ( 1 ) //pi > 0 )
		{
			double bound ;
			if ( meanBound[1] != -1 )  // the EM for ratio
			{
				bound = ( theta[1] * k[1] > 1 ) ? 1 : ( theta[1] * k[1] ) / ( 1 + tries );
				GradientDescentGammaDistribution( nk[0], ntheta[0], k[0], theta[0], k[1], -1, -1, bound, x, z, n ) ; // It seems setting an upper bound 1 for k[0] is not a good idea.
			}
			else
			{
				bound = ( theta[1] * k[1] > 1 ) ? 1 : ( theta[1] * k[1] ) / ( 1 + tries ) ;
				GradientDescentGammaDistribution( nk[0], ntheta[0], k[0], theta[0], k[1], -1, meanBound[0], bound, x, z, n ) ; // It seems setting an upper bound 1 for k[0] is not a good idea.
			}
			GradientDescentGammaDistribution( nk[1], ntheta[1], k[1], theta[1], -1, k[0], theta[0] * k[0], maxX, x, oneMinusZ,  n ) ;
		}
		else
		{
			GradientDescentGammaDistribution( nk[1], ntheta[1], k[1], theta[1], 0, 0, 0, 0, x, oneMinusZ,  n ) ;
		}

		double diff ;
		if ( isnan( npi ) || isnan( nk[0] ) || isnan( nk[1] ) || isnan( ntheta[0] ) || isnan( ntheta[1] ) )
		{
			delete[] z ;
			delete[] oneMinusZ ;
			return -1 ;
		}
		diff = ABS( nk[0] - k[0] ) + ABS( nk[1] - k[1] )
			+ ABS( ntheta[0] - theta[0] ) + ABS( ntheta[1] - theta[1] ) ; // pi is fully determined by these 4 parameters.
		if ( diff < 1e-4 )
			break ;
		diff = ABS( nk[0] - history[1] ) + ABS( nk[1] - history[2] )
			+ ABS( ntheta[0] - history[3]  ) + ABS( ntheta[1] - history[4] ) ; // pi is fully determined by these 4 parameters.
		if ( diff < 1e-4 )
			break ;

		history[0] = pi ;
		history[1] = k[0] ;
		history[2] = k[1] ;
		history[3] = theta[0] ;
		history[4] = theta[1] ;
		
		pi = npi ;
		k[0] = nk[0] ;
		k[1] = nk[1] ;
		theta[0] = ntheta[0] ;
		theta[1] = ntheta[1] ;

		/*double logLikelihood = 0 ;
		for ( i = 0 ; i < n ; ++i )
			logLikelihood += log( pi * exp( LogGammaDensity( x[i], k[0], theta[0]) ) + 
				(1 - pi ) * exp( LogGammaDensity( x[i], k[1], theta[1] ) ) ) ;*/

		//printf( "%d: %lf %lf %lf %lf %lf\n", t, pi, k[0], theta[0], k[1], theta[1] ) ;
		
		++t ;
		if ( iter != -1 && t >= iter )
			break ;
	}
	delete[] z ;
	delete[] oneMinusZ ;
	return 0 ;
}

bool IsParametersTheSame( double *k, double *theta )
{
	if ( ABS( k[0] - k[1] ) < 1e-2 && ABS( theta[0] - theta[1] ) < 1e-2 )	
		return true ;
	return false ;
}

int RatioAndCovEM( double *covRatio, double *cov, int n, double &piRatio, double kRatio[2], 
	double thetaRatio[2], double &piCov, double kCov[2], double thetaCov[2] )
{
	int i ;
	piRatio = 0.6 ; // mixture coefficient for model 0 and 1
	kRatio[0] = 0.9 ;
	kRatio[1] = 0.45 ;
	thetaRatio[0] = 0.05 ;
	thetaRatio[1] = 1 ;
	double meanBound[2] = {-1, 1} ; // [0] is for the lower bound of the noise model, [1] is for the upper bound of the true model
	
	/*double *filteredCovRatio = new double[n] ;// ignore the ratio that is greater than 5.
	int m = 0 ;
	for ( i = 0 ; i < n ; ++i )
		if ( covRatio[i] < 1.0 )
		{
			filteredCovRatio[m] = covRatio[i] ;
			++m ;
		}*/
	srand( 17 ) ;
	int maxTries = 10 ;
	int t = 0 ;
	double *buffer = new double[n] ;
	for ( i = 0 ; i < n ; ++i )
		buffer[i] = covRatio[i] ;
	qsort( buffer, n, sizeof( double ), CompDouble ) ;
	//covRatio = buffer ;
	while ( 1 )
	{
		//printf( "EM\n" )  ;
		MixtureGammaEM( covRatio, n, piRatio, kRatio, thetaRatio, t, meanBound ) ;
		//printf( "%lf %lf %lf %lf %lf\n", piRatio, kRatio[0], kRatio[1], thetaRatio[0], thetaRatio[1] ) ;
		if ( piRatio > 0.999 || piRatio < 0.001 || IsParametersTheSame( kRatio, thetaRatio ) )
		{
			++t ;
			if ( t > maxTries )
				break ;
			piRatio = 0.6 ;
			kRatio[0] += ( ( rand() * 0.5 - RAND_MAX ) / (double)RAND_MAX * 0.1 ) ;
			if ( kRatio[0] <= 0 )
				kRatio[0] = 0.9 ;
			kRatio[1] += ( ( rand() * 0.5 - RAND_MAX ) / (double)RAND_MAX * 0.1 ) ;
			if ( kRatio[1] <= 0 )
				kRatio[1] = 0.45 ;
			thetaRatio[0] += ( ( rand() * 0.5 - RAND_MAX ) / (double)RAND_MAX * 0.1 ) ;
			if ( thetaRatio[0] <= 0 )
				thetaRatio[0] = 0.05 ;
			thetaRatio[1] += ( ( rand() * 0.5 - RAND_MAX ) / (double)RAND_MAX * 0.1 ) ;
			if ( thetaRatio[1] <= 0 )
				thetaRatio[1] = 1 ;
			if ( kRatio[0] < kRatio[1] )
			{	
				if ( rand() & 1 )
					kRatio[0] = kRatio[1] ;
				else
					kRatio[1] = kRatio[0] ;
			}
			if ( kRatio[0] * thetaRatio[0] > kRatio[1] * thetaRatio[1] )
			{
				thetaRatio[0] = kRatio[1] * thetaRatio[1] / kRatio[0] ;
			}
			//printf( "%lf %lf %lf %lf %lf\n", piRatio, kRatio[0], kRatio[1], thetaRatio[0], thetaRatio[1] ) ;

			continue ;
		}
		
		break ;
	}
	//delete[] filteredCovRatio ;
	if ( t > maxTries && piRatio > 0.999 )
	{
		/*piRatio = 0.6 ; // mixture coefficient for model 0 and 1
		kRatio[0] = 0.9 ;
		kRatio[1] = 0.45 ;
		thetaRatio[0] = 0.05 ;
		thetaRatio[1] = 1 ;*/
		piRatio = 0.999 ;
	}
	if ( IsParametersTheSame( kRatio, thetaRatio ) || piRatio <= 1e-3 )
		piRatio = 1e-3 ;	
	
	piCov = piRatio ; // mixture coefficient for model 0 and 1
	kCov[0] = 0.9 ;
	kCov[1] = 0.45 ;
	thetaCov[0] = 3 ;
	thetaCov[1] = 12 ;
	
	// only do one iteration of EM, so that pi does not change?
	// But it seems it still better to run full EM.
	meanBound[0] = 1.01 ;
	meanBound[1] = -1 ;

	//printf( "for coverage:\n" ) ;
	//piCov = 0.001000 ;
	MixtureGammaEM( cov, n, piCov, kCov, thetaCov, 0, meanBound ) ;	
	//printf( "for coverage done\n" ) ;
	piCov = piRatio ;

	delete []buffer ;

	return 0 ;
}

double GetPValue( double x, double *k, double *theta )
{
	int fault ;
	double p ;
	p = 1 - gammad( x / theta[0], k[0], &fault ) ;
	return p ;
}

// if x's value is less than the average of (k0-1)*theta0, then we force x=(k0-1)*theta0,
//    the mode of the model 0. Of course, it does not affect when k0<=1 already.
double MixtureGammaAssignmentAdjust( double x, double pi, double* k, double *theta ) 
{
	if ( x < ( k[0] - 1 ) * theta[0] )
	{
		x = ( k[0] - 1 ) * theta[0] ;
	}
	return MixtureGammaAssignment( x, pi, k, theta ) ;
}


// Transform the cov number for better fitting 
double TransformCov( double c )
{
	double ret ;
	// original it is c-1.
	// Use -2 instead of -1 is that many region covered to 1 reads will be filtered when
	// build the subexons.
	//
	//ret = sqrt( c ) - 1 ;
	if ( c <= 2 + 1e-6 )
		ret = 1e-6 ;
	else
		ret = c - 2 ;
	return ret ;
	//return log( c ) / log( 2.0 ) ;
}

int main( int argc, char *argv[] )
{
	int i, j ;
	bool noStats = false ;
	if ( argc < 3 )
	{
		fprintf( stderr, usage ) ;
		exit( 1 ) ;
	}

	gMinDepth = 2 ;

	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "--noStats" ) )
		{
			noStats = true ;
			continue ;
		}
		else if ( !strcmp( argv[i], "--minDepth" ) )
		{
			gMinDepth = atoi( argv[i + 1] ) ;
			++i ;
			continue ;
		}
		else
		{
			fprintf( stderr, "Unknown argument: %s\n", argv[i] ) ;
			return 0 ;
		}
	}

	Alignments alignments ;
	alignments.Open( argv[1] ) ;
	std::vector<struct _splitSite> splitSites ; // only compromised the 
	std::vector<struct _splitSite> allSplitSites ;

	// read in the splice site
	FILE *fp ;
	fp = fopen( argv[2], "r" ) ;
	char chrom[50] ;
	int64_t start, end ;
	int support ;
	char strand[3] ;
	int uniqSupport, secondarySupport, uniqEditDistance, secondaryEditDistance ;
	while ( fscanf( fp, "%s %" PRId64 " %" PRId64 " %d %s %d %d %d %d", chrom, &start, &end, &support, strand, 
				&uniqSupport, &secondarySupport, &uniqEditDistance, &secondaryEditDistance ) != EOF )
	{
		if ( support <= 0 )
			continue ;
		//if ( !( uniqSupport >= 1 
		//	|| secondarySupport > 10 ) )
		//if ( uniqSupport <= 0.01 * ( uniqSupport + secondarySupport ) || ( uniqSupport == 0 && secondarySupport < 20 ) )
		//if ( uniqSupport == 0 && secondarySupport <= 10 )
		//	continue ;
		int chrId = alignments.GetChromIdFromName( chrom ) ; 
		struct _splitSite ss ;
		--start ;
		--end ;
		ss.pos = start ;
		ss.chrId = chrId ;
		ss.type = 2 ;
		ss.oppositePos = end ;
		ss.strand = strand[0] ;
		ss.support = support ;
		ss.uniqSupport = uniqSupport ;
		ss.mismatchSum = uniqEditDistance + secondaryEditDistance ;
		splitSites.push_back( ss ) ;

		ss.pos = end ; 
		ss.type = 1 ;
		ss.oppositePos = start ;
		ss.strand = strand[0] ;
		ss.support = support ;
		ss.uniqSupport = uniqSupport ;
		ss.mismatchSum = uniqEditDistance + secondaryEditDistance ;
		splitSites.push_back( ss ) ;
	}
	fclose( fp ) ;
	//printf( "ss:%d\n", splitSites.size() ) ;
	
	//printf( "ss:%d\n", splitSites.size() ) ;
	
	alignments.GetGeneralInfo( true ) ;
	// Build the blocks
	Blocks regions ;
	alignments.Rewind() ;
	regions.BuildExonBlocks( alignments ) ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;
	
	FilterAndSortSplitSites( splitSites ) ; 
	FilterNearSplitSites( splitSites ) ;
	FilterRepeatSplitSites( splitSites ) ;
	regions.FilterSplitSitesInRegions( splitSites ) ;
	regions.FilterGeneMergeSplitSites( splitSites ) ;


	allSplitSites = splitSites ;
	KeepUniqSplitSites( splitSites ) ;
	
	//for ( i = 0 ; i < splitSites.size() ; ++i )
	//	printf( "%d %d\n", splitSites[i].pos + 1, splitSites[i].oppositePos + 1 ) ;
	// Split the blocks using split site
	regions.SplitBlocks( alignments, splitSites ) ;
	//printf( "%d\n", regions.exonBlocks.size() ) ;
	/*for ( i = 0 ; i < regions.exonBlocks.size() ; ++i )
	{
		struct _block &e = regions.exonBlocks[i] ;
		printf( "%s %" PRId64 " %" PRId64 " %d %d\n", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1,  e.leftType, e.rightType ) ;
	}
	return 0 ;*/
	// Recompute the coverage for each block. 
	alignments.Rewind() ;
	//printf( "Before computeDepth: %d\n", regions.exonBlocks.size() ) ;

	regions.ComputeDepth( alignments ) ;
	//printf( "After computeDepth: %d\n", regions.exonBlocks.size() ) ;

	// Merge blocks that may have a hollow coverage by accident.
	regions.MergeNearBlocks() ;
	//printf( "After merge: %d\n", regions.exonBlocks.size() ) ;
	
	// Put the intron informations
	regions.AddIntronInformation( allSplitSites, alignments ) ;
	//printf( "After add information.\n" ) ;

	// Compute the average ratio against the left and right connected subexons.
	regions.ComputeRatios() ;
	//printf( "After compute ratios.\n" ) ;	
	
	//printf( "Finish building regions.\n" ) ;	
	if ( noStats ) 
	{ 
		// just output the subexons.
		if ( realpath( argv[1], buffer ) == NULL )
		{
			strcpy( buffer, argv[1] ) ;
		}
		printf( "#%s\n", buffer ) ;
		printf( "#fitted_ir_parameter_ratio: pi: -1 k0: -1 theta0: -1 k1: -1 theta1: -1\n" ) ;
		printf( "#fitted_ir_parameter_cov: pi: -1 k0: -1 theta0: -1 k1: -1 theta1: -1\n" ) ;
		
		int blockCnt = regions.exonBlocks.size() ;
		for ( int i = 0 ; i < blockCnt ; ++i )
		{
			struct _block &e = regions.exonBlocks[i] ;
			double avgDepth = (double)e.depthSum / ( e.end - e.start + 1 ) ;
			printf( "%s %" PRId64 " %" PRId64 " %d %d %lf -1 -1 -1 -1 ", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.leftType, e.rightType, avgDepth ) ;
			int prevCnt = e.prevCnt ;
			if ( i > 0 && e.start == regions.exonBlocks[i - 1].end + 1 &&
					e.leftType == regions.exonBlocks[i - 1].rightType )
			{
				printf( "%d ", prevCnt + 1 ) ;
				for ( j = 0 ; j < prevCnt ; ++j )
					printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
				printf( "%" PRId64 " ", regions.exonBlocks[i - 1].end + 1 ) ;
			}
			else
			{
				printf( "%d ", prevCnt ) ;
				for ( j = 0 ; j < prevCnt ; ++j )
					printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
			}

			int nextCnt = e.nextCnt ;
			if ( i < blockCnt - 1 && e.end == regions.exonBlocks[i + 1].start - 1 &&
					e.rightType == regions.exonBlocks[i + 1].leftType )
			{
				printf( "%d %" PRId64 " ", nextCnt + 1, regions.exonBlocks[i + 1].start + 1 ) ;
			}
			else
				printf( "%d ", nextCnt ) ;
			for ( j = 0 ; j < nextCnt ; ++j )
				printf( "%" PRId64 " ", regions.exonBlocks[ e.next[j] ].start + 1 ) ;
			printf( "\n" ) ;

		}
		return 0 ;
	}

	// Extract the blocks for different events.
	int blockCnt = regions.exonBlocks.size() ;
	std::vector<struct _block> irBlocks ; // The regions corresponds to intron retention events.
	double *leftClassifier = new double[ blockCnt ] ; 
	double *rightClassifier = new double[ blockCnt ] ;
	std::vector<struct _block> overhangBlocks ; //blocks like (...[ or ]...)
	std::vector<struct _block> islandBlocks ; // blocks like (....) 

	for ( i = 0 ; i < blockCnt ; ++i )	
	{
		int ltype = regions.exonBlocks[i].leftType ;
		int rtype = regions.exonBlocks[i].rightType ;
		leftClassifier[i] = -1 ;
		rightClassifier[i] = -1 ;
		
		//double avgDepth = (double)regions.exonBlocks[i].depthSum / ( regions.exonBlocks[i].end - regions.exonBlocks[i].start + 1 ) ;
		
		if ( ltype == 2 && rtype == 1 )
		{
			// candidate intron retention.
			// Note that when I compute the ratio, it is already made sure that the avgDepth>1.
			double ratio = regions.PickLeftAndRightRatio( regions.exonBlocks[i] ) ;

			//printf( "%lf %lf\n", regions.exonBlocks[i].leftRatio, regions.exonBlocks[i].rightRatio ) ;
			if ( ratio > 0 )
			{
				regions.exonBlocks[i].ratio = ratio ;
				irBlocks.push_back( regions.exonBlocks[i] ) ;
				irBlocks[ irBlocks.size() - 1 ].contigId = i ;
			}
		}
		else if ( ( ltype == 0 && rtype == 1 ) || ( ltype == 2 && rtype == 0 ) )
		{
			// subexons like (...[ or ]...)
			double ratio ;
			if ( ltype == 0 )
			{
				ratio = regions.exonBlocks[i].rightRatio ;
			}
			else if ( ltype == 2 )
			{
				ratio = regions.exonBlocks[i].leftRatio ;	
			}
			if ( ratio > 0 )
			{
				regions.exonBlocks[i].ratio = ratio ;
				overhangBlocks.push_back( regions.exonBlocks[i] ) ;
				overhangBlocks[ overhangBlocks.size() - 1].contigId = i ;
			}
		}
		else if ( ltype == 0 && rtype == 0 )
		{
			islandBlocks.push_back( regions.exonBlocks[i] ) ;
			islandBlocks[ islandBlocks.size() - 1].contigId = i ;
		}
	}
	
	// Compute the histogram for each intron.
	int irBlockCnt = irBlocks.size() ;
	double *cov = new double[irBlockCnt] ;
	double *covRatio = new double[ irBlockCnt ] ;
	double piRatio, kRatio[2], thetaRatio[2] ;
	double piCov, kCov[2], thetaCov[2] ;
	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		double avgDepth = regions.GetAvgDepth( irBlocks[i] ) ;
		//cov[i] = ( avgDepth - 1 ) / ( flankingAvg - 1 ) ;
		cov[i] = TransformCov( avgDepth ) ;

		covRatio[i] = regions.PickLeftAndRightRatio( irBlocks[i] ) ; 
		//cov[i] = avgDepth / anchorAvg ;
		//printf( "%"PRId64" %d %d: %lf %lf\n", irBlocks[i].depthSum, irBlocks[i].start, irBlocks[i].end, avgDepth, cov[i] ) ;
	}

	/*fp = fopen( "ratio.out", "r" ) ;
	int irBlockCnt = 0 ;
	double *cov = new double[10000] ;
	while ( 1 ) 
	{
		double r ;
		if ( fscanf( fp, "%lf", &r ) == EOF )
			break ;
		cov[ irBlockCnt ] = r ;
		++irBlockCnt ;
	}
	fclose( fp ) ;*/
	RatioAndCovEM( covRatio, cov, irBlockCnt, piRatio, kRatio, thetaRatio, piCov, kCov, thetaCov ) ;

	for ( i = 0 ; i < irBlockCnt ; ++i )
	{
		//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ) ;
		//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ) ;

		double p1, p2, p ;
		
		p1 = MixtureGammaAssignmentAdjust( covRatio[i], piRatio, kRatio, thetaRatio ) ;
		p2 = MixtureGammaAssignmentAdjust( cov[i], piCov, kCov, thetaCov  ) ;

		/*p1 = GetPValue( covRatio[i], kRatio, thetaRatio ) ; //1 - gammad( covRatio[i] / thetaRatio[0], kRatio[0], &fault ) ;
		if ( piRatio != 0 )
			p2 = GetPValue( cov[i], kCov, thetaCov ) ;//1 - gammad( cov[i] / thetaCov[0], kCov[0], &fault ) ;
		else
			p2 = p1 ;*/
		//printf( "%lf %lf: %lf %lf\n", covRatio[i], cov[i], p1, p2 ) ;
		p = p1 > p2 ? p1 : p2 ;
		
		
		//printf( "%d %d: avg: %lf ratio: %lf p: %lf\n", irBlocks[i].start, irBlocks[i].end, irBlocks[i].depthSum / (double)( irBlocks[i].end - irBlocks[i].start + 1 ), covRatio[i], 
		//		p ) ;	
		leftClassifier[ irBlocks[i].contigId ] = p ;
		rightClassifier[ irBlocks[i].contigId ] = p ;
	}

	// Process the classifier for overhang subexons and the subexons to see whether we need soft boundary
	int overhangBlockCnt = overhangBlocks.size() ;
	delete []cov ;
	delete []covRatio ;
	
	cov = new double[ overhangBlockCnt ] ;
	covRatio = new double[overhangBlockCnt] ;
	double overhangPiRatio, overhangKRatio[2], overhangThetaRatio[2] ;
	double overhangPiCov, overhangKCov[2], overhangThetaCov[2] ;

	for ( i = 0 ; i < overhangBlockCnt ; ++i )	
	{
		covRatio[i] = overhangBlocks[i].ratio ;
		cov[i] = TransformCov( regions.GetAvgDepth( overhangBlocks[i] ) ) ;
	}
	RatioAndCovEM( covRatio, cov, overhangBlockCnt, overhangPiRatio, overhangKRatio, overhangThetaRatio, overhangPiCov, overhangKCov, overhangThetaCov ) ;

	for ( i = 0 ; i < overhangBlockCnt ; ++i )
	{
		//double lf0 = -k[0] * log( theta[0] ) + ( k[0] - 1 ) * log( cov[i]) - cov[i] / theta[0] - lgamma( k[0] ) ;
		//double lf1 = -k[1] * log( theta[1] ) + ( k[1] - 1 ) * log( cov[i]) - cov[i] / theta[1] - lgamma( k[1] ) ;

		double p1, p2, p ;
		p1 = MixtureGammaAssignmentAdjust( covRatio[i], overhangPiRatio, overhangKRatio, overhangThetaRatio ) ;
		p2 = MixtureGammaAssignmentAdjust( cov[i], overhangPiCov, overhangKCov, overhangThetaCov  ) ;
			
		/*p1 = GetPValue( covRatio[i], overhangKRatio, overhangThetaRatio ) ; //1 - gammad( covRatio[i] / thetaRatio[0], kRatio[0], &fault ) ;
		if ( overhangPiRatio != 0)
			p2 = GetPValue( cov[i], overhangKCov, overhangThetaCov ) ;//1 - gammad( cov[i] / thetaCov[0], kCov[0], &fault ) ;
		else
			p2 = p1 ;*/

		//p = p1 < p2 ? p1 : p2 ;
		p = sqrt( p1 * p2 ) ;
		
		int idx = overhangBlocks[i].contigId ;
		if ( regions.exonBlocks[idx].rightType == 0 )
			leftClassifier[ idx ] = rightClassifier[ idx ] = p ;
		else
			leftClassifier[ idx ] = rightClassifier[ idx ] = p ;
	}

	for ( i = 0 ; i < blockCnt ; ++i )		
	{
		struct _block &e = regions.exonBlocks[i] ;
		//if ( ( e.leftType == 0 && e.rightType == 1 ) || 
		// variance-stabailizing transformation of poisson distribution. But we are more conservative here.
		// The multiply 2 before that is because we ignore the region below 0, so we need to somehow renormalize the distribution.
		if ( e.leftType == 1 )
		{
			if ( e.leftRatio >= 0 )
				leftClassifier[i] = 2 * alnorm( e.leftRatio * 2.0 , true ) ;
			else
				leftClassifier[i] = 1 ;
		}
		/*else if ( e.leftType == 0 )
		{
			// If this region is a start of a gene, the other sample might introduce
			//    a new intron before it. So we want to test whether this region can still
			//    be a start of a gene even there is an intron before it.
			for ( j = i + 1 ; j < blockCnt ; ++j )
			{
				if ( regions.exonBlocks[j].contigId != regions.exonBlocks[i].contigId )
					break ;
			}

			for ( k = i ; k < j ; ++k )
				if ( regions.exonBlocks[j].leftType == 1 )
					break ;
			if ( k >= j )
			{
				leftClassifier[i] = alnorm( )
			}
		}*/

		//if ( ( e.rightType == 0 && e.leftType == 2 ) || 
		if ( e.rightType == 2 )
		{
			if ( e.rightRatio >= 0 )
				rightClassifier[i] = 2 * alnorm( e.rightRatio * 2.0, true ) ;
			else
				rightClassifier[i] = 1 ;
		}
	}

	// Process the result for subexons seems like single-exon transcript (...)
	int islandBlockCnt = islandBlocks.size() ;
	//std::sort( islandBlocks.begin(), islandBlocks.end(), CompBlocksByAvgDepth ) ;
	for ( i = 0, j = 0 ; i < islandBlockCnt ; ++i )
	{
		/*if ( regions.GetAvgDepth( islandBlocks[i] ) != regions.GetAvgDepth( islandBlocks[j] ) )
			j = i ;
		leftClassifier[ islandBlocks[i].contigId ] = 1 - (j + 1) / (double)( islandBlockCnt + 1 ) ;
		rightClassifier[ islandBlocks[i].contigId ] = 1 - (j + 1) / (double)( islandBlockCnt + 1 ) ;*/
		double p = GetPValue( TransformCov( regions.GetAvgDepth( islandBlocks[i] ) ), kCov, thetaCov ) ;
		leftClassifier[ islandBlocks[i].contigId ] = p ;
		rightClassifier[ islandBlocks[i].contigId ] = p ;
	}

	
	// Output the result.
	if ( realpath( argv[1], buffer ) == NULL )
	{
		strcpy( buffer, argv[1] ) ;
	}
	printf( "#%s\n", buffer ) ;
	// TODO: higher precision.
	printf( "#fitted_ir_parameter_ratio: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piRatio, kRatio[0], thetaRatio[0], kRatio[1], thetaRatio[1] ) ;
	printf( "#fitted_ir_parameter_cov: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", piCov, kCov[0], thetaCov[0], kCov[1], thetaCov[1] ) ;
	
	printf( "#fitted_overhang_parameter_ratio: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", overhangPiRatio, overhangKRatio[0], overhangThetaRatio[0], overhangKRatio[1], overhangThetaRatio[1] ) ;
	printf( "#fitted_overhang_parameter_cov: pi: %lf k0: %lf theta0: %lf k1: %lf theta1: %lf\n", overhangPiCov, overhangKCov[0], overhangThetaCov[0], overhangKCov[1], overhangThetaCov[1] ) ;


	for ( int i = 0 ; i < blockCnt ; ++i )
	{
		struct _block &e = regions.exonBlocks[i] ;
		double avgDepth = regions.GetAvgDepth( e ) ;
		printf( "%s %" PRId64 " %" PRId64 " %d %d %c %c %lf %lf %lf %lf %lf ", alignments.GetChromName( e.chrId ), e.start + 1, e.end + 1, e.leftType, e.rightType, 
			e.leftStrand, e.rightStrand, avgDepth, 
			e.leftRatio, e.rightRatio, leftClassifier[i], rightClassifier[i] ) ;
		int prevCnt = e.prevCnt ;
		if ( i > 0 && e.start == regions.exonBlocks[i - 1].end + 1 )
			//&& e.leftType == regions.exonBlocks[i - 1].rightType )
		{
			printf( "%d ", prevCnt + 1 ) ;
			for ( j = 0 ; j < prevCnt ; ++j )
				printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
			printf( "%" PRId64 " ", regions.exonBlocks[i - 1].end + 1 ) ;
		}
		else
		{
			printf( "%d ", prevCnt ) ;
			for ( j = 0 ; j < prevCnt ; ++j )
				printf( "%" PRId64 " ", regions.exonBlocks[ e.prev[j] ].end + 1 ) ;
		}

		int nextCnt = e.nextCnt ;
		if ( i < blockCnt - 1 && e.end == regions.exonBlocks[i + 1].start - 1 ) 
			//&& e.rightType == regions.exonBlocks[i + 1].leftType )
		{
			printf( "%d %" PRId64 " ", nextCnt + 1, regions.exonBlocks[i + 1].start + 1 ) ;
		}
		else
			printf( "%d ", nextCnt ) ;
		for ( j = 0 ; j < nextCnt ; ++j )
			printf( "%" PRId64 " ", regions.exonBlocks[ e.next[j] ].start + 1 ) ;
		printf( "\n" ) ;
	}

	delete[] cov ;
	delete[] covRatio ;
	delete[] leftClassifier ;
	delete[] rightClassifier ;
}
