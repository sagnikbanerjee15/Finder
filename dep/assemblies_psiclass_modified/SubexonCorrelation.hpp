#ifndef _MOURISL_CLASSES_SUBEXONCORRELATION_HEADER
#define _MOURISL_CLASSES_SUBEXONCORRELATION_HEADER

#include "SubexonGraph.hpp"

#include <stdio.h>
#include <math.h>

class SubexonCorrelation
{
private:
	struct _subexon *lastSubexons ;
	std::vector<FILE *> fileList ;
	int offset ;
	int prevSeCnt ;
	int seCnt ;

	double **correlation ;
		
	int OverlapSize( int s0, int e0, int s1, int e1 )
	{
		int s = -1, e = -1 ;
		if ( e0 < s1 || s0 > e1 )
			return 0 ;
		s = s0 > s1 ? s0 : s1 ;
		e = e0 < e1 ? e0 : e1 ;
		return e - s + 1 ;
	}
public:
	SubexonCorrelation() 
	{
		offset = 0 ;
		lastSubexons = NULL ;
		correlation = NULL ;
		prevSeCnt = 0 ;
	}

	~SubexonCorrelation()
	{
		int cnt = fileList.size() ;
		int i ;
		for ( i = 0 ; i < cnt ; ++i )
			fclose( fileList[i] ) ;
		delete[] lastSubexons ;
	}

	void Initialize( char *f )
	{
		FILE *fpSl ;
		fpSl = fopen( f, "r" ) ;
		char buffer[1024] ;
		while ( fgets( buffer, sizeof( buffer ), fpSl ) != NULL )
		{
			int len = strlen( buffer ) ;
			if ( buffer[len - 1] == '\n' )
			{
				buffer[len - 1] = '\0' ;
				--len ;

			}
			FILE *fp = fopen( buffer, "r" ) ;
			fileList.push_back( fp ) ;
		}

		lastSubexons = new struct _subexon[ fileList.size() ] ;
		int i, cnt ;
		cnt = fileList.size() ;
		for ( i = 0 ; i < cnt ; ++i )
			lastSubexons[i].chrId = -1 ;
	
	}
	
	// We assume the subexons only contains the subexons we are interested in.
	// And we assume that each time we call this function, the subexons we are interested are 
	// sorted in order.
	void ComputeCorrelation( struct _subexon *subexons, int cnt, Alignments &alignments )
	{
		int sampleCnt = fileList.size() ;
		if ( sampleCnt <= 1 )
			return ;
		int i, j, k ;
		seCnt = cnt ;

		// Obtain the depth matrix.
		double **depth ;
		depth = new double* [sampleCnt] ; 
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			depth[i] = new double[ cnt ] ; 
			memset( depth[i], 0, sizeof( double ) * cnt ) ;
		}
		
		char buffer[2048] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			struct _subexon se ;
			bool useLastSubexon = false ;
			if ( lastSubexons[i].chrId != -1 )
			{
				se = lastSubexons[i] ;
				useLastSubexon = true ;
			}
			else
				se.chrId = -1 ;

			// Locate the subexon
			while ( 1 )
			{
				if ( se.chrId < subexons[0].chrId 
					|| ( se.chrId == subexons[0].chrId && se.end < subexons[0].start ) )
				{
					if ( fgets( buffer, sizeof( buffer ), fileList[i] ) == NULL )
					{
						buffer[0] = '\0' ;
						break ;
					}
					if ( buffer[0] == '#' )
						continue ;
					
					SubexonGraph::InputSubexon( buffer, alignments, se ) ;
					--se.start ; --se.end ;
					useLastSubexon = false ;
					continue ;
				}
				break ;
			}
			if ( buffer[0] == '\0' )
			{
				lastSubexons[i] = se ;	
				continue ;
			}

			// Process the subexon.
			int tag = 0 ;
			while ( 1 )
			{
				lastSubexons[i] = se ;
				
				if ( useLastSubexon == false )
				{
					SubexonGraph::InputSubexon( buffer, alignments, se ) ; 
					--se.start ; --se.end ;
				}

				if ( se.chrId > subexons[cnt - 1].chrId || se.start > subexons[cnt - 1].end )
					break ;

				while ( 1 )	
				{
					if ( tag > cnt || subexons[tag].end >= se.start )	
						break ;
					++tag ;
				}

				for ( j = tag ; j < cnt && subexons[j].start <= se.end ; ++j )
				{
					int overlap = OverlapSize( se.start, se.end, subexons[j].start, subexons[j].end ) ; 
					depth[i][j] += overlap * se.avgDepth ;		
				}
				
				if ( fgets( buffer, sizeof( buffer ), fileList[i] ) == NULL )
					break ;
			}
		}
		// Normalize the depth
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			for ( j = 0 ; j < cnt ; ++j )
				depth[i][j] /= ( subexons[j].end - subexons[j].start + 1 ) ;
		}

		// Compute the correlation.
		double *avg = new double[cnt] ;
		double *var = new double[cnt] ;
		if ( correlation != NULL )
		{
			for ( i = 0 ; i < prevSeCnt ; ++i )
				delete[] correlation[i] ;
			delete[] correlation ;
		}

		correlation = new double*[cnt] ;
		for ( i = 0 ; i < cnt ; ++i )
			correlation[i] = new double[cnt] ;
		
		memset( avg, 0, sizeof( double ) * cnt ) ;
		memset( var, 0, sizeof( double ) * cnt ) ;
		for ( i = 0 ; i < cnt ; ++i ) 
		{
			//avg[i] = 0 ;
			//var[i] = 0 ;
			memset( correlation[i], 0, sizeof( double ) * cnt ) ;
			correlation[i][i] = 1 ;
		}

		for ( i = 0 ; i < sampleCnt ; ++i )
			for ( j = 0 ; j < cnt ; ++j )
			{
				avg[j] += depth[i][j] ;
				var[j] += depth[i][j] * depth[i][j] ;
			}
		for ( i = 0 ; i < cnt ; ++i )
		{
			avg[i] /= sampleCnt ;
			var[i] = var[i] / sampleCnt - avg[i] * avg[i] ;
		}
		
		for ( i = 0 ; i < sampleCnt ; ++i )
			for ( j = 0 ; j < cnt ; ++j )
				for ( k = j + 1 ; k < cnt ; ++k )
				{
					//if ( subexons[j].geneId == subexons[k].geneId )
					correlation[j][k] += depth[i][j] * depth[i][k] ;	
				}
		for ( j = 0 ; j < cnt ; ++j )
			for ( k = j + 1 ; k < cnt ; ++k )
			{
				if ( var[j] > 1e-6 && var[k] > 1e-6 )
				{
					correlation[j][k] = ( correlation[j][k] / sampleCnt - avg[j] * avg[k] ) / sqrt( var[j] * var[k] ) ;
				}
				else
					correlation[j][k] = 0 ;
				//printf( "%d %d %d\n", j, k, cnt ) ;
				correlation[k][j] = correlation[j][k] ;
			}
		//printf( "%lf\n", correlation[0][1] ) ;
		// Release the memory.
		for ( i = 0 ; i < sampleCnt ; ++i )
			delete[] depth[i] ;
		delete[] depth ;
		delete[] avg ;
		delete[] var ;

		prevSeCnt = cnt ;
	}

	double Query( int i, int j )
	{
		if ( fileList.size() > 1 )
			return correlation[i][j] ;
		else
			return 0 ;
	}

	void Assign( const SubexonCorrelation &c )
	{
		int i ;
		fileList = c.fileList ;
		if ( fileList.size() <= 1 )
			return ;

		if ( correlation != NULL )
		{
			for ( i = 0 ; i < seCnt ; ++i )
				delete[] correlation[i] ;
			delete[] correlation ;
		}
		
		seCnt = c.seCnt ;
		correlation = new double*[seCnt] ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			correlation[i] = new double[seCnt] ;
			memcpy( correlation[i], c.correlation[i], sizeof( double ) * seCnt ) ;
		}
	}
} ;

#endif
