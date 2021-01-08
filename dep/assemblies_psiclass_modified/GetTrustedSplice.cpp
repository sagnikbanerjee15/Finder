#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include "alignments.hpp"

#define MAX(x, y) (((x)<(y))?(y):(x))
#define MIN(x, y) (((x)<(y))?(x):(y))

char usage[] = "Usage: ./trust-splice splice_file_list one_bam_file [OPTIONS]\n"
		"Options:\n"
		"\t-a FLOAT: average number of supported reads from the samples (default: 0.5)\n" ;

struct _intron
{
	int chrId ;
	int start, end ;
	double support ;
	int uniqSupport ;
	int secSupport ;
	char strand ;
	int editDist ;

	int sampleSupport ;
} ;

struct _site
{
	int chrId ;
	int pos ;
	double support ;

	char strand ;
	int associatedIntronCnt ; 
} ;

bool CompIntrons( struct _intron a, struct _intron b )
{
	if ( a.chrId != b.chrId )
		return a.chrId < b.chrId ;
	else if ( a.start != b.start )
		return a.start < b.start ;
	else if ( a.end != b.end )
		return a.end < b.end ;
	return false ;
}

bool CompSites( struct _site a, struct _site b )
{
	if ( a.chrId != b.chrId )
		return a.chrId < b.chrId ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	return false ;
}

void CoalesceIntrons( std::vector<struct _intron> &introns )
{
	std::sort( introns.begin(), introns.end(), CompIntrons ) ;
	int size = introns.size() ;
	int i, k ;
	k = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( introns[i].chrId == introns[k].chrId && introns[i].start == introns[k].start 
			&& introns[i].end == introns[k].end )
		{
			introns[k].support += introns[i].support ;
			introns[k].uniqSupport += introns[i].uniqSupport ;
			introns[k].secSupport += introns[i].secSupport ;
			introns[k].sampleSupport += introns[i].sampleSupport ;

			if ( introns[k].strand == '?' )
				introns[k].strand = introns[i].strand ;	
		}
		else
		{
			++k ;
			introns[k] = introns[i] ;
		}
	}
	introns.resize( k + 1 ) ;
}

void CoalesceSites( std::vector<struct _site> &sites )
{
	std::sort( sites.begin(), sites.end(), CompSites ) ;
	int size = sites.size() ;
	int i, k ;
	k = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( sites[i].chrId == sites[k].chrId && sites[i].pos == sites[k].pos ) 
		{
			sites[k].support += sites[i].support ;
			sites[k].associatedIntronCnt += sites[i].associatedIntronCnt ;
		}
		else
		{
			++k ;
			sites[k] = sites[i] ;
		}
	}
	sites.resize( k + 1 ) ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	FILE *fpList ;	
	FILE *fp ;
	char spliceFile[2048] ;
	char chrName[1024], strand[3] ;
	int start, end, uniqSupport, secSupport, uniqEditDistance, secEditDistance ;
	double support ;
	Alignments alignments ;
	std::vector<struct _intron> introns ;
	std::vector<struct _site> sites ;
	double averageSupportThreshold = 0.5 ;

	if ( argc <= 1 )
	{
		printf( "%s", usage ) ;
		exit( 1 ) ;
	}

	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[ i ], "-a" ) )
		{
			averageSupportThreshold = atof( argv[i + 1] ) ;
			++i ;
		}
		else
		{
			printf( "Unknown option: %s", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	alignments.Open( argv[2] ) ;

	// Get the number of samples.
	fpList = fopen( argv[1], "r" ) ;
	int sampleCnt = 0 ;
	while ( fscanf( fpList, "%s", spliceFile ) != EOF )
	{
		++sampleCnt ;
	}
	fclose( fpList ) ;

	// Get all the introns.
	fpList = fopen( argv[1], "r" ) ;
	while ( fscanf( fpList, "%s", spliceFile ) != EOF )
	{
		fp = fopen( spliceFile, "r" ) ;	
		while (  fscanf( fp, "%s %d %d %lf %s %d %d %d %d", chrName, &start, &end, &support,
				strand, &uniqSupport, &secSupport, &uniqEditDistance, &secEditDistance ) != EOF )
		{
			
			if ( support <= 0 )
				support = 0.1 ;
			else if ( support == 1 && sampleCnt > 5 )
				support = 0.75 ;

			struct _intron ni ;
			ni.chrId = alignments.GetChromIdFromName( chrName ) ;
			ni.start = start ;
			ni.end = end ;
			ni.support = support ;
			ni.strand = strand[0] ;
			ni.uniqSupport = uniqSupport ;
			ni.secSupport = secSupport ;
			ni.sampleSupport = 1 ;
			ni.editDist = uniqEditDistance + secEditDistance ;
			introns.push_back( ni ) ;
		}
		fclose( fp ) ;

		CoalesceIntrons( introns ) ;
	}
	fclose( fpList ) ;

	// Obtain the split sites.
	int intronCnt = introns.size() ;
	for ( i = 0 ; i < intronCnt ; ++i )
	{
		struct _site ns ;
		ns.chrId = introns[i].chrId ;
		ns.associatedIntronCnt = 1 ;
		ns.support = introns[i].support ;
		ns.strand = introns[i].strand ;

		ns.pos = introns[i].start ;
		sites.push_back( ns ) ;
		ns.pos = introns[i].end ;
		sites.push_back( ns ) ;
	}
	CoalesceSites( sites ) ;

	// Get the chromsomes with too many split sites.
	int siteCnt = sites.size() ;
	std::vector<bool> badChrom ;

	badChrom.resize( alignments.GetChromCount() ) ;
	int size = alignments.GetChromCount() ;
	for ( i = 0 ; i < size ; ++i )
		badChrom[i] = false ;

	for ( i = 0 ; i < siteCnt ; ) 	
	{
		for ( j = i + 1 ; sites[j].chrId == sites[i].chrId ; ++j )	
			;
		//printf( "%s %d %d:\n", alignments.GetChromName( sites[i].chrId ), alignments.GetChromLength( sites[i].chrId ), j - i ) ;
		if ( ( j - i ) * 20 > alignments.GetChromLength( sites[i].chrId ) )
			badChrom[ sites[i].chrId ] = true ;
		i = j ;
	}
	
	// Output the valid introns.
	k = 0 ;
	double unit = sampleCnt / 50 ;
	if ( unit < 1 )
		unit = 1 ;
	
	int longIntronSize ; 
	std::vector<int> intronSizes ;
	for (i = 0 ; i < intronCnt ; ++i)
		intronSizes.push_back( introns[i].end - introns[i].start + 1 ) ;
	std::sort( intronSizes.begin(), intronSizes.end() ) ;
	longIntronSize = intronSizes[ int(intronCnt * 0.99) ] ;
	if ( longIntronSize > 100000 )
		longIntronSize = 100000 ;
	for ( i = 0 ; i < intronCnt ; ++i )
	{
		if ( introns[i].support / sampleCnt < averageSupportThreshold )
			continue ;

		if ( badChrom[ introns[i].chrId ] )
		{
			if ( introns[i].sampleSupport <= sampleCnt / 2 )	
				continue ;
		}
		
		if ( introns[i].uniqSupport <= 2 && ( introns[i].support / sampleCnt < 1 || introns[i].sampleSupport < MIN( 2, sampleCnt ) ) )
			continue ;
		
		//Locate the two split sites.
		while ( sites[k].chrId < introns[i].chrId || ( sites[k].chrId == introns[i].chrId && sites[k].pos < introns[i].start ) )
		{
			++k ;
		}
		int a, b ;
		a = k ;
		for ( b = a + 1 ; b < siteCnt ; ++b )
		{
			if ( sites[b].chrId == introns[i].chrId && sites[b].pos == introns[i].end )
				break ;
		}
		double siteSupport = MAX( sites[a].support, sites[b].support ) ;
		//if ( introns[i].start == 100233371 && introns[i].end == 100236850 )
		//	printf( "test: %lf %lf\n", introns[i].support, siteSupport) ;
		if ( introns[i].support < siteSupport / 10.0 )
		{
			double needSample = MIN( ( -log( introns[i].support  / siteSupport ) / log( 10.0 ) + 1 ) * unit, sampleCnt ) ;
			if ( introns[i].sampleSupport < needSample )
				continue ;
		}
		
		if ( sampleCnt >= 100 ) //&& introns[i].end - introns[i].start + 1 >= 50000 )
		{
			if ( introns[i].sampleSupport <= sampleCnt * 0.01 )
				continue ;
		}
		/*if ( sampleCnt >= 50 )
		{
			// just some randomly intron.
			if ( introns[i].sampleSupport == 1 
				&& ( sites[a].associatedIntronCnt == 1 || sites[b].associatedIntronCnt == 1 ) )
				continue ;
		}*/

		/*if (  introns[i].end - introns[i].start + 1 < 50 )
		{
			int needSample = MIN( ( 5 - ( introns[i].end - introns[i].start + 1 ) / 10 ) * unit, sampleCnt ) ;
			int flag = 0 ;

			if ( introns[i].sampleSupport < needSample )
				flag = 1 ;
			if ( flag == 1 )
				continue ;
		}*/
		if ( introns[i].strand == '?' && sampleCnt == 1 && 
			( introns[i].support < 5 || introns[i].uniqSupport == 0 || introns[i].support < 2 * introns[i].editDist ) )
			continue ;
		
		// Since the strand is uncertain, the alinger may make different decision sample from sample. 
		//	To keep this intron, one of its splice sites must be more supported than adjacent splice sites.
		if ( introns[i].strand == '?' && introns[i].sampleSupport <= 0.5 * sampleCnt ) 
		{
			int s, e ;
			int l ;
			int cnt = 0 ;
			for ( l = 0 ; l < 2 ; ++l )
			{
				int ind = ( l == 0 ) ? a : b ;
				double max = sites[ind].support ;
				int maxTag = ind ;
				for ( s = ind - 1 ; s >= 0 && sites[s].chrId == sites[ind].chrId ; --s )
				{
					if ( sites[s].pos + 7 < sites[s + 1].pos )
						break ;
					if ( sites[s].support >= max )
					{
						max = sites[s].support ;
						maxTag = s ;
					}
				}

				for ( e = ind + 1 ; e < siteCnt && sites[e].chrId == sites[ind].chrId ; ++e )
				{
					if ( sites[e].pos - 7 > sites[e - 1].pos )
						break ;
					if ( sites[e].support >= max )
					{
						max = sites[e].support ;
						maxTag = e ;
					}
				}

				if ( maxTag == ind )
					++cnt ;
			}
			if ( cnt == 0 )
				continue ;
		}

		// Test whether this a intron coming from a wrong strand
		/*if ( b - a + 1 >= 10 && introns[i].strand != '?' && introns[i].sampleSupport <= 0.5 * sampleCnt )
		{
			int plusStrand = 0 ;
			int minusStrand = 0 ;
			int l ;
			int s, e ;

			if ( introns[i].strand == '+' )
				plusStrand = 2 ;
			else
				minusStrand = 2 ;

			for ( l = 0 ; l < 2 ; ++l )
			{
				int ind = ( l == 0 ) ? a : b ;
				for ( s = ind - 1 ; s >= 0 && sites[s].chrId == sites[ind].chrId ; --s )
				{
					if ( sites[s].pos + 10000 < sites[s + 1].pos )
						break ;
					
					if ( sites[s].strand == '+' )
						++plusStrand ; 
					else if ( sites[s].strand == '-' ) 
						++minusStrand ; 
				}

				for ( e = ind + 1 ; e < siteCnt && sites[e].chrId == sites[ind].chrId ; ++e )
				{
					if ( sites[e].pos - 10000 > sites[e - 1].pos )
						break ;
					
					if ( sites[e].strand == '+' )
						++plusStrand ; 
					else if ( sites[e].strand == '-' ) 
						++minusStrand ;
				}
			}

			if ( introns[i].start == 161517978 )
				printf( "capture: %d %d %d %d\n", a, b, minusStrand, plusStrand) ;
			
			if ( introns[i].strand == '+' && minusStrand >= 20  && plusStrand == 2 )
				continue ;
			else if ( introns[i].strand == '-' && plusStrand >= 20 && minusStrand == 2 )
				continue ;

		}*/
	
		// Filter a intron if one of its splice site associated with too many introns.
		/*if ( sites[a].associatedIntronCnt >= 10 )
		{
			int needSample = MIN( sites[a].associatedIntronCnt / 100 * sampleCnt, sampleCnt ) ;
			if ( introns[i].support < sites[a].support / 100 && 
				introns[i].sampleSupport < needSample )
			{
				continue ;		
			}
		}
		if ( sites[b].associatedIntronCnt >= 10 )
		{
			int needSample = MIN( sites[b].associatedIntronCnt / 100 * sampleCnt, sampleCnt ) ;
			if ( introns[i].support < sites[b].support / 100 && 
				introns[i].sampleSupport < needSample )
			{
				continue ;
			}
		}*/


		// Test for long intron
		if ( introns[i].end - introns[i].start + 1 >= longIntronSize )
		{
			int needSample = MIN( ( ( introns[i].end - introns[i].start + 1 ) / longIntronSize + 1 ) * unit, sampleCnt ) ;
			int flag = 0 ;
			if ( (double)introns[i].uniqSupport / ( introns[i].uniqSupport + introns[i].secSupport ) < 0.1 
				|| introns[i].uniqSupport / sampleCnt < 1 
				|| introns[i].sampleSupport < needSample )
				flag = 1 ;
			if ( flag == 1 && introns[i].end - introns[i].start + 1 >= 3 * longIntronSize )
				continue ;
			if ( flag == 1 && ( sites[a].associatedIntronCnt > 1 || sites[b].associatedIntronCnt > 1 || introns[i].sampleSupport <= 1 ) ) // an intron may connect two genes
				continue ;
		}
		
				
		printf( "%s %d %d 10 %c 10 0 0 0\n", alignments.GetChromName( introns[i].chrId ), introns[i].start, introns[i].end, introns[i].strand ) ;
	}

	return 0 ;
}
