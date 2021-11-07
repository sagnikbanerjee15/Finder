#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include "alignments.hpp"
#include "blocks.hpp"
#include "stats.hpp"
#include "SubexonGraph.hpp"

char usage[] = "combineSubexons [options]\n"
	       "Required options:\n"
	       "\t-s STRING: the path to the predicted subexon information. Can use multiple -s to specify multiple subexon prediction files\n" 
	       "\t\tor\n"
	       "\t--ls STRING: the path to file of the list of the predicted subexon information.\n" 
	       "Optional options:\n"
	       "\t-q FLOAT: the quantile of samples to determine the extension of subexon soft boundaries. (default: 0.5)"
	       ;

struct _overhang
{
	int cnt ; //  the number of samples support this subexon.
	int validCnt ; // The number of samples that are used for compute probability.
	int length ;
	double classifier ;
} ;

struct _intronicInfo
{
	int chrId ;
	int start, end ;
	int leftSubexonIdx, rightSubexonIdx ;
	double irClassifier ;
	int irCnt ;
	int validIrCnt ;
	struct _overhang leftOverhang, rightOverhang ; // leftOverhangClassifier is for the overhang subexon at the left side of this intron.
} ;

struct _seInterval
{
	int chrId ;
	int start, end ;
	int type ; // 0-subexon, 1-intronicInfo
	int idx ;
} ;

struct _subexonSplit
{
	int chrId ;
	int pos ;
	int type ; //1-start of a subexon. 2-end of a subexon 
	int splitType ; //0-soft boundary, 1-start of an exon, 2-end of an exon.
	int strand ;

	int weight ;
} ;

struct _interval // exon or intron
{
	int chrId ;
	int start, end ;
	int strand ;
	int sampleSupport ;
} ;

struct _subexonSupplement // supplement the subexon structure defined in SubexonGraph.
{
	int *nextSupport ;		
	int *prevSupport ;
} ;

char buffer[4096] ;

bool CompSubexonSplit( struct _subexonSplit a, struct _subexonSplit b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.pos != b.pos )
		return a.pos < b.pos ;
	else if ( a.type != b.type )
	{
		// the split site with no strand information should come first.
		/*if ( a.splitType != b.splitType )
		{
			if ( a.splitType == 0 ) 
				return true ;
			else if ( b.splitType == 0 )
				return false ;
		}*/
		return a.type < b.type ;
	}
	else if ( a.splitType != b.splitType )
	{
		//return a.splitType < b.splitType ;
		if ( a.splitType == 0 )
			return true ;
		else if ( b.splitType == 0 )
			return false ;

		if ( a.type == 1 )
			return a.splitType > b.splitType ;
		else
			return a.splitType < b.splitType ;
	}
	
	return false ;
}

bool CompInterval( struct _interval a, struct _interval b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.start != b.start )
		return a.start < b.start ;
	else if ( a.end != b.end )
		return a.end < b.end ;
	return false ;
}

bool CompSeInterval( struct _seInterval a, struct _seInterval b )
{
	if ( a.chrId < b.chrId )
		return true ;
	else if ( a.chrId > b.chrId )
		return false ;
	else if ( a.start < b.start )
		return true ;
	else if ( a.start > b.start )
		return false ;
	else if ( a.end < b.end )
		return true ;
	else
		return false ;
}

// Keep this the same as in SubexonInfo.cpp.
double TransformCov( double c )
{
	double ret ;
	//return sqrt( c ) - 1 ;

	if ( c <= 2 + 1e-6 )
		ret = 1e-6 ;
	else
		ret = c - 2 ;
	
	return ret ;
}

double GetUpdateMixtureGammaClassifier( double ratio, double cov, double piRatio, double kRatio[2], double thetaRatio[2],
	double piCov, double kCov[2], double thetaCov[2], bool conservative )
{
	double p1 = 0, p2 ;

	cov = TransformCov( cov ) ;
	if ( cov < ( kCov[0] - 1 ) * thetaCov[0] )
		cov = ( kCov[0] - 1 ) * thetaCov[0] ;

	if ( ratio > 0 )
		p1 = MixtureGammaAssignment( ratio, piRatio, kRatio, thetaRatio ) ;
	// Make sure cov > 1?	
	p2 = MixtureGammaAssignment( cov, piCov, kCov, thetaCov ) ;
	double ret = 0 ;

	if ( conservative )
	{
		if ( p1 >= p2 ) // we should use ratio.
			ret = LogGammaDensity( ratio, kRatio[1], thetaRatio[1] ) 
				- LogGammaDensity( ratio, kRatio[0], thetaRatio[0] ) ;
		else
			ret = LogGammaDensity( cov, kCov[1], thetaCov[1] ) 	
				- LogGammaDensity( cov, kCov[0], thetaCov[0] ) ;
	}
	else
	{
		if ( p1 >= p2 ) // we should use ratio.
			ret = LogGammaDensity( ratio, kRatio[1], thetaRatio[1] ) 
				- LogGammaDensity( ratio, kRatio[0], thetaRatio[0] ) ;
		else
			ret = LogGammaDensity( cov, kCov[1], thetaCov[1] ) 	
				- LogGammaDensity( cov, kCov[0], thetaCov[0] ) ;
	}
	return ret ;
}

double GetPValueOfGeneEnd( double cov )
{
	if ( cov <= 2.0 )
		return 1.0 ;
	double tmp = 2.0 * ( sqrt( cov ) - log( cov ) ) ;
	if ( tmp <= 0 )
		return 1.0 ;
	return 2.0 * alnorm( tmp, true ) ;
}

char StrandNumToSymbol( int strand )
{
	if ( strand > 0 )
		return '+' ;
	else if ( strand < 0 )
		return '-' ;
	else
		return '.' ;
}

int StrandSymbolToNum( char c )
{
	if ( c == '+' )
		return 1 ;
	else if ( c == '-' )
		return -1 ;
	else
		return 0 ;
}

int *MergePositions( int *old, int ocnt, int *add, int acnt, int &newCnt ) //, int **support )
{
	int i, j, k ;
	int *ret ;
	if ( acnt == 0 )	
	{
		newCnt = ocnt ;
		return old ;
	}
	if ( ocnt == 0 )
	{
		newCnt = acnt ;
		ret = new int[acnt] ;
		//*support = new int[acnt] ;
		for ( i = 0 ; i < acnt ; ++i )
		{
			//(*support)[i] = 1 ;
			ret[i] = add[i] ;		
		}
		return ret ;
	}
	newCnt = 0 ;
	for ( i = 0, j = 0 ; i < ocnt && j < acnt ; )
	{
		if ( old[i] < add[j] )
		{
			++i ;
			++newCnt ;
		}
		else if ( old[i] == add[j] )
		{
			++i ; ++j ;
			++newCnt ;
		}
		else 
		{
			++j ;
			++newCnt ;
		}
	}
	newCnt = newCnt + ( ocnt - i ) + ( acnt - j ) ;
	// no new elements.
	if ( newCnt == ocnt )
	{
		/*i = 0 ;
		for ( j = 0 ; j < acnt ; ++j )
		{
			for ( ; old[i] < add[j] ; ++i )
				;
			++(*support)[i] ;
		}*/
		return old ;
	}
	k = 0 ;
	//delete []old ;
	ret = new int[ newCnt ] ;
	//int *bufferSupport = new int[newCnt] ;
	for ( i = 0, j = 0 ; i < ocnt && j < acnt ; )
	{
		if ( old[i] < add[j] )
		{
			ret[k] = old[i] ;
			//bufferSupport[k] = (*support)[i] ;
			++i ;
			++k ;
		}
		else if ( old[i] == add[j] )
		{
			ret[k] = old[i] ;
			//bufferSupport[k] = (*support)[i] + 1 ;
			++i ; ++j ;
			++k ;
		}
		else 
		{
			ret[k] = add[j] ;
			//bufferSupport[k] = 1 ;
			++j ;
			++k ;
		}
	}
	for ( ; i < ocnt ; ++i, ++k )
	{
		ret[k] = old[i] ;
		//bufferSupport[k] = (*support)[i] ;
	}
	for ( ; j < acnt ; ++j, ++k )
	{
		ret[k] = add[j] ;
		//bufferSupport[k] = 1 ;
	}
	delete[] old ;
	//delete[] *support ;
	//*support = bufferSupport ;
	return ret ;
}

void CoalesceSubexonSplits( std::vector<struct _subexonSplit> &splits, int mid )  
{
	int i, j, k ;
	int cnt = splits.size() ;
	//std::sort( splits.begin(), splits.end(), CompSubexonSplit ) ;
	
	std::vector<struct _subexonSplit> sortedSplits ;
	sortedSplits.resize( cnt ) ;
	
	k = 0 ;
	for ( i = 0, j = mid ; i < mid && j < cnt ; ++k )
	{
		if ( CompSubexonSplit( splits[i], splits[j] ) )		
		{
			sortedSplits[k] = splits[i] ;
			++i ;
		}
		else
		{
			sortedSplits[k] = splits[j] ;
			++j ;
		}
	}
	for ( ; i < mid ; ++i, ++k )
		sortedSplits[k] = splits[i] ;
	for ( ; j < cnt ; ++j, ++k )
		sortedSplits[k] = splits[j] ;
	splits = sortedSplits ;	

	k = 0 ;
	for ( i = 1 ; i < cnt ; ++i )
	{
		if ( splits[i].chrId == splits[k].chrId && splits[i].pos == splits[k].pos && splits[i].type == splits[k].type && splits[i].splitType == splits[k].splitType 
			&& splits[i].strand == splits[k].strand )	
		{
			splits[k].weight += splits[i].weight ;
		}
		else
		{
			++k ;
			splits[k] = splits[i] ;
		}
	}
	splits.resize( k + 1 ) ;
}

void CoalesceDifferentStrandSubexonSplits( std::vector<struct _subexonSplit> &splits )
{
	int i, j, k, l ;
	int cnt = splits.size() ;
	k = 0 ;
	for ( i = 0 ; i < cnt ;  )
	{
		for ( j = i + 1 ; j < cnt ; ++j )
		{
			if ( splits[i].chrId == splits[j].chrId && splits[i].pos == splits[j].pos && splits[i].type == splits[j].type && splits[i].splitType == splits[j].splitType )
				continue ;
			break ;
		}

		int maxWeight = -1 ;
		int weightSum = 0 ;
		int strand = splits[i].strand ;
		for ( l = i ; l < j ; ++l )
		{
			weightSum += splits[l].weight ;
			if ( splits[l].strand != 0 && splits[l].weight > maxWeight )
			{
				strand = splits[l].strand ;
				maxWeight = splits[l].weight ;
			}
		}
		//printf( "%d: %d %d %d\n", splits[i].pos, i, j, strand ) ;
		splits[k] = splits[i] ;
		splits[k].strand = strand ;
		splits[k].weight = weightSum ;
		++k ;

		i = j ;
	}
	splits.resize( k ) ;
}


void CoalesceIntervals( std::vector<struct _interval> &intervals )
{
	int i, k ;
	std::sort( intervals.begin(), intervals.end(), CompInterval ) ;
	int cnt = intervals.size() ;
	k = 0 ;
	for ( i = 1 ; i < cnt ; ++i )
	{
		if ( intervals[i].chrId == intervals[k].chrId && intervals[i].start == intervals[k].start && intervals[i].end == intervals[k].end )
			intervals[k].sampleSupport += intervals[i].sampleSupport ;
		else
		{
			++k ;
			intervals[k] = intervals[i] ;
		}
	}
	intervals.resize( k + 1 ) ;
}

void CleanIntervalIrOverhang( std::vector<struct _interval> &irOverhang )
{
	int i, j, k ;
	std::sort( irOverhang.begin(), irOverhang.end(), CompInterval ) ;

	int cnt = irOverhang.size() ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( irOverhang[i].start == -1 )
			continue ;


		// locate the longest interval start at the same coordinate.
		int tag = i ;
	
		for ( j = i + 1 ; j < cnt ; ++j )	
		{
			if ( irOverhang[j].chrId != irOverhang[i].chrId || irOverhang[j].start != irOverhang[i].start )
				break ;
			if ( irOverhang[j].start == -1 )
				continue ;
			tag = j ;
		}
		
		for ( k = i ; k < tag ; ++k )
		{
			irOverhang[k].start = -1 ;
		}
		
		for ( k = tag + 1 ; k < cnt ; ++k )
		{
			if ( irOverhang[k].chrId != irOverhang[tag].chrId || irOverhang[k].start > irOverhang[tag].end )
				break ;
			if ( irOverhang[k].end <= irOverhang[tag].end )
			{
				irOverhang[k].start = -1 ;
			}
		}
	}
	
	k = 0 ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		if ( irOverhang[i].start == -1 )
			continue ;
		irOverhang[k] = irOverhang[i] ;
		++k ;
	}
	irOverhang.resize( k ) ;
}

// Remove the connection that does not match the boundary
//  of subexons.
void CleanUpSubexonConnections( std::vector<struct _subexon> &subexons )
{
	int seCnt = subexons.size() ;
	int i, j, k, m ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].prevCnt > 0 )	
		{
			for ( k = i ; k >= 0 ; --k )
				if ( subexons[k].chrId != subexons[i].chrId || subexons[k].end <= subexons[i].prev[0] )
					break ;
			if ( subexons[k].chrId != subexons[i].chrId )
				++k ;
			m = 0 ;
			for ( j = 0 ; j < subexons[i].prevCnt ; ++j )
			{
				for ( ; k <= i ; ++k )
					if ( subexons[k].end >= subexons[i].prev[j] )
						break ;
				if ( subexons[k].end == subexons[i].prev[j] 
					&& ( SubexonGraph::IsSameStrand( subexons[k].rightStrand, subexons[i].leftStrand ) || subexons[k].end + 1 == subexons[i].start ) )
				{
					subexons[i].prev[m] = subexons[i].prev[j] ;
					++m ;
				}
			}
			subexons[i].prevCnt = m ;
		}

		m = 0 ;
		k = i ;
		for ( j = 0 ; j < subexons[i].nextCnt ; ++j )
		{
			for ( ; k < seCnt ; ++k )
				if ( subexons[k].chrId != subexons[i].chrId || subexons[k].start >= subexons[i].next[j] )			
					break ;
				if ( subexons[k].start == subexons[i].next[j] 
					&& ( SubexonGraph::IsSameStrand( subexons[i].rightStrand, subexons[k].leftStrand ) || subexons[i].end + 1 == subexons[k].start ) )
				{
					subexons[i].next[m] = subexons[i].next[j] ;
					++m ;
				}
		}
		subexons[i].nextCnt = m ;
	}
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	FILE *fp ;
	std::vector<char *> files ;

	Blocks regions ;
	Alignments alignments ;

	double exonSoftBoundaryMergeQuantile = 0.5 ;

	if ( argc == 1 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-s" ) )
		{
			files.push_back( argv[i + 1] ) ;
			++i ;
			continue ;
		}
		else if ( !strcmp( argv[i], "--ls" ) )
		{
			FILE *fpLs = fopen( argv[i + 1], "r" ) ;
			char buffer[1024] ;
			while ( fgets( buffer, sizeof( buffer ), fpLs ) != NULL )
			{
				int len = strlen( buffer ) ;
				if ( buffer[len - 1] == '\n' )
				{
					buffer[len - 1] = '\0' ;
					--len ;

				}
				char *fileName = strdup( buffer ) ;
				files.push_back( fileName ) ;
			}
		}
		else if ( !strcmp( argv[i], "-q" ) )
		{
			sscanf( argv[i + 1], "%lf", &exonSoftBoundaryMergeQuantile ) ;
			++i ;
		}
	}
	int fileCnt = files.size() ;
	// Obtain the chromosome ids through bam file.
	fp = fopen( files[0], "r" ) ;		
	if ( fgets( buffer, sizeof( buffer ), fp ) != NULL )
	{
		int len = strlen( buffer ) ;
		buffer[len - 1] = '\0' ;
		alignments.Open( buffer + 1 ) ;
	}
	fclose( fp ) ;

	// Collect the split sites of subexons.
	std::vector<struct _subexonSplit> subexonSplits ;
	std::vector<struct _interval> intervalIrOverhang ; // intervals contains ir and overhang.
	std::vector<struct _interval> introns ;
	std::vector<struct _interval> exons ;


	for ( k = 0 ; k < fileCnt ; ++k )
	{
		fp = fopen( files[k], "r" ) ;		
		struct _subexon se ;
		struct _subexonSplit sp ;
		char chrName[50] ;
		int origSize = subexonSplits.size() ;
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
				continue ;

			SubexonGraph::InputSubexon( buffer, alignments, se, false) ;
			// Record all the intron rentention, overhang from the samples
			if ( ( se.leftType == 2 && se.rightType == 1 ) 
				|| ( se.leftType == 2 && se.rightType == 0 )
				|| ( se.leftType == 0 && se.rightType == 1 ) )  
			{
				struct _interval si ;
				si.chrId = se.chrId ;
				si.start = se.start ;
				si.end = se.end ;

				intervalIrOverhang.push_back( si ) ;
			}

			// Ignore overhang subexons and ir subexons for now.
			if ( ( se.leftType == 0 && se.rightType == 1 ) 
				|| ( se.leftType == 2 && se.rightType == 0 ) 
				||  ( se.leftType == 2 && se.rightType == 1 ) )
				continue ;

			if ( se.leftType == 0 && se.rightType == 0 && ( se.leftClassifier == -1 || se.leftClassifier == 1 ) ) // ignore noisy single-exon island
				continue ;
			if ( se.leftType == 0 && se.rightType == 0 && ( fileCnt >= 10 && se.leftClassifier > 0.99 ) )  
				continue ;

			if ( se.leftType == 1 && se.rightType == 2 ) // a full exon, we allow mixtured strand here.
			{
				struct _interval ni ;
				ni.chrId = se.chrId ;
				ni.start = se.start ;
				ni.end = se.end ;
				ni.strand = se.rightStrand ;  
				ni.sampleSupport = 1 ;
				exons.push_back( ni ) ;
			}


			/*for ( i = 0 ; i < se.nextCnt ; ++i )
			{
				struct _interval ni ;
				ni.chrId = se.chrId ;
				ni.start = se.end ;
				ni.end = se.next[i] ;
				ni.strand = se.rightStrand ;  
				ni.sampleSupport = 1 ;
				if ( ni.start + 1 < ni.end )
					introns.push_back( ni ) ;
			}*/

			sp.chrId = se.chrId ;
			sp.pos = se.start ;
			sp.type = 1 ;
			sp.splitType = se.leftType ;			
			sp.strand = se.leftStrand ;
			sp.weight = 1 ;
			subexonSplits.push_back( sp ) ;

			sp.chrId = se.chrId ;
			sp.pos = se.end ;
			sp.type = 2 ;
			sp.splitType = se.rightType ;				
			sp.strand = se.rightStrand ;
			sp.weight = 1 ;
			subexonSplits.push_back( sp ) ;
			
			/*if ( se.prevCnt > 0 )
				delete[] se.prev ;
			if ( se.nextCnt > 0 )
				delete[] se.next ;*/
		}
		CoalesceIntervals( exons ) ;
		CoalesceIntervals( introns ) ;
		CoalesceSubexonSplits( subexonSplits, origSize ) ;
		CleanIntervalIrOverhang( intervalIrOverhang ) ;
		
		fclose( fp ) ;
	}

	CoalesceDifferentStrandSubexonSplits( subexonSplits ) ;
	
	// Obtain the split sites from the introns.
	int intronCnt = introns.size() ;
	std::vector<struct _subexonSplit> intronSplits ;
	for ( i = 0 ; i < intronCnt ; ++i )
	{
		/*if ( introns[i].sampleSupport < 0.05 * fileCnt )
		{
			continue ;
		}*/
		struct _interval &it = introns[i] ;
		struct _subexonSplit sp ;
		sp.chrId = it.chrId ;
		sp.pos = it.start ;
		sp.type = 2 ;
		sp.splitType = 2 ;			
		sp.strand = it.strand ;
		intronSplits.push_back( sp ) ;

		sp.chrId = it.chrId ;
		sp.pos = it.end ;
		sp.type = 1 ;
		sp.splitType = 1 ;				
		sp.strand = it.strand ;
		intronSplits.push_back( sp ) ;
	}

	// Pair up the split sites to get subexons
	std::sort( intronSplits.begin(), intronSplits.end(), CompSubexonSplit ) ;
	//std::sort( subexonSplits.begin(), subexonSplits.end(), CompSubexonSplit ) ;
	
	// Convert the hard boundary to soft boundary if the split sites is filtered from the introns
	// Seems NO need to do this now.
	int splitCnt = subexonSplits.size() ;
	int intronSplitCnt = intronSplits.size() ;
	k = 0 ;
	//for ( i = 0 ; i < splitCnt ; ++i )
	while ( 0 )
	{
		if ( subexonSplits[i].type != subexonSplits[i].splitType )
			continue ;
			
		while ( k < intronSplitCnt && ( intronSplits[k].chrId < subexonSplits[i].chrId 
			|| ( intronSplits[k].chrId == subexonSplits[i].chrId && intronSplits[k].pos < subexonSplits[i].pos ) ) )			
			++k ;
		j = k ;
		while ( j < intronSplitCnt && intronSplits[j].chrId == subexonSplits[i].chrId 
			&& intronSplits[j].pos == subexonSplits[i].pos && intronSplits[j].splitType != subexonSplits[i].splitType )			
			++j ;
		
		// the split site is filtered.
		if ( j >= intronSplitCnt || intronSplits[j].chrId != subexonSplits[i].chrId ||
			intronSplits[j].pos > subexonSplits[i].pos )
		{
			//printf( "%d %d. %d %d\n", subexonSplits[i].pos, intronSplits[j].pos, intronSplits[j].chrId , subexonSplits[i].chrId ) ;
			subexonSplits[i].splitType = 0 ;
			
			// Convert the adjacent subexon split.
			for ( int l = i + 1 ; i < splitCnt && subexonSplits[l].chrId == subexonSplits[i].chrId 
				&& subexonSplits[l].pos == subexonSplits[i].pos + 1 ; ++l )
			{
				if ( subexonSplits[l].type != subexonSplits[i].type 
					&& subexonSplits[l].splitType == subexonSplits[i].splitType )
				{
					subexonSplits[l].splitType = 0 ;	
				}
			}

			// And the other direction
			for ( int l = i - 1 ; l >= 0 && subexonSplits[l].chrId == subexonSplits[i].chrId 
				&& subexonSplits[l].pos == subexonSplits[i].pos - 1 ; --l )
			{
				if ( subexonSplits[l].type != subexonSplits[i].type 
					&& subexonSplits[l].splitType == subexonSplits[i].splitType )
				{
					subexonSplits[l].splitType = 0 ;	
				}
			}
		}
	}
	intronSplits.clear() ;
	std::vector<struct _subexonSplit>().swap( intronSplits ) ;
	
	// Force the soft boundary that collides with hard boundaries to be hard boundary.
	for ( i = 0 ; i < splitCnt ; ++i )
	{
		if ( subexonSplits[i].splitType != 0 )
			continue ;
		int newSplitType = 0 ;
		int newStrand = subexonSplits[i].strand ;
		for ( j = i + 1 ; j < splitCnt ; ++j )
		{
			if ( subexonSplits[i].type != subexonSplits[j].type || subexonSplits[i].pos != subexonSplits[j].pos ||
					subexonSplits[i].chrId != subexonSplits[j].chrId )
				break ;
			if ( subexonSplits[j].splitType != 0 )
			{
				newSplitType = subexonSplits[j].splitType ;
				newStrand = subexonSplits[j].strand ;
				break ;
			}
		}

		if ( newSplitType == 0 )
		{
			for ( j = i - 1 ; j >= 0 ; --j )
			{
				if ( subexonSplits[i].type != subexonSplits[j].type || subexonSplits[i].pos != subexonSplits[j].pos ||
						subexonSplits[i].chrId != subexonSplits[j].chrId )
					break ;
				if ( subexonSplits[j].splitType != 0 )
				{
					newSplitType = subexonSplits[j].splitType ;
					newStrand = subexonSplits[j].strand ;
					break ;
				}
			}

		}
		/*if ( subexonSplits[i].pos == 154464157 )
		{
			printf( "force conversion: %d %d %d. %d %d\n", subexonSplits[i].pos, subexonSplits[i].splitType, subexonSplits[i].weight, subexonSplits[i + 1].pos, subexonSplits[i + 1].splitType ) ;
		}*/
		subexonSplits[i].splitType = newSplitType ;
		subexonSplits[i].strand = newStrand ;
	}

	/*for ( i = 0 ; i < splitCnt ; ++i )
	{
		printf( "%d: type=%d splitType=%d weight=%d\n", subexonSplits[i].pos, subexonSplits[i].type, subexonSplits[i].splitType, subexonSplits[i].weight ) ;
	}*/
	
	// Build subexons from the collected split sites.
	
	std::vector<struct _subexon> subexons ;
	int diffCnt = 0 ; // |start of subexon split| - |end of subexon split|
	int seCnt = 0 ;
	for ( i = 0 ; i < splitCnt - 1 ; ++i )	
	{
		struct _subexon se ;
		/*if ( subexonSplits[i + 1].pos == 144177260 )
		{
			printf( "%d %d %d: %d %d %d. %d\n", subexonSplits[i].pos, subexonSplits[i].type, subexonSplits[i].splitType, 
				subexonSplits[i + 1].pos, subexonSplits[i + 1].type, subexonSplits[i + 1].splitType, diffCnt ) ;
		}*/

		if ( subexonSplits[i].type == 1 )
			diffCnt += subexonSplits[i].weight ;
		else
			diffCnt -= subexonSplits[i].weight ;

		if ( subexonSplits[i + 1].chrId != subexonSplits[i].chrId )
		{
			diffCnt = 0 ;		
			continue ;
		}

		if ( diffCnt == 0 ) // the interval between subexon
			continue ;

		se.chrId = subexonSplits[i].chrId ;
		se.start = subexonSplits[i].pos ;
		se.leftType = subexonSplits[i].splitType ;
		se.leftStrand = subexonSplits[i].strand ;
		if ( subexonSplits[i].type == 2 )	
		{
			se.leftStrand = 0 ;
			++se.start ;
		}

		se.end = subexonSplits[i + 1].pos ;
		se.rightType = subexonSplits[i + 1].splitType ;
		se.rightStrand = subexonSplits[i + 1].strand ;
		if ( subexonSplits[i + 1].type == 1 )
		{
			se.rightStrand = 0 ;
			--se.end ;
		}
			
		/*if ( se.end == 24613649 )
		{
			//printf( "%d %d %d: %d %d %d. %d\n", subexonSplits[i].pos, subexonSplits[i].type, subexonSplits[i].splitType, 
			//	subexonSplits[i + 1].pos, subexonSplits[i + 1].type, subexonSplits[i + 1].splitType, diffCnt ) ;
			printf( "%d %d %d. %d %d %d\n", se.start, se.leftType, se.leftStrand, se.end, se.rightType, se.rightStrand ) ;
		}*/

		if ( se.start > se.end ) //Note: this handles the case of repeated subexon split.
		{
			// handle the case in sample 0: [...[..]
			// in sample 1:                 [..]...]
			if ( seCnt > 0 && se.end == subexons[seCnt - 1].end && subexons[seCnt - 1].rightType < se.rightType ) 
			{
				subexons[seCnt - 1].rightType = se.rightType ;
				subexons[seCnt - 1].rightStrand = se.rightStrand ;
			}
			continue ;
		}
		se.leftClassifier = se.rightClassifier = 0 ;
		se.lcCnt = se.rcCnt = 0 ;
		
		/*if ( 1 ) //se.chrId == 25 )	
		{
			printf( "%d: %d-%d: %d %d\n", se.chrId, se.start, se.end, se.leftType, se.rightType ) ;
		}*/


		se.next = se.prev = NULL ;
		se.nextCnt = se.prevCnt = 0 ;
		subexons.push_back( se ) ;
		++seCnt ;
	}
	subexonSplits.clear() ;
	std::vector<struct _subexonSplit>().swap( subexonSplits ) ;
	
	// Adjust the split type.
	seCnt = subexons.size() ;
	for ( i = 1 ; i < seCnt ; ++i )
	{
		if ( subexons[i - 1].end + 1 == subexons[i].start )
		{
			if ( subexons[i - 1].rightType == 0 )
				subexons[i - 1].rightType = subexons[i].leftType ;
			if ( subexons[i].leftType == 0 )
				subexons[i].leftType = subexons[i - 1].rightType ;
		}
	}

	// Merge the adjacent soft boundaries 
	std::vector<struct _subexon> rawSubexons = subexons ;
	int exonCnt = exons.size() ;
	subexons.clear() ;

	k = 0 ; // hold index for exon.
	for ( i = 0 ; i < seCnt ;  )
	{
		/*if ( rawSubexons[k].rightType == 0 && rawSubexons[i].leftType == 0 
			&& rawSubexons[k].end + 1 == rawSubexons[i].start )			
		{
			rawSubexons[k].end = rawSubexons[i].end ;
			rawSubexons[k].rightType = rawSubexons[i].rightType ;
			rawSubexons[k].rightStrand = rawSubexons[i].rightStrand ;
		}
		else
		{
			subexons.push_back( rawSubexons[k] ) ;		
			k = i ;
		}*/

		while ( k < exonCnt && ( exons[k].chrId < rawSubexons[i].chrId 
				|| ( exons[k].chrId == rawSubexons[i].chrId && exons[k].start < rawSubexons[i].start ) ) )
			++k ;

		for ( j = i + 1 ; j < seCnt ; ++j )
		{
			if ( rawSubexons[j - 1].chrId != rawSubexons[j].chrId || rawSubexons[j - 1].rightType != 0 || rawSubexons[j].leftType != 0 
				|| ( fileCnt > 1 && rawSubexons[j - 1].end + 1 != rawSubexons[j].start )
				|| ( fileCnt == 1 && rawSubexons[j - 1].end + 50 < rawSubexons[j].start ) )
				break ;
		}
		// rawsubexons[i...j-1] will be merged.

		/*if ( rawSubexons[i].start == 119625875 )
		{
			printf( "merge j-1: %d %d %d %d\n", rawSubexons[j - 1].end, rawSubexons[j - 1].rightType,
				rawSubexons[j].start, rawSubexons[j].leftType ) ;
		}*/
		bool merge = true ;
		if ( rawSubexons[i].leftType == 1 && rawSubexons[j - 1].rightType == 2 && j - i > 1 
			&& rawSubexons[j - 1].end - rawSubexons[i].start >= 1000 )
		{
			merge = false ;
			int sampleSupport = 0 ;
			for ( int l = k ; l < exonCnt ; ++l )
			{
				if ( exons[l].chrId != rawSubexons[i].chrId || exons[l].start > rawSubexons[i].start )
					break ;
				if ( exons[l].end == rawSubexons[j - 1].end )
				{
					merge = true ;
					sampleSupport = exons[l].sampleSupport ;
					break ;
				}
			}

			if ( merge == true && rawSubexons[j - 1].end - rawSubexons[i].start >= 1000 )
			{
				if ( sampleSupport <= 0.2 * fileCnt )
				{
					merge = false ;
				}
			}
			
			if ( merge == false )
			{
				if ( j - i >= 3 )
				{
					double adjustQuantile = exonSoftBoundaryMergeQuantile ;
					if ( adjustQuantile > 0.5 ) 
						adjustQuantile = 0.5 ;
					rawSubexons[i].end = rawSubexons[ i + ( j - 1 - i ) * adjustQuantile ].start ;
					rawSubexons[j - 1].start = rawSubexons[ i + ( j - 1 - i ) * (1 - adjustQuantile) ].end ;
				}

				if ( rawSubexons[i].end + 1 == rawSubexons[j - 1].start )
				{
					--rawSubexons[i].end ;
					++rawSubexons[j - 1].start ;
				}
				subexons.push_back( rawSubexons[i] ) ;
				subexons.push_back( rawSubexons[j - 1] ) ;
			}
		}

		if ( merge )
		{
			rawSubexons[i].end = rawSubexons[j - 1].end ;
			rawSubexons[i].rightType = rawSubexons[j - 1].rightType ;
			rawSubexons[i].rightStrand = rawSubexons[j - 1].rightStrand ;
			
			if ( rawSubexons[i].leftType == 0 && rawSubexons[i].rightType != 0 )
			{
				rawSubexons[i].start = rawSubexons[ i + ( j - 1 - i ) * ( 1 - exonSoftBoundaryMergeQuantile ) ].start ;
			}
			else if ( rawSubexons[i].rightType == 0 && rawSubexons[i].leftType != 0 )
			{
				rawSubexons[i].end = rawSubexons[ i + ( j - 1 - i ) * exonSoftBoundaryMergeQuantile ].end ;
			}

			subexons.push_back( rawSubexons[i] ) ;		
		}

		i = j ;
	}
	exons.clear() ;
	std::vector<struct _interval>().swap( exons ) ;

	// Remove overhang, ir subexons intron created after putting multiple sample to gether.
	// eg: s0: [......)
	//     s1: [...]--------[....]
	//     s2: [...]..)-----[....]
	// Though the overhang from s2 is filtered in readin, there will a new overhang created combining s0,s1.
	// 	But be careful about how to compute the classifier for the overhang part contributed from s0.
	// Furthermore, note that the case of single-exon island showed up in intron retention region after combining is not possible when get here.
	//    eg: s0:[...]-----[...]
	//        s1:      (.)
	//        s2:[.............]
	//  After merge adjacent soft boundaries, the single-exon island will disappear.
	rawSubexons = subexons ;
	seCnt = subexons.size() ;
	subexons.clear() ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( ( rawSubexons[i].leftType == 2 && rawSubexons[i].rightType == 1 )		// ir
			|| ( rawSubexons[i].leftType == 2 && rawSubexons[i].rightType == 0 )    // overhang	
			|| ( rawSubexons[i].leftType == 0 && rawSubexons[i].rightType == 1 ) )  
			continue ;
		subexons.push_back( rawSubexons[i] ) ;
	}
	
	// Remove the single-exon island if it overlaps with intron retentioned or overhang.
	rawSubexons = subexons ;
	seCnt = subexons.size() ;
	subexons.clear() ;
	k = 0 ;
	std::sort( intervalIrOverhang.begin(), intervalIrOverhang.end(), CompInterval ) ;
	int irOverhangCnt = intervalIrOverhang.size() ;

	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( rawSubexons[i].leftType != 0 || rawSubexons[i].rightType != 0 )
		{
			subexons.push_back( rawSubexons[i] ) ;
			continue ;
		}
		
		while ( k < irOverhangCnt )
		{
			// Locate the interval that before the island
			if ( intervalIrOverhang[k].chrId < rawSubexons[i].chrId 
				|| ( intervalIrOverhang[k].chrId == rawSubexons[i].chrId && intervalIrOverhang[k].end < rawSubexons[i].start ) )
			{
				++k ;
				continue ;
			}
			break ;
		}
		bool overlap = false ;
		for ( j = k ; j < irOverhangCnt ; ++j )
		{
			if ( intervalIrOverhang[j].chrId > rawSubexons[i].chrId || intervalIrOverhang[j].start > rawSubexons[i].end )
				break ;
			if ( ( intervalIrOverhang[j].start <= rawSubexons[i].start && intervalIrOverhang[j].end >= rawSubexons[i].start ) 
				|| ( intervalIrOverhang[j].start <= rawSubexons[i].end && intervalIrOverhang[j].end >= rawSubexons[i].end ) )
			{
				overlap = true ;
				break ;
			}
		}

		if ( !overlap )
			subexons.push_back( rawSubexons[i] ) ;
	}
	rawSubexons.clear() ;
	std::vector<struct _subexon>().swap( rawSubexons ) ;

	intervalIrOverhang.clear() ;
	std::vector<struct _interval>().swap( intervalIrOverhang ) ;
	
	// Create the dummy intervals.
	seCnt = subexons.size() ;
	std::vector<struct _intronicInfo> intronicInfos ;
	std::vector<struct _seInterval> seIntervals ;
	std::vector<struct _subexonSupplement> subexonInfo ;
	
	//subexonInfo.resize( seCnt ) ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		struct _seInterval ni ; // new interval
		ni.start = subexons[i].start ;
		ni.end = subexons[i].end ;
		ni.type = 0 ;
		ni.idx = i ;
		ni.chrId = subexons[i].chrId ;
		seIntervals.push_back( ni ) ;

		/*if ( subexons[i].end == 42671717 )	
		{
			printf( "%d: %d-%d: %d\n", subexons[i].chrId, subexons[i].start, subexons[i].end, subexons[i].rightType ) ;
		}*/
		//subexonInfo[i].prevSupport = subexonInfo[i].nextSupport = NULL ;
		
		/*int nexti ;
		for ( nexti = i + 1 ; nexti < seCnt ; ++nexti )
			if ( subexons[ nexti ].leftType == 0 && subexons[nexti].rightType == 0 )*/

		if ( i < seCnt - 1 && subexons[i].chrId == subexons[i + 1].chrId && 
			subexons[i].end + 1 < subexons[i + 1].start &&
			subexons[i].rightType + subexons[i + 1].leftType != 0 )
		{
			// Only consider the intervals like ]..[,]...(, )...[
			// The case like ]...] is actaully things like ][...] in subexon perspective,
			// so they won't pass the if-statement
			struct _intronicInfo nii ; // new intronic info
			ni.start = subexons[i].end + 1 ;
			ni.end = subexons[i + 1].start - 1 ;
			ni.type = 1 ;
			ni.idx = intronicInfos.size() ;
			seIntervals.push_back( ni ) ;
			
			nii.chrId = subexons[i].chrId ;
			nii.start = ni.start ;
			nii.end = ni.end ; 
			nii.leftSubexonIdx = i ;
			nii.rightSubexonIdx = i + 1 ;
			nii.irClassifier = 0 ;
			nii.irCnt = 0 ;
			nii.validIrCnt = 0 ;
			nii.leftOverhang.cnt = 0 ;
			nii.leftOverhang.validCnt = 0 ;
			nii.leftOverhang.length = 0 ;
			nii.leftOverhang.classifier = 0 ;
			nii.rightOverhang.cnt = 0 ;
			nii.rightOverhang.validCnt = 0 ;
			nii.rightOverhang.length = 0 ;
			nii.rightOverhang.classifier = 0 ;
			intronicInfos.push_back( nii ) ;
			/*if ( nii.end == 23667 )
			{
				printf( "%d %d. %d (%d %d %d)\n", nii.start, nii.end, subexons[i].rightType, subexons[i+1].start, subexons[i + 1].end, subexons[i + 1].leftType ) ;
			}*/
		}
	}
	
	// Go through all the files to get some statistics number
	double avgIrPiRatio = 0 ;
	double avgIrPiCov = 0 ;
	double irPiRatio, irKRatio[2], irThetaRatio[2] ; // Some statistical results
	double irPiCov, irKCov[2], irThetaCov[2] ;
	
	double avgOverhangPiRatio = 0 ;
	double avgOverhangPiCov = 0 ;
	double overhangPiRatio, overhangKRatio[2], overhangThetaRatio[2] ; // Some statistical results
	double overhangPiCov, overhangKCov[2], overhangThetaCov[2] ;
	
	for ( k = 0 ; k < fileCnt ; ++k )
	{
		fp = fopen( files[k], "r" ) ;		
		
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
			{
				char buffer2[100] ;
				sscanf( buffer, "%s", buffer2 ) ;	
				if ( !strcmp( buffer2, "#fitted_ir_parameter_ratio:" ) )
				{
					// TODO: ignore certain samples if the coverage seems wrong.
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiRatio, buffer2, &irKRatio[0], buffer2, &irThetaRatio[0],
								buffer2, &irKRatio[1], buffer2, &irThetaRatio[1] ) ;	
					avgIrPiRatio += irPiRatio ;
				}
				else if ( !strcmp( buffer2, "#fitted_ir_parameter_cov:" ) )
				{
				}
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiRatio, buffer2, &overhangKRatio[0], buffer2, &overhangThetaRatio[0],
								buffer2, &overhangKRatio[1], buffer2, &overhangThetaRatio[1] ) ;	
					avgOverhangPiRatio += overhangPiRatio ;
				}
			}
			else
				break ;
		}
		fclose( fp ) ;
	}
	avgIrPiRatio /= fileCnt ;
	avgOverhangPiRatio /= fileCnt ;

	// Go through all the files to put statistical results into each subexon.
	std::vector< struct _subexon > sampleSubexons ;
	int subexonCnt = subexons.size() ;
	for ( k = 0 ; k < fileCnt ; ++k )
	{
		//if ( k == 220 )
		//	exit( 1 ) ;
		fp = fopen( files[k], "r" ) ;		
		struct _subexon se ;
		struct _subexonSplit sp ;
		char chrName[50] ;
		
		sampleSubexons.clear() ;

		int tag = 0 ;
		while ( fgets( buffer, sizeof( buffer), fp  ) != NULL )
		{
			if ( buffer[0] == '#' )
			{
				char buffer2[200] ;
				sscanf( buffer, "%s", buffer2 ) ;	
				if ( !strcmp( buffer2, "#fitted_ir_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiRatio, buffer2, &irKRatio[0], buffer2, &irThetaRatio[0],
								buffer2, &irKRatio[1], buffer2, &irThetaRatio[1] ) ;	
				}
				else if ( !strcmp( buffer2, "#fitted_ir_parameter_cov:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &irPiCov, buffer2, &irKCov[0], buffer2, &irThetaCov[0],
								buffer2, &irKCov[1], buffer2, &irThetaCov[1] ) ;	
				}
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_ratio:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiRatio, buffer2, &overhangKRatio[0], buffer2, &overhangThetaRatio[0],
								buffer2, &overhangKRatio[1], buffer2, &overhangThetaRatio[1] ) ;	
				}	
				else if ( !strcmp( buffer2, "#fitted_overhang_parameter_cov:" ) )
				{
					sscanf( buffer, "%s %s %lf %s %lf %s %lf %s %lf %s %lf", 
								buffer2, buffer2, &overhangPiCov, buffer2, &overhangKCov[0], buffer2, &overhangThetaCov[0],
								buffer2, &overhangKCov[1], buffer2, &overhangThetaCov[1] ) ;	
				}
				continue ;
			}
			else 
				break ;

			//SubexonGraph::InputSubexon( buffer, alignments, se, true ) ;
			//sampleSubexons.push_back( se ) ;
		}
		
		//int sampleSubexonCnt = sampleSubexons.size() ;
		int intervalCnt = seIntervals.size() ;
		//for ( i = 0 ; i < sampleSubexonCnt ; ++i )	
		int iterCnt = 0 ;
		while ( 1 )
		{
			if ( iterCnt > 0 && fgets( buffer, sizeof( buffer), fp  ) == NULL)
				break ;
			++iterCnt ;

			struct _subexon se ;
			SubexonGraph::InputSubexon( buffer, alignments, se, true ) ;

			while ( tag < intervalCnt )	
			{
				if ( seIntervals[tag].chrId < se.chrId || 
					( seIntervals[tag].chrId == se.chrId && seIntervals[tag].end < se.start ) )
				{
					++tag ;
					continue ;
				}
				else
					break ;
			}
			
			for ( j = tag ; j < intervalCnt ; ++j )
			{
				if ( seIntervals[j].start > se.end || seIntervals[j].chrId > se.chrId ) // terminate if no overlap.
					break ;
				int idx ;	
				
				if ( seIntervals[j].type == 0 )
				{
					idx = seIntervals[j].idx ;
					if ( subexons[idx].leftType == 1 && se.leftType == 1 && subexons[idx].start == se.start )
					{
						double tmp = se.leftClassifier ;
						if ( se.leftClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].leftClassifier -= 2.0 * log( tmp ) ;		
						++subexons[idx].lcCnt ;
						subexons[idx].prev = MergePositions( subexons[idx].prev, subexons[idx].prevCnt, 
										se.prev, se.prevCnt, subexons[idx].prevCnt ) ;

						if ( se.rightType == 0 ) // a gene end here
						{
							for ( int l = idx ; l < subexonCnt ; ++l )
							{
								if ( l > idx && ( subexons[l].end > subexons[l - 1].start + 1 
									|| subexons[l].chrId != subexons[l - 1].chrId ) )				
									break ;
								if ( subexons[l].rightType == 2 )
								{
									double adjustAvgDepth = se.avgDepth ;
									if ( se.end - se.start + 1 >= 100 )
										adjustAvgDepth += se.avgDepth * 100.0 / ( se.end - se.start + 1 ) ;
									else
										adjustAvgDepth *= 2 ;
									double p = GetPValueOfGeneEnd( adjustAvgDepth ) ;
									//if ( se.end - se.start + 1 >= 500 && p > 0.001 )
									//	p = 0.001 ;
									
									subexons[l].rightClassifier -= 2.0 * log( p ) ; 			
									++subexons[l].rcCnt ;
									break ;
								}
							}
						}
					}
					//if ( se.chrId == 25 )	
					//	printf( "(%d %d %d. %d) (%d %d %d. %d)\n", se.chrId, se.start, se.end, se.rightType, subexons[idx].chrId, subexons[idx].start, subexons[idx].end, subexons[idx].rightType ) ;
					if ( subexons[idx].rightType == 2 && se.rightType == 2 && subexons[idx].end == se.end )
					{
						double tmp = se.rightClassifier ;
						if ( se.rightClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].rightClassifier -= 2.0 * log( tmp ) ;
						++subexons[idx].rcCnt ;

						subexons[idx].next = MergePositions( subexons[idx].next, subexons[idx].nextCnt, 
										se.next, se.nextCnt, subexons[idx].nextCnt ) ;

						if ( se.leftType == 0 )
						{
							for ( int l = idx ; l >= 0 ; --l )
							{
								if ( l < idx && ( subexons[l].end < subexons[l + 1].start - 1 
											|| subexons[l].chrId != subexons[l + 1].chrId ) )				
									break ;
								if ( subexons[l].leftType == 1 )
								{
									double adjustAvgDepth = se.avgDepth ;
									if ( se.end - se.start + 1 >= 100 )
										adjustAvgDepth += se.avgDepth * 100.0 / ( se.end - se.start + 1 ) ;
									else
										adjustAvgDepth *= 2 ;
									double p = GetPValueOfGeneEnd( adjustAvgDepth ) ;
									//if ( se.end - se.start + 1 >= 500 && p >= 0.001 )
									//	p = 0.001 ;
									subexons[l].leftClassifier -= 2.0 * log( p ) ; 			
									++subexons[l].lcCnt ;
									break ;
								}
							}
						}
					}

					if ( subexons[idx].leftType == 0 && subexons[idx].rightType == 0
						&& se.leftType == 0 && se.rightType == 0 ) // the single-exon island.
					{
						double tmp = se.leftClassifier ;
						if ( se.leftClassifier == 0 )
							tmp = 1e-7 ;
						subexons[idx].leftClassifier -= 2.0 * log( tmp ) ;
						subexons[idx].rightClassifier = subexons[idx].leftClassifier ;
						++subexons[idx].lcCnt ;
						++subexons[idx].rcCnt ;
					}
				}
				else if ( seIntervals[j].type == 1 )
				{
					idx = seIntervals[j].idx ;
					// Overlap on the left part of intron
					if ( se.start <= intronicInfos[idx].start && se.end < intronicInfos[idx].end 
						&& subexons[ intronicInfos[idx].leftSubexonIdx ].rightType != 0 )
					{
						int len = se.end - intronicInfos[idx].start + 1 ;
						intronicInfos[idx].leftOverhang.length += len ;
						++intronicInfos[idx].leftOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at right hand side, 
						// otherwise, this part is not an intron and won't show up in intronic Info.
						if ( se.leftType == 2 )
						{
							if ( se.leftRatio > 0 && se.avgDepth > 1 )
							{
								++intronicInfos[idx].leftOverhang.validCnt ;

								double update = GetUpdateMixtureGammaClassifier( se.leftRatio, se.avgDepth, 
										overhangPiRatio, overhangKRatio, overhangThetaRatio, 
										overhangPiCov, overhangKCov, overhangThetaCov, false ) ;
								intronicInfos[idx].leftOverhang.classifier += update ;				
							}
						}
						else if ( se.leftType == 1 )
						{
							++intronicInfos[idx].leftOverhang.validCnt ;
							double update = GetUpdateMixtureGammaClassifier( 1.0, se.avgDepth, 
									overhangPiRatio, overhangKRatio, overhangThetaRatio, 
									overhangPiCov, overhangKCov, overhangThetaCov, true ) ;
							intronicInfos[idx].leftOverhang.classifier += update ;			
							
							int seIdx = intronicInfos[idx].leftSubexonIdx ;
							subexons[seIdx].rightClassifier -= 2.0 * log( GetPValueOfGeneEnd( se.avgDepth ) ) ;
							++subexons[ seIdx ].rcCnt ; 
						}
						// ignore the contribution of single-exon island here?
					}
					// Overlap on the right part of intron
					else if ( se.start > intronicInfos[idx].start && se.end >= intronicInfos[idx].end 
							&& subexons[ intronicInfos[idx].rightSubexonIdx ].leftType != 0 )
					{
						int len = intronicInfos[idx].end - se.start + 1 ;
						intronicInfos[idx].rightOverhang.length += len ;
						++intronicInfos[idx].rightOverhang.cnt ;
						
						// Note that the sample subexon must have a soft boundary at left hand side, 
						// otherwise, this won't show up in intronic Info
						if ( se.rightType == 1 )
						{
							if ( se.rightRatio > 0 && se.avgDepth > 1 )
							{
								++intronicInfos[idx].rightOverhang.validCnt ;

								double update = GetUpdateMixtureGammaClassifier( se.rightRatio, se.avgDepth, 
										overhangPiRatio, overhangKRatio, overhangThetaRatio, 
										overhangPiCov, overhangKCov, overhangThetaCov, false ) ;
								intronicInfos[idx].rightOverhang.classifier += update ;				
							}
						}
						else if ( se.rightType == 2 )
						{
							++intronicInfos[idx].rightOverhang.validCnt ;

							double update = GetUpdateMixtureGammaClassifier( 1, se.avgDepth, 
									overhangPiRatio, overhangKRatio, overhangThetaRatio, 
									overhangPiCov, overhangKCov, overhangThetaCov, true ) ;
							intronicInfos[idx].rightOverhang.classifier += update ;				

							int seIdx = intronicInfos[idx].rightSubexonIdx ;
							/*if ( subexons[ seIdx ].start == 6873648 )
							{
								printf( "%lf %lf: %lf %lf %lf\n", subexons[seIdx].leftClassifier, GetPValueOfGeneEnd( se.avgDepth ), se.avgDepth, sqrt( se.avgDepth ), log( se.avgDepth ) )  ;
							}*/
							subexons[seIdx].leftClassifier -= 2.0 * log( GetPValueOfGeneEnd( se.avgDepth ) ) ;
							++subexons[ seIdx ].lcCnt ;
						}
					}
					// Intron is fully contained in this sample subexon, then it is a ir candidate
					else if ( se.start <= intronicInfos[idx].start && se.end >= intronicInfos[idx].end )
					{
						if ( se.leftType == 2 && se.rightType == 1 )		
						{
							double ratio = regions.PickLeftAndRightRatio( se.leftRatio, se.rightRatio ) ;
							++intronicInfos[idx].irCnt ;
							if ( ratio > 0 && se.avgDepth > 1 )
							{
								double update = GetUpdateMixtureGammaClassifier( ratio, se.avgDepth,
										irPiRatio, irKRatio, irThetaRatio,
										irPiCov, irKCov, irThetaCov, true ) ;
								//if ( intronicInfos[idx].start == 37617368 )
								//	printf( "hi %lf %d %d: %d %d\n", update, se.start, se.end, intronicInfos[idx].start, intronicInfos[idx].end ) ;
								intronicInfos[idx].irClassifier += update ;
								++intronicInfos[idx].validIrCnt ;
							}
						}
						else if ( se.leftType == 1 || se.rightType == 2 )
						{
							//intronicInfos[idx].irClassifier += LogGammaDensity( 4.0, irKRatio[1], irThetaRatio[1] )
							//                                         - LogGammaDensity( 4.0, irKRatio[0], irThetaRatio[0] ) ;
							/*if ( se.start == 37617368 )
							{
								printf( "%lf: %lf %lf\n", se.avgDepth, MixtureGammaAssignment( ( irKCov[0] - 1 ) * irThetaCov[0], irPiRatio, irKCov, irThetaCov ),
									MixtureGammaAssignment( TransformCov( 4.0 ), irPiRatio, irKCov, irThetaCov ) ) ;
							}*/
							if ( se.avgDepth > 1 )
							{
								// let the depth be the threshold to determine.
								double update = GetUpdateMixtureGammaClassifier( 4.0, se.avgDepth,
										irPiRatio, irKRatio, irThetaRatio,
										irPiCov, irKCov, irThetaCov, true ) ;
								//if ( intronicInfos[idx].start == 36266630 )
								//	printf( "hi %lf %d %d: %d %d\n", update, se.start, se.end, intronicInfos[idx].start, intronicInfos[idx].end ) ;
								intronicInfos[idx].irClassifier += update ;
								++intronicInfos[idx].irCnt ;
								++intronicInfos[idx].validIrCnt ;
							}
						}
						else
						{
							// the intron is contained in a overhang subexon from the sample or single-exon island
						}
					}
					// sample subexon is contained in the intron.
					else
					{
						// Do nothing.			
					}
				}
			}

			//if ( se.nextCnt > 0 )
				delete[] se.next ;
			//if ( se.prevCnt > 0 )
				delete[] se.prev ;
		}
		fclose( fp ) ;
		
		/*for ( i = 0 ; i < sampleSubexonCnt ; ++i )
		{
			if ( sampleSubexons[i].nextCnt > 0 )
				delete[] sampleSubexons[i].next ;
			if ( sampleSubexons[i].prevCnt > 0 )
				delete[] sampleSubexons[i].prev ;
		}*/
	}

	CleanUpSubexonConnections( subexons ) ;

	// Convert the temporary statistics number into formal statistics result.
	for ( i = 0 ; i < subexonCnt ; ++i ) 
	{
		struct _subexon &se = subexons[i] ;
		if ( se.leftType == 0 && se.rightType == 0 ) // single-exon txpt.
		{
			se.leftClassifier = se.rightClassifier = 1 - chicdf( se.rightClassifier, 2 * se.rcCnt ) ;
		}
		else
		{
			if ( se.leftType == 1 )
			{
				se.leftClassifier = 1 - chicdf( se.leftClassifier, 2 * se.lcCnt ) ; 	
			}
			else
				se.leftClassifier = -1 ;

			if ( se.rightType == 2 )
				se.rightClassifier = 1 - chicdf( se.rightClassifier, 2 * se.rcCnt ) ;
			else
				se.rightClassifier = -1 ;
		}
	}

	int iiCnt = intronicInfos.size() ; //intronicInfo count
	for ( i = 0 ; i < iiCnt ; ++i )
	{
		struct _intronicInfo &ii = intronicInfos[i] ;
		if ( ii.validIrCnt > 0 )
		{
			for ( j = 0 ; j < fileCnt - ii.validIrCnt ; ++j )
			{
				ii.irClassifier -= log( 10.0 ) ;
			}
			/*if ( ii.validIrCnt < fileCnt * 0.15 )
				ii.irClassifier -= log( 1000.0 ) ;
			else if ( ii.validIrCnt < fileCnt * 0.5 )
				ii.irClassifier -= log( 100.0 ) ;*/
			ii.irClassifier = (double)1.0 / ( 1.0 + exp( ii.irClassifier + log( 1 - avgIrPiRatio ) - log( avgIrPiRatio ) ) ) ;
		}
		else
			ii.irClassifier = -1 ;
		
		if ( ii.leftOverhang.validCnt > 0 )
			ii.leftOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.leftOverhang.classifier + 
						log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
		else
			ii.leftOverhang.classifier = -1 ;

		if ( ii.rightOverhang.validCnt > 0 )
			ii.rightOverhang.classifier = (double)1.0 / ( 1.0 + exp( ii.rightOverhang.classifier + 
						log( 1 - avgOverhangPiRatio ) - log( avgOverhangPiRatio ) ) ) ;
		else
			ii.rightOverhang.classifier = -1 ;
	}

	// Change the classifier for the hard boundaries if its adjacent intron has intron retention classifier
	//    which collide with overhang subexon.
	int intervalCnt = seIntervals.size() ;
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		if ( seIntervals[i].type == 1 && intronicInfos[ seIntervals[i].idx ].irCnt > 0 )
		{
			int idx = seIntervals[i].idx ;
			if ( intronicInfos[idx].leftOverhang.cnt > 0 )
			{
				int k = seIntervals[i - 1].idx ;
				// Should aim for more conservative?
				if ( subexons[k].rightClassifier > intronicInfos[idx].leftOverhang.classifier )
					subexons[k].rightClassifier = intronicInfos[idx].leftOverhang.classifier ;
			}

			if ( intronicInfos[idx].rightOverhang.cnt > 0 )
			{
				int k = seIntervals[i + 1].idx ;
				if ( subexons[k].leftClassifier > intronicInfos[idx].rightOverhang.classifier )
					subexons[k].leftClassifier = intronicInfos[idx].rightOverhang.classifier ;
			}
		}
	}
	
	// Output the result.
	for ( i = 0 ; i < intervalCnt ; ++i )
	{
		if ( seIntervals[i].type == 0 )
		{
			struct _subexon &se = subexons[ seIntervals[i].idx ] ;
			
			char ls, rs ;

			ls = StrandNumToSymbol( se.leftStrand ) ;
			rs = StrandNumToSymbol( se.rightStrand ) ;

			printf( "%s %d %d %d %d %c %c -1 -1 -1 %lf %lf ", alignments.GetChromName( se.chrId ), se.start, se.end,
					se.leftType, se.rightType, ls, rs, se.leftClassifier, se.rightClassifier ) ;
			if ( i > 0 && seIntervals[i - 1].chrId == seIntervals[i].chrId 
				&& seIntervals[i - 1].end + 1 == seIntervals[i].start 
				&& !( seIntervals[i - 1].type == 0 && 
					subexons[ seIntervals[i - 1].idx ].rightType != se.leftType ) 
				&& !( seIntervals[i - 1].type == 1 && intronicInfos[ seIntervals[i - 1].idx ].irCnt == 0
					&& intronicInfos[ seIntervals[i - 1].idx ].rightOverhang.cnt == 0 ) 
				&& ( se.prevCnt == 0 || se.start - 1 != se.prev[ se.prevCnt - 1 ] ) ) // The connection showed up in the subexon file.
			{
				printf( "%d ", se.prevCnt + 1 ) ;
				for ( j = 0 ; j < se.prevCnt ; ++j )
					printf( "%d ", se.prev[j] ) ;
				printf( "%d ", se.start - 1 ) ;
			}
			else
			{
				printf( "%d ", se.prevCnt ) ;
				for ( j = 0 ; j < se.prevCnt ; ++j )
					printf( "%d ", se.prev[j] ) ;
			}

			if ( i < intervalCnt - 1 && seIntervals[i].chrId == seIntervals[i + 1].chrId 
				&& seIntervals[i].end == seIntervals[i + 1].start - 1
				&& !( seIntervals[i + 1].type == 0 &&
					subexons[ seIntervals[i + 1].idx ].leftType != se.rightType ) 
				&& !( seIntervals[i + 1].type == 1 && intronicInfos[ seIntervals[i + 1].idx ].irCnt == 0
					&& intronicInfos[ seIntervals[i + 1].idx ].leftOverhang.cnt == 0 ) 
				&& ( se.nextCnt == 0 || se.end + 1 != se.next[0] ) )
			{
				printf( "%d %d ", se.nextCnt + 1, se.end + 1 ) ;
			}
			else
				printf( "%d ", se.nextCnt ) ;
			for ( j = 0 ; j < se.nextCnt ; ++j )
				printf( "%d ", se.next[j] ) ;
			printf( "\n" ) ;
		}
		else if ( seIntervals[i].type == 1 )
		{
			struct _intronicInfo &ii = intronicInfos[ seIntervals[i].idx ] ;
			if ( ii.irCnt > 0 )
			{
				printf( "%s %d %d 2 1 . . -1 -1 -1 %lf %lf 1 %d 1 %d\n",
					alignments.GetChromName( ii.chrId ), ii.start, ii.end, 
					ii.irClassifier, ii.irClassifier,
					seIntervals[i - 1].end, seIntervals[i + 1].start ) ;
			}
			else
			{
				// left overhang.
				if ( ii.leftOverhang.cnt > 0 )
				{
					printf( "%s %d %d 2 0 . . -1 -1 -1 %lf %lf 1 %d 0\n",
						alignments.GetChromName( ii.chrId ), ii.start, 
						ii.start + ( ii.leftOverhang.length /  ii.leftOverhang.cnt ) - 1,
						ii.leftOverhang.classifier, ii.leftOverhang.classifier,
						ii.start - 1 ) ; 
				}

				// right overhang.
				if ( ii.rightOverhang.cnt > 0 )
				{
					printf( "%s %d %d 0 1 . . -1 -1 -1 %lf %lf 0 1 %d\n",
						alignments.GetChromName( ii.chrId ), 
						ii.end - ( ii.rightOverhang.length / ii.rightOverhang.cnt ) + 1, ii.end,
						ii.rightOverhang.classifier, ii.rightOverhang.classifier,
						ii.end + 1 ) ;
				}

			}
		}
	}

	return 0 ;
}

