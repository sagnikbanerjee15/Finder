#ifndef _MOURISL_CLASSES_SUBEXONGRAPH_HEADER
#define _MOURISL_CLASSES_SUBEXONGRAPH_HEADER

#include "alignments.hpp"
#include "blocks.hpp"

struct _subexon
{
	int chrId ;
	int geneId ;
	int start, end ;
	int leftType, rightType ;
	double avgDepth ;
	//double ratio, classifier ;
	double leftRatio, rightRatio ;
	double leftClassifier, rightClassifier ;
	int lcCnt, rcCnt ;
	int leftStrand, rightStrand ;
	
	int nextCnt, prevCnt ;
	int *next, *prev ;
	
	bool canBeStart, canBeEnd ;
} ;

struct _geneInterval
{
	int startIdx, endIdx ;
	int start, end ; // The start and end of a gene interval might be adjusted, so it does not 
			// need to be match with the corresponding subexons
} ;

class SubexonGraph
{
private:
	int *visit ;
	double classifierThreshold ;

	int usedGeneId ;
	int baseGeneId ;

	// The function to assign gene ids to subexons.
	void SetGeneId( int tag, int strand, struct _subexon *subexons, int seCnt, int id ) ;
	void GetGeneBoundary( int tag, int &boundary, int timeStamp ) ;
	void UpdateGeneId( struct _subexon *subexons, int seCnt ) ;
public:
	std::vector<struct _subexon> subexons ;
	std::vector<struct _geneInterval> geneIntervals ;

	~SubexonGraph() 
	{
		int i ;
		int size = subexons.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( subexons[i].next )
				delete[] subexons[i].next ;
			if ( subexons[i].prev )
				delete[] subexons[i].prev ;
		}
	} 

	SubexonGraph( double classifierThreshold, Alignments &bam, FILE *fpSubexon ) 
	{ 
		// Read in the subexons
		rewind( fpSubexon ) ; 
		char buffer[2048] ;
		int subexonCnt ;
		int i, j, k ;
		while ( fgets( buffer, sizeof( buffer ), fpSubexon ) != NULL )
		{
			if ( buffer[0] == '#' )
				continue ;

			struct _subexon se ;
			InputSubexon( buffer, bam, se, true ) ;

			// filter.
			if ( ( se.leftType == 0 && se.rightType == 0 ) 
				|| ( se.leftType == 0 && se.rightType == 1 ) 	// overhang
				|| ( se.leftType == 2 && se.rightType == 0 ) // overhang
				|| ( se.leftType == 2 && se.rightType == 1 ) ) // ir
			{
				if ( ( se.leftType == 0 && se.rightType == 1 ) 
					|| ( se.leftType == 2 && se.rightType == 0 ) ) // if the overhang is too small
				{
					if ( se.end - se.start + 1 <= 7 )
					{
						if ( se.next )
							delete[] se.next ;
						if ( se.prev )
							delete[] se.prev ;
						continue ;
					}
				}

				if ( se.leftClassifier >= classifierThreshold || se.leftClassifier < 0 )
				{
					if ( se.next )
						delete[] se.next ;
					if ( se.prev )
						delete[] se.prev ;
					continue ;
				}
			}
			
			// Adjust the coordinate.
			subexons.push_back( se ) ;	
		}

		// Convert the coordinate to index
		// Note that each coordinate can only associate with one subexon.
		subexonCnt = subexons.size() ;
		for ( i = 0 ; i < subexonCnt ; ++i )
		{	
			struct _subexon &se = subexons[i] ;
			//printf( "hi1 %d: %d %d\n", i, se.prevCnt, se.prev[0] ) ;
			int cnt = 0 ;

			// due to filter, we may not fully match the coordinate and the subexon
			int bound = 0 ;
			if ( se.prevCnt > 0 )
				bound = se.prev[0] ;
			for ( j = i - 1, k = 0 ; k < se.prevCnt && j >= 0 && subexons[j].end >= bound ; --j )
			{
				//printf( " %d %d: %d %d\n", j, k, se.prev[ se.prevCnt - 1 - k], subexons[j].end ) ;
				if ( subexons[j].end == se.prev[se.prevCnt - 1 - k] ) // notice the order is reversed
				{
					se.prev[se.prevCnt - 1 - cnt] = j ;
					++k ;
					++cnt ;
				}
				else if ( subexons[j].end < se.prev[ se.prevCnt - 1 - k ] ) // the corresponding subexon gets filtered.
				{
					++k ;
					++j ; // counter the --j in the loop
				}
			}
			//printf( "hi2 %d : %d\n", i, se.prevCnt ) ;
			// shft the list
			for ( j = 0, k = se.prevCnt - cnt ; j < cnt ; ++j, ++k )
			{
				se.prev[j] = se.prev[k] ;
			}
			se.prevCnt = cnt ;
			cnt = 0 ;
			if ( se.nextCnt > 0 )
				bound = se.next[ se.nextCnt - 1] ;
			for ( j = i + 1, k = 0 ; k < se.nextCnt && j < subexonCnt && subexons[j].start <= bound ; ++j )
			{
				if ( subexons[j].start == se.next[k] )
				{
					se.next[cnt] = j ; // cnt is always less than k, so we don't need to worry about overwrite.
					++k ;
					++cnt ;
				}
				else if ( subexons[j].start > se.next[k] ) 
				{
					++k ;
					--j ;
				}
			}
			se.nextCnt = cnt ;
		}

		// Adjust the coordinate
		int seCnt = subexons.size() ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			--subexons[i].start ;
			--subexons[i].end ;
		}
		rewind( fpSubexon ) ;
		
		// Adjust the classifier for hard boundary, if there is a overhang attached to that region.
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( subexons[i].leftType == 1 && subexons[i].leftClassifier < 1 )
			{
				for ( j = i - 1 ; j >= 0 ; --j )
					if ( subexons[j].end < subexons[j + 1].start - 1 )
						break ;
				if ( subexons[j + 1].leftType == 0 )
					subexons[i].leftClassifier = 1 ;
			}
			if ( subexons[i].rightType == 2 && subexons[i].rightClassifier < 1 )
			{
				for ( j = i + 1 ; j < seCnt ; ++j )
					if ( subexons[j].start > subexons[j - 1].end + 1 )
						break ;
				if ( subexons[j - 1].rightType == 0 )
					subexons[i].rightClassifier = 1 ;
			}
		}

		// For the region of mixture of plus and minus strand subexons, if there is
		// no overhang attached to it, we need to let the hard boundary be a candidate terminal sites.
		for ( i = 0 ; i < seCnt ; )
		{
			// [i,j) is a region
			int support[2] = {0, 0} ; // the index, 0 is for minus strand, 1 is for plus strand
			for ( j = i + 1 ; j < seCnt ; ++j )
			{
				if ( subexons[j].start > subexons[j - 1].end + 1 )	
					break ;
			}

			for ( k = i ; k < j ; ++k )
			{
				if ( subexons[k].leftStrand != 0 )
					++support[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
				if ( subexons[k].rightStrand != 0 )
					++support[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
			}
			if ( support[0] == 0 || support[1] == 0 )
			{
				i = j ;
				continue ;
			}
			// a mixture region. 
			// We force a terminal site if we have only coming-in and no going-out introns.
			int leftSupport[2] = {0, 0}, rightSupport[2] = {0, 0};
			int l ;
			for ( k = i ; k < j ; ++k )
			{
				int cnt = subexons[k].prevCnt ; 
				if ( subexons[k].leftStrand != 0 )
					for ( l = 0 ; l < cnt ; ++l )
						if ( subexons[k].prev[l] < i )
						{
							++leftSupport[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
							break ;
						}
				cnt = subexons[k].nextCnt ; 
				if ( subexons[k].rightStrand != 0 )
					for ( l = 0 ; l < cnt ; ++l )
						if ( subexons[k].next[l] >= j )
						{
							++rightSupport[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
							break ;
						}
			}

			if ( ( ( leftSupport[0] > 0 && rightSupport[0] == 0 ) || 
				( leftSupport[1] > 0 && rightSupport[1] == 0 ) ) &&
				subexons[j - 1].rightType != 0 )
			{
				subexons[j - 1].rightClassifier = 0 ;
			}

			if ( ( ( leftSupport[0] == 0 && rightSupport[0] > 0 ) || 
				( leftSupport[1] == 0 && rightSupport[1] > 0 ) ) &&
				subexons[j - 1].leftType != 0 )
			{
				subexons[j - 1].leftClassifier = 0 ;
			}

			i = j ;
		}

		this->classifierThreshold = classifierThreshold ;

		usedGeneId = baseGeneId = 0 ;
	} 
	
	static bool IsSameStrand( int a, int b )
	{
		if ( a == 0 || b == 0 )
			return true ;
		if ( a != b )
			return false ;
		return true ;
	}
	// Parse the input line
	static int InputSubexon( char *in, Alignments &alignments, struct _subexon &se, bool needPrevNext = false )
	{
		int i ;
		char chrName[50] ;
		char ls[3], rs[3] ;	
		sscanf( in, "%s %d %d %d %d %s %s %lf %lf %lf %lf %lf", chrName, &se.start, &se.end, &se.leftType, &se.rightType, ls, rs,
				&se.avgDepth, &se.leftRatio, &se.rightRatio, 
				&se.leftClassifier, &se.rightClassifier ) ;	
		se.chrId = alignments.GetChromIdFromName( chrName ) ;
		se.nextCnt = se.prevCnt = 0 ;
		se.next = se.prev = NULL ;
		se.lcCnt = se.rcCnt = 0 ;

		if ( ls[0] == '+' )
			se.leftStrand = 1 ;
		else if ( ls[0] == '-' )
			se.leftStrand = -1 ;
		else
			se.leftStrand = 0 ;

		if ( rs[0] == '+' )
			se.rightStrand = 1 ;
		else if ( rs[0] == '-' )
			se.rightStrand = -1 ;
		else
			se.rightStrand = 0 ;

		if ( needPrevNext )
		{
			char *p = in ;
			// Locate the offset for prevCnt
			for ( i = 0 ; i <= 11 ; ++i )
			{
				p = strchr( p, ' ' ) ;
				++p ;
			}

			sscanf( p, "%d", &se.prevCnt ) ;
			p = strchr( p, ' ' ) ;
			++p ;
			se.prev = new int[ se.prevCnt ] ;
			for ( i = 0 ; i < se.prevCnt ; ++i )
			{
				sscanf( p, "%d", &se.prev[i] ) ;
				p = strchr( p, ' ' ) ;
				++p ;
			}

			sscanf( p, "%d", &se.nextCnt ) ;
			p = strchr( p, ' ' ) ;
			++p ;
			se.next = new int[ se.nextCnt ] ;
			for ( i = 0 ; i < se.nextCnt ; ++i )
			{
				sscanf( p, "%d", &se.next[i] ) ;
				p = strchr( p, ' ' ) ;
				++p ;
			}
			
		}
		return 1 ;
	}
	
	int GetGeneIntervalIdx( int startIdx, int &endIdx, int timeStamp ) ;

	//@return: the number of intervals found
	int ComputeGeneIntervals() ;
	
	// Return a list of subexons in that interval and in retList the id of subexon 
	// should be adjusted to start from 0.
	int ExtractSubexons( int startIdx, int endIdx, struct _subexon *retList ) ;
} ;

#endif
