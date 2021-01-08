#include "SubexonGraph.hpp"

void SubexonGraph::GetGeneBoundary( int tag, int &boundary, int timeStamp )
{
	if ( visit[tag] == timeStamp )	
		return ;
	//printf( "%d %d\n", tag, timeStamp ) ;
	visit[tag] = timeStamp ;
	if ( subexons[tag].end > boundary )
		boundary = subexons[tag].end ;
	//if ( subexons[tag].start == 2858011 )
	//	printf( "%d: %d %d\n", tag, subexons[tag].nextCnt, subexons[tag].prevCnt) ;
	int i ;
	int cnt = subexons[tag].nextCnt ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		//printf( "next of %d: %d %d\n", tag, i, subexons[tag].next[i] ) ;
		GetGeneBoundary( subexons[tag].next[i], boundary, timeStamp ) ;
	}
}

int SubexonGraph::GetGeneIntervalIdx( int startIdx, int &endIdx, int timeStamp )
{
	int i ;
	int seCnt = subexons.size() ;
	if ( startIdx >= seCnt )
		return -1 ;
	int farthest = -1 ;
	GetGeneBoundary( startIdx, farthest, timeStamp ) ;

	for ( i = startIdx + 1 ; i < seCnt ; ++i )
	{
		if ( subexons[i].start > farthest || subexons[i].chrId != subexons[ startIdx ].chrId )								
			break ;
		
		GetGeneBoundary( i, farthest, timeStamp ) ;
	}
	endIdx = i - 1 ;

	return endIdx ;
}

int SubexonGraph::ComputeGeneIntervals()
{
	int i, cnt ;
	int seCnt = subexons.size() ;
	visit = new int[seCnt] ;
	memset( visit, -1, sizeof( int ) * seCnt ) ;
	int tag = 0 ;
	cnt = 0 ;
	while ( 1 )		
	{
		struct _geneInterval ngi ;
		//printf( "%d %d %d\n", tag, subexons[tag].start + 1, subexons[tag].end + 1 ) ;
		if ( GetGeneIntervalIdx( tag, ngi.endIdx, cnt ) == -1 )
			break ;
		++cnt ;
		ngi.startIdx = tag ;
		ngi.start = subexons[ ngi.startIdx ].start ;
		ngi.end = subexons[ ngi.endIdx ].end ;
		
		tag = ngi.endIdx + 1 ;
		// Adjust the extent
		// Adjust the start 
		if ( subexons[ ngi.startIdx ].leftStrand != 0 
			&& subexons[ngi.startIdx].leftStrand != subexons[ngi.startIdx ].rightStrand )
			// We should make sure that rightstrand is non-zero whenever left-strand is non-zero for the startIdx.
		{
			for ( i = ngi.startIdx ; i >= 0 ; --i )	
			{
				if ( ( subexons[i].leftType == 1 && subexons[i].leftClassifier < classifierThreshold ) // an end within the subexon
					|| ( subexons[i].leftType == 0 ) // probably a overhang subexon. It should be a subset of the criterion following.
				  	|| ( i > 0 && subexons[i - 1].end + 1 < subexons[i].start ) ) // a gap. 
					break ;
			}
			ngi.start = subexons[i].start ;
		}

		// Adjust the end. 
		// And here, we also need to decide wether we need to adjust "tag" or not, 
		// because the next interval might be overlap with current interval by the last subexon.
		// We solve the overlap genes now, so we DON'T need to adjust tag.
		if ( subexons[ ngi.endIdx ].rightStrand != 0 
			&& subexons[ngi.endIdx].leftStrand != subexons[ngi.endIdx ].rightStrand )
		{
			for ( i = ngi.endIdx ; i < seCnt ; ++i )	
			{
				if ( ( subexons[i].rightType == 2 && subexons[i].rightClassifier < classifierThreshold ) // an end within the subexon
					|| ( subexons[i].rightType == 0 ) // probably a overhang subexon.
				  	|| ( i < seCnt - 1 && subexons[i].end + 1 < subexons[i + 1].start ) ) // a gap
					break ;
			}
			ngi.end = subexons[i].end ;
			
			/*if ( subexons[ ngi.endIdx ].rightType == 2 )
			{
				for ( i = ngi.endIdx ; i >= ngi.startIdx ; --i )
				{	
					if ( subexons[i].leftType == 1 )
						break ;
				}
				// The last region overlapps.
				if ( i >= ngi.startIdx && subexons[i].leftStrand != subexons[ ngi.endIdx ].rightStrand )
					--tag ;
			}*/
		}
		geneIntervals.push_back( ngi ) ;
	}
	delete[] visit ;

	return cnt ;
}

int SubexonGraph::ExtractSubexons( int startIdx, int endIdx, struct _subexon *retList )
{
	int i, j, k ;
	int cnt = endIdx - startIdx + 1 ;
	//printf( "%s: %d %d %d\n", __func__, startIdx, endIdx, cnt ) ;
	for ( i = 0 ; i < cnt ; ++i )
	{
		retList[i] = subexons[i + startIdx] ;
		retList[i].geneId = -1 ;
		retList[i].prev = new int[ retList[i].prevCnt ]	;
		retList[i].next = new int[ retList[i].nextCnt ] ;
		
		for ( j = 0 ; j < retList[i].prevCnt ; ++j )
			retList[i].prev[j] = subexons[i + startIdx].prev[j] - startIdx ;
		for ( j = 0 ; j < retList[i].nextCnt ; ++j )
			retList[i].next[j] = subexons[i + startIdx].next[j] - startIdx ;
		
		for ( j = 0, k = 0 ; j < retList[i].prevCnt ; ++j )
			if ( retList[i].prev[j] >= 0 && retList[i].prev[j] < cnt )
			{
				retList[i].prev[k] = retList[i].prev[j] ;
				++k ;
			}
		retList[i].prevCnt = k ;

		for ( j = 0, k = 0 ; j < retList[i].nextCnt ; ++j )
			if ( retList[i].next[j] >= 0 && retList[i].next[j] < cnt )
			{
				retList[i].next[k] = retList[i].next[j] ;
				++k ;
			}
		retList[i].nextCnt = k ;
	}
	UpdateGeneId( retList, cnt ) ;
	return cnt ;	
}

void SubexonGraph::SetGeneId( int tag, int strand, struct _subexon *subexons, int seCnt, int id )
{
	if ( subexons[tag].geneId != -1 && subexons[tag].geneId != -2 )
	{
		if ( subexons[tag].geneId != id ) // a subexon may belong to more than one gene.
		{
			//printf( "Set -2, %d: %d %d %d %d\n", id, tag, subexons[tag].geneId, subexons[tag].start + 1, strand ) ;
			subexons[tag].geneId = -2 ; 
		}
		else
			return ;
		// There is no need to terminate at the ambiguous exon, the strand will prevent
		//    us from overwriting previous gene ids.
		//return ; 
	}
	else if ( subexons[tag].geneId == -2 )
		return ;
	//printf( "%d: %d %d %d %d\n", id, tag, subexons[tag].geneId, subexons[tag].start + 1, strand ) ;
	int i ;
	if ( subexons[tag].geneId != -2 )
		subexons[ tag ].geneId = id ;
	int cnt = subexons[tag].nextCnt ;
	// Set through the introns.
	if ( IsSameStrand( strand, subexons[tag].rightStrand ) )
	{
		for ( i = 0 ; i < cnt ; ++i )
			if ( subexons[ subexons[tag].next[i] ].start > subexons[tag].end + 1 )
				SetGeneId( subexons[tag].next[i], strand, subexons, seCnt, id ) ;
	}

	cnt = subexons[tag].prevCnt ;
	if ( IsSameStrand( strand, subexons[tag].leftStrand ) )
	{
		for ( i = 0 ; i < cnt ; ++i )
			if ( subexons[ subexons[tag].prev[i] ].end < subexons[tag].start - 1 )
				SetGeneId( subexons[tag].prev[i], strand, subexons, seCnt, id ) ;
	}

	// Set through the adjacent subexons.
	if ( tag < seCnt - 1 && subexons[tag + 1].start == subexons[tag].end + 1 )
	{
		SetGeneId( tag + 1, strand, subexons, seCnt, id ) ;
	}

	if ( tag > 0 && subexons[tag].start - 1 == subexons[tag - 1].end )
	{
		SetGeneId( tag - 1, strand, subexons, seCnt, id ) ;
	}
	
}

void SubexonGraph::UpdateGeneId( struct _subexon *subexons, int seCnt )
{
	int i ;
	baseGeneId = usedGeneId ;
	int lastMinusStrandGeneId = -1 ;
	for ( int strand = -1 ; strand <= 1 ; strand +=2 )
	{
		for ( i = 0 ; i < seCnt ; ++i )
		{
			//printf( "%d (%d %d) %d.\n", i, subexons[i].start + 1, subexons[i].end + 1, subexons[i].geneId ) ;
			if ( ( subexons[i].geneId == -1 && ( ( strand == 1 && subexons[i].rightStrand == 0 ) || subexons[i].rightStrand == strand ) )
					|| ( strand == 1 && baseGeneId <= subexons[i].geneId && subexons[i].geneId <= lastMinusStrandGeneId && subexons[i].rightStrand == strand ) )
			{
				SetGeneId( i, strand, subexons, seCnt, usedGeneId ) ;
				if ( strand == -1 )
					lastMinusStrandGeneId = usedGeneId ;
				++usedGeneId ;
			}
		}
	}

	for ( i = 0 ; i < seCnt ; ++i )
		if ( subexons[i].leftType == 0 && subexons[i].rightType == 0 )
		{
			subexons[i].geneId = usedGeneId ;
			++usedGeneId ;
		}
	// Put base and usedGeneId in lcCnt, rcCnt field.
	for ( i = 0 ; i < seCnt ; ++i )	
	{
		subexons[i].lcCnt = baseGeneId ;
		subexons[i].rcCnt = usedGeneId ;
	}

	/*for ( i = 0 ; i < seCnt ; ++i ) 
	  {
	  printf( "geneId %d: %d-%d %d\n", i, subexons[i].start + 1, subexons[i].end + 1, subexons[i].geneId ) ;
	  }
	printf("%d %d\n", baseGeneId, usedGeneId ) ;*/
}
