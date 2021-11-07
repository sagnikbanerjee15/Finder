#include "TranscriptDecider.hpp"

void TranscriptDecider::OutputTranscript( int sampleId, struct _subexon *subexons, struct _transcript &transcript )
{
	int i, j ;
	// determine the strand
	std::vector<int> subexonInd ;
	transcript.seVector.GetOnesIndices( subexonInd ) ;

	// Determine the strand
	char strand[2] = "." ;
	int size = subexonInd.size() ;
	if ( size > 1 )
	{
		// locate the intron showed up in this transcript.
		for ( i = 0 ; i < size - 1 ; ++i )
		{
			/*int nextCnt = subexons[ subexonInd[i] ].nextCnt ;
			if ( nextCnt == 0 )
				continue ;

			for ( j = 0 ; j < nextCnt ; ++j )
			{
				int a = subexons[ subexonInd[i] ].next[j] ;
				if ( subexonInd[i + 1] == a 
					&& subexons[ subexonInd[i] ].end + 1 < subexons[a].start ) // avoid the case like ..(...[...
				{
					break ;
				}
			}
			if ( j < nextCnt )*/

			if ( subexons[ subexonInd[i] ].end + 1 < subexons[ subexonInd[i + 1] ].start )
			{
				if ( subexons[ subexonInd[i] ].rightStrand == 1 )
					strand[0] = '+' ;
				else if (  subexons[ subexonInd[i] ].rightStrand == -1 )
					strand[0] = '-' ;
				break ;
			}
		}
	}

	// TODO: transcript_id
	char *chrom = alignments.GetChromName( subexons[0].chrId ) ;
	char prefix[10] = "" ;
	struct _subexon *catSubexons = new struct _subexon[ size + 1 ] ;
	// Concatenate adjacent subexons 
	catSubexons[0] = subexons[ subexonInd[0] ] ;
	j = 1 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( subexons[ subexonInd[i] ].start == catSubexons[j - 1].end + 1 )
		{
			catSubexons[j - 1].end = subexons[ subexonInd[i] ].end ;
		}
		else		
		{
			catSubexons[j] = subexons[ subexonInd[i] ] ;
			++j ;
		}
	}
	size = j ;
	
	int gid = GetTranscriptGeneId( subexonInd, subexons ) ;
	if ( 0 ) //numThreads <= 1 )
	{
		fprintf( outputFPs[sampleId], "%s\tCLASSES\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; Abundance \"%.6lf\";\n",
				chrom, catSubexons[0].start + 1, catSubexons[size - 1].end + 1, strand,
				prefix, chrom, gid,
				prefix, chrom, gid, transcriptId[ gid - baseGeneId ], transcript.FPKM ) ;
		for ( i = 0 ; i < size ; ++i )
		{
			fprintf( outputFPs[ sampleId ], "%s\tCLASSES\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s%s.%d\"; "
					"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; Abundance \"%.6lf\"\n",
					chrom, catSubexons[i].start + 1, catSubexons[i].end + 1, strand,
					prefix, chrom, gid,
					prefix, chrom, gid, transcriptId[ gid - baseGeneId ],
					i + 1, transcript.FPKM ) ;
		}
	}
	else
	{
		struct _outputTranscript t ;	
		int len = 0 ;
		t.chrId = subexons[0].chrId ;
		t.geneId = gid ;
		t.transcriptId = transcriptId[ gid - baseGeneId ] ;
		t.FPKM = transcript.FPKM ;
		t.sampleId = sampleId ;
		t.exons = new struct _pair32[size] ;
		for ( i = 0 ; i < size ; ++i )
		{
			t.exons[i].a = catSubexons[i].start + 1 ;
			t.exons[i].b = catSubexons[i].end + 1 ;
			len += t.exons[i].b - t.exons[i].a + 1 ;
		}
		t.cov = transcript.abundance * alignments.readLen / len ;
		t.ecnt = size ;
		t.strand = strand[0] ;
		//printf( "%lf\n", transcript.correlationScore ) ;

		if ( numThreads > 1 )
			outputHandler->Add( t ) ;
		else
			outputHandler->Add_SingleThread( t ) ;
	}
	++transcriptId[ gid - baseGeneId ] ;

	delete[] catSubexons ;
}

int TranscriptDecider::GetFather( int f, int *father ) 
{
	if ( father[f] != f )
		return father[f] = GetFather( father[f], father ) ;
	return f ;
}

int TranscriptDecider::GetTranscriptGeneId( std::vector<int> &subexonInd, struct _subexon *subexons )
{
	int i ;
	int size = subexonInd.size() ;

	for ( i = 0 ; i < size ; ++i )
		if ( subexons[ subexonInd[i] ].geneId != -2  )
			return  subexons[ subexonInd[i] ].geneId ;
	
	// Some extreme case, where all the regions are mixture regions.
	for ( i = 0 ; i < size - 1 ; ++i )
		if ( subexons[ subexonInd[i] ].end + 1 < subexons[ subexonInd[i + 1] ].start )
		{
			return defaultGeneId[ ( subexons[ subexonInd[i] ].rightStrand + 1 ) / 2 ] ;
		}
	return defaultGeneId[0] ;
}

int TranscriptDecider::GetTranscriptGeneId( struct _transcript &t, struct _subexon *subexons )
{
	if ( subexons[ t.first ].geneId != -2 )
		return subexons[ t.first ].geneId ;
	if ( subexons[ t.last ].geneId != -2 )
		return subexons[ t.last ].geneId ;
	std::vector<int> subexonInd ;
	t.seVector.GetOnesIndices( subexonInd ) ;
	return GetTranscriptGeneId( subexonInd, subexons ) ;
}

void TranscriptDecider::InitTranscriptId()
{
	int i ;
	for ( i = 0 ; i < usedGeneId - baseGeneId ; ++i )
		transcriptId[i] = 0 ;
}

bool TranscriptDecider::IsStartOfMixtureStrandRegion( int tag, struct _subexon *subexons, int seCnt ) 
{
	int j, k ;
	int leftStrandCnt[2] = {0, 0}, rightStrandCnt[2] = {0, 0};
	for ( j = tag + 1 ; j < seCnt ; ++j )
		if ( subexons[j].start > subexons[j - 1].end + 1 )
			break ;
	
	for ( k = tag ; k < j ; ++k )
	{
		if ( subexons[k].leftStrand != 0 )
			++leftStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
		if ( subexons[k].rightStrand != 0 )
			++rightStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
	}

	if ( rightStrandCnt[0] > 0 && leftStrandCnt[0] == 0 && leftStrandCnt[1] > 0 )
		return true ;
	if ( rightStrandCnt[1] > 0 && leftStrandCnt[1] == 0 && leftStrandCnt[0] > 0 )
		return true ;
	return false ;
}

// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
// the partial compatible case (return 2) mostly likely happen in DP where we have partial transcript.
int TranscriptDecider::IsConstraintInTranscript( struct _transcript transcript, struct _constraint &c ) 
{
	//printf( "%d %d, %d %d\n", c.first, c.last, transcript.first, transcript.last ) ;
	if ( c.first < transcript.first || c.first > transcript.last 
		|| !transcript.seVector.Test( c.first ) 
		|| ( !transcript.partial && !transcript.seVector.Test( c.last ) ) ) // no overlap or starts too early or some chosen subexons does not compatible
		return 0 ; 
	
	// Extract the subexons we should focus on.
	int s, e ;
	s = c.first ;
	e = c.last ;
	bool returnPartial = false ;
	if ( e > transcript.last ) // constraints ends after the transcript.
	{
		if ( transcript.partial )	
		{
			e = transcript.last ;
			returnPartial = true ;
		}
		else
			return 0 ;
	}
	/*printf( "%s: %d %d: (%d %d) (%d %d)\n", __func__, s, e,
	  transcript.seVector.Test(0), transcript.seVector.Test(1), 
	  c.vector.Test(0), c.vector.Test(1) ) ;*/

	compatibleTestVectorT.Assign( transcript.seVector ) ;
	//compatibleTestVectorT.MaskRegionOutsideInRange( s, e, transcript.first, transcript.last ) ;
	compatibleTestVectorT.MaskRegionOutside( s, e ) ;

	compatibleTestVectorC.Assign( c.vector ) ;
	if ( c.last > transcript.last )
	{
		//compatibleTestVectorC.MaskRegionOutsideInRange( s, e, c.first, c.last ) ;
		//compatibleTestVectorC.MaskRegionOutside( s, e ) ;
		compatibleTestVectorC.MaskRegionOutside( 0, e ) ; // Because the bits before s are already all 0s in C. 
	}
	/*printf( "after masking %d %d. %d %d %d %d:\n", s, e, transcript.first, transcript.last, c.first, c.last ) ;
	compatibleTestVectorT.Print() ;
	compatibleTestVectorC.Print() ; */
	// Test compatible.
	int ret = 0 ;
	if ( compatibleTestVectorT.IsEqual( compatibleTestVectorC ) )
	{
		if ( returnPartial )
			ret = 2 ;
		else
			ret = 1 ;
	}

	return ret ;
}

int TranscriptDecider::IsConstraintInTranscriptDebug( struct _transcript transcript, struct _constraint &c ) 
{
	//printf( "%d %d, %d %d\n", c.first, c.last, transcript.first, transcript.last ) ;
	if ( c.first < transcript.first || c.first > transcript.last ) // no overlap or starts too early.
		return 0 ; 
	printf( "hi\n" ) ;
	// Extract the subexons we should focus on.
	int s, e ;
	s = c.first ;
	e = c.last ;
	bool returnPartial = false ;
	if ( e > transcript.last ) // constraints ends after the transcript.
	{
		if ( transcript.partial )	
		{
			e = transcript.last ;
			returnPartial = true ;
		}
		else
			return 0 ;
	}
	/*printf( "%s: %d %d: (%d %d) (%d %d)\n", __func__, s, e,
	  transcript.seVector.Test(0), transcript.seVector.Test(1), 
	  c.vector.Test(0), c.vector.Test(1) ) ;*/

	compatibleTestVectorT.Assign( transcript.seVector ) ;
	compatibleTestVectorT.MaskRegionOutside( s, e ) ;

	compatibleTestVectorC.Assign( c.vector ) ;
	if ( e > transcript.last )
		compatibleTestVectorC.MaskRegionOutside( s, e ) ;
	/*printf( "after masking: (%d %d) (%d %d)\n", 
	  compatibleTestVectorT.Test(0), compatibleTestVectorT.Test(1), 
	  compatibleTestVectorC.Test(0), compatibleTestVectorC.Test(1) ) ;*/

	// Test compatible.
	int ret = 0 ;
	if ( compatibleTestVectorT.IsEqual( compatibleTestVectorC ) )
	{
		if ( returnPartial )
			ret = 2 ;
		else
			ret = 1 ;
	}
	compatibleTestVectorT.Print() ;
	compatibleTestVectorC.Print() ;
	printf( "ret=%d\n", ret ) ;
	return ret ;
}
int TranscriptDecider::SubTranscriptCount( int tag, struct _subexon *subexons, int *f )
{
	if ( f[tag] != -1 )
		return f[tag] ;

	int ret = 0 ;
	int i ;
	if ( subexons[tag].canBeEnd )
		ret = 1 ;
	for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
	{
		ret += SubTranscriptCount( subexons[tag].next[i], subexons, f ) ;
	}

	if ( ret == 0 )
		ret = 1 ;
	return f[tag] = ret ;
}

void TranscriptDecider::CoalesceSameTranscripts( std::vector<struct _transcript> &t ) 
{
	int i, k ;
	if ( t.size() == 0 )
		return ;

	std::sort( t.begin(), t.end(), CompSortTranscripts ) ;

	int size = t.size() ;
	k = 0 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( t[k].seVector.IsEqual( t[i].seVector ) )
		{
			t[k].abundance += t[i].abundance ;
			t[i].seVector.Release() ;
		}
		else
		{
			++k ;
			if ( i != k )
				t[k] = t[i] ;
		}
	}
	t.resize( k + 1 ) ;
}

void TranscriptDecider::EnumerateTranscript( int tag, int strand, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, std::vector<struct _transcript> &alltranscripts, int &atcnt )
{
	int i ;
	visit[ vcnt ] = tag ;
	//printf( "%s: %d %d %d %d. %d %d\n", __func__, vcnt, tag, subexons[tag].nextCnt, strand, subexons[tag].start, subexons[tag].end ) ;
	// Compute the correlation score
	double minCor = correlationScore ;
	for ( i = 0 ; i < vcnt ; ++i )
	{
		double tmp = correlation.Query( visit[i], visit[vcnt] ) ;
		if ( tmp < minCor ) 
			minCor = tmp ;
	}

	if ( subexons[tag].canBeEnd )
	{
		struct _transcript &txpt = alltranscripts[atcnt] ;
		for ( i = 0 ; i <= vcnt ; ++i )
			txpt.seVector.Set( visit[i] ) ;

		txpt.first = visit[0] ;
		txpt.last = visit[vcnt] ;
		txpt.partial = false ;
		txpt.correlationScore = minCor ;
		
		//printf( "%lf %d %d ", txpt.correlationScore, vcnt, visit[0] ) ;
		//txpt.seVector.Print() ;
		++atcnt ;
	}

	for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
	{
		int a = subexons[tag].next[i] ;
		if ( !SubexonGraph::IsSameStrand( subexons[tag].rightStrand, strand ) 
			&& subexons[a].start > subexons[tag].end + 1 )
			continue ;
		int backupStrand = strand ;
		if ( subexons[a].start > subexons[tag].end + 1 && strand == 0 )
			strand = subexons[tag].rightStrand ;
		EnumerateTranscript( subexons[tag].next[i], strand, visit, vcnt + 1, subexons, correlation, minCor, alltranscripts, atcnt ) ;
		strand = backupStrand ;
	}
}

void TranscriptDecider::SearchSubTranscript( int tag, int strand, int parents[], int pcnt, struct _dp &pdp, int visit[], int vcnt, int extends[], int extendCnt,
std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) 
{
	int i ;
	int size ;
	double cover ;
	bool keepSearch = true ;
	bool belowMin = false ;
	
	struct _subexon *subexons = attr.subexons ;

	visit[vcnt] = tag ;
	++vcnt ;
	struct _dp visitdp ;
	
	visitdp.cover = -1 ;
	
	struct _transcript &subTxpt = attr.bufferTxpt ;
	subTxpt.seVector.Reset() ;
	for ( i = 0 ; i < pcnt ; ++i ) 
		subTxpt.seVector.Set( parents[i] ) ;
	subTxpt.first = parents[0] ;
	subTxpt.last = parents[ pcnt - 1] ;
	for ( i = 0 ; i < vcnt ; ++i )
		subTxpt.seVector.Set( visit[i] ) ;
	subTxpt.last = visit[ vcnt - 1 ] ;
	subTxpt.partial = true ;

	// Adjust the extendsCnt
	/*printf( "%s: %d %d %d\n", __func__, vcnt , extendCnt, extends[ extendCnt - 1] ) ;
	subTxpt.seVector.Print() ;
	tc[extends[extendCnt - 1]].vector.Print() ;
	printf( "Adjust extend:\n") ;*/
	for ( i = extendCnt - 1 ; i >= 0 ; --i )
	{
		if ( tc[ extends[i] ].last <= tag || ( tc[ extends[i] ].vector.Test( tag ) && IsConstraintInTranscript( subTxpt, tc[ extends[i] ] ) != 0 ) )
			break ;
	}
	extendCnt = i + 1 ;
	
	// If the extension ends.
	subTxpt.partial = false ;
	if ( subexons[tag].nextCnt > 0 && ( extendCnt == 0 || tag >= tc[ extends[ extendCnt - 1 ] ].last ) )
	{
		// Solve the subtranscript beginning with visit.
		// Now we got the optimal transcript for visit. 
		visitdp = SolveSubTranscript( visit, vcnt, strand, tc, tcStartInd, attr ) ;	
		keepSearch = false ;
	}
	//printf( "%s %d %d: visitdp.cover=%lf\n", __func__, parents[0], tag, visitdp.cover ) ;

	// the constraints across the parents and visit.
	size = tc.size() ;
	if ( visitdp.cover >= 0 )
	{
		cover = visitdp.cover ;
		// Reset the subTxpt, since its content is modofitied in SolveSubTxpt called above.
		subTxpt.seVector.Reset() ;
		for ( i = 0 ; i < pcnt ; ++i ) 
			subTxpt.seVector.Set( parents[i] ) ;
		subTxpt.seVector.Or( visitdp.seVector ) ;
		subTxpt.first = parents[0] ;
		subTxpt.last = visitdp.last ;
		subTxpt.partial = false ;

		if ( !attr.forAbundance && attr.minAbundance > 0 )
		{
			for ( i = 0 ; i < pcnt - 1 ; ++i )
			{
				if ( attr.uncoveredPair.find( parents[i] * attr.seCnt + parents[i + 1] ) != attr.uncoveredPair.end() )
					belowMin = true ;
			}
			for ( i = -1 ; i < vcnt - 1 ; ++i )
			{
				if ( i == -1 && pcnt >= 1 )
				{
					if ( attr.uncoveredPair.find( parents[pcnt - 1]  * attr.seCnt + visit[0] ) != attr.uncoveredPair.end() )
						belowMin = true ;
				}
				else
				{
					if ( attr.uncoveredPair.find( visit[i] * attr.seCnt + visit[i + 1] ) != attr.uncoveredPair.end() )
						belowMin = true ;	
				}
			}
			if ( attr.forAbundance && belowMin )
				cover = 1e-6 ;	
		}

		for ( i = tcStartInd ; i < size ; ++i )		
		{
			if ( tc[i].first > parents[ pcnt - 1] )
				break ;

			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					belowMin = true ;
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;

				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )
						cover = tc[i].normAbund ;
				}
				else
				{
					++cover ;
				}
			}
		}
		if ( belowMin && pdp.cover == -1 )	
		{
			pdp.cover = -2 ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
		else if ( cover > pdp.cover )
		{
			pdp.cover = cover ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
	}
	else if ( visitdp.cover == -2 && pdp.cover == -1 ) // no valid extension from visit 
	{
		subTxpt.seVector.Reset() ;
		for ( i = 0 ; i < pcnt ; ++i ) 
			subTxpt.seVector.Set( parents[i] ) ;
		subTxpt.seVector.Or( visitdp.seVector ) ;
		subTxpt.first = parents[0] ;
		subTxpt.last = visitdp.last ;

		pdp.cover = -2 ;
		pdp.seVector.Assign( subTxpt.seVector ) ;
		pdp.first = subTxpt.first ;
		pdp.last = subTxpt.last ;
		pdp.strand = strand ;
	}

	if ( subexons[tag].canBeEnd && ( visitdp.cover < 0 || attr.forAbundance ) ) 
	// This works is because that the extension always covers more constraints. So we only go this branch if the extension does not work
	// and it goes this branch if it violates minAbundance
	// But we need to go here when we want to compute the maxAbundance transcript.
	// This part also works as the exit point of the recurive function.
	{
		bool belowMin = false ;
		subTxpt.seVector.Reset() ;
		for ( i = 0 ; i < pcnt ; ++i ) 
			subTxpt.seVector.Set( parents[i] ) ;
		for ( i = 0 ; i < vcnt ; ++i )
			subTxpt.seVector.Set( visit[i] ) ;
		subTxpt.first = parents[0] ;
		subTxpt.last = visit[ vcnt - 1] ;
		subTxpt.partial = false ;

		cover = 0 ;
		if ( attr.forAbundance || attr.minAbundance > 0 )
		{
			for ( i = 0 ; i < pcnt - 1 ; ++i )
			{
				if ( attr.uncoveredPair.find( parents[i] * attr.seCnt + parents[i + 1] ) != attr.uncoveredPair.end() )
					belowMin = true ;
			}
			for ( i = -1 ; i < vcnt - 1 ; ++i )
			{
				if ( i == -1 && pcnt >= 1 )
				{
					if ( attr.uncoveredPair.find( parents[pcnt - 1]  * attr.seCnt + visit[0] ) != attr.uncoveredPair.end() )
						belowMin = true ;
				}
				else
				{
					if ( attr.uncoveredPair.find( visit[i] * attr.seCnt + visit[i + 1] ) != attr.uncoveredPair.end() )
						belowMin = true ;	
				}
			}

			//if ( belowMin == true )
			//	printf( "turned belowMin. %d. %d %d: %d %d %d\n", attr.uncoveredPair.size(), pcnt, vcnt, parents[0], visit[0], visit[ vcnt - 1] ) ;

			if ( attr.forAbundance && belowMin )
				cover = 1e-6 ;	
		}

		for ( i = tcStartInd ; i < size ; ++i )		
		{ 
			// note that the value is parents[ pcnt - 1], because  
			// in above the part of "visit" is computed in SolveSubTranscript( visit ).
			if ( tc[i].first > visit[ vcnt - 1] )
				break ;
			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					belowMin = true ;
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;
				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )	
						cover = tc[i].normAbund ;
				}
				else
				{
					++cover ;
				}
			}
		}

		if ( belowMin && pdp.cover == -1 )	
		{
			pdp.cover = -2 ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
		else if ( cover > pdp.cover )
		{
			pdp.cover = cover ;
			pdp.seVector.Assign( subTxpt.seVector ) ;
			pdp.first = subTxpt.first ;
			pdp.last = subTxpt.last ;
			pdp.strand = strand ;
		}
	}
	//printf( "%s %d: pdp.cover=%lf\n", __func__, tag, pdp.cover ) ;

	// keep searching.
	if ( keepSearch )	
	{
		for ( i = 0 ; i < subexons[tag].nextCnt ; ++i )
		{
			int b = subexons[tag].next[i] ;
			if ( ( SubexonGraph::IsSameStrand( subexons[tag].rightStrand, strand ) 
				&& SubexonGraph::IsSameStrand( subexons[b].leftStrand, strand ) ) ||
					subexons[b].start == subexons[tag].end + 1 )		
			{
				int backupStrand = strand ;
				if ( subexons[b].start > subexons[tag].end + 1 ) 
					strand = subexons[tag].rightStrand ;

				SearchSubTranscript( subexons[tag].next[i], strand, parents, pcnt, pdp, visit, vcnt, 
						extends, extendCnt, tc, tcStartInd, attr ) ;
				strand = backupStrand ;
			}
		}

	}

	return ;
}

struct _dp TranscriptDecider::SolveSubTranscript( int visit[], int vcnt, int strand, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) 
{
	int i ;
	int size ;
	/*printf( "%s: ", __func__ ) ;	
	for ( i = 0 ; i < vcnt ; ++i )
		printf( "%d ", visit[i] ) ;
	printf( ": %lf %d %d", attr.f1[ visit[0] ].cover, attr.f1[ visit[0] ].timeStamp, attr.timeStamp ) ;
	printf( "\n" ) ;*/
	// Test whether it is stored in dp 
	if ( vcnt == 1 )
	{
		if ( attr.f1[ visit[0] ].cover != -1 && attr.f1[ visit[0] ].strand == strand && ( attr.f1[ visit[0] ].timeStamp == attr.timeStamp  || 
			( attr.f1[ visit[0] ].minAbundance < attr.minAbundance && attr.f1[visit[0]].cover == -2 ) ) ) //even given lower minAbundance threshold, it fails
		{
			return attr.f1[ visit[0] ] ;
		}
	}
	else if ( vcnt == 2 && attr.f2 )
	{
		int a = visit[0] ;
		int b = visit[1] ;
		
		if ( attr.f2[a][b].cover != -1 && attr.f2[a][b].strand == strand && ( attr.f2[a][b].timeStamp == attr.timeStamp || 
			( attr.f2[a][b].minAbundance < attr.minAbundance && attr.f2[a][b].cover == -2 ) ) )
		{
			return attr.f2[a][b] ;
		}
	}
	else
	{
		int key = 0 ;	
		for ( i = 0 ; i < vcnt ; ++i )
			key = ( key * attr.seCnt + visit[i] ) % hashMax ;
		if ( key < 0 )
			key += hashMax ;

		if ( attr.hash[key].cover != -1 && attr.hash[key].cnt == vcnt && attr.hash[key].strand == strand && 
			( attr.hash[key].first == visit[0] )  &&
			( attr.hash[key].timeStamp == attr.timeStamp || 
				( attr.hash[key].minAbundance < attr.minAbundance && attr.hash[key].cover == -2 ) ) )
		{
			struct _transcript subTxpt = attr.bufferTxpt ;
			subTxpt.seVector.Reset() ;
			for ( i = 0 ; i < vcnt ; ++i )
				subTxpt.seVector.Set( visit[i] ) ;
			//subTxpt.seVector.Print() ;
			//attr.hash[key].seVector.Print() ;
			subTxpt.seVector.Xor( attr.hash[key].seVector ) ;
			subTxpt.seVector.MaskRegionOutside( visit[0], visit[ vcnt - 1] ) ;
			//printf( "hash test: %d %d\n", key, subTxpt.seVector.IsAllZero() ) ;
			if ( subTxpt.seVector.IsAllZero() )
			{
				return attr.hash[key] ;
			}
			
			// Can't use the code below, because vcnt is the header of subexons.
			/*for ( i = 0 ; i < vcnt ; ++i )
				if ( !attr.hash[key].seVector.Test( visit[i] ) )
					break ;
			if ( i >= vcnt )
				return attr.hash[key] ;*/
				
		}
	}
	// adjust tcStartInd 
	size = tc.size() ;
	for ( i = tcStartInd ; i < size ; ++i )
		if ( tc[i].first >= visit[0] )
			break ;
	tcStartInd = i ;


	struct _subexon *subexons = attr.subexons ;
	struct _dp visitdp ;
	visitdp.seVector.Init( attr.seCnt ) ;
	visitdp.cover = -1 ;

	struct _transcript &subTxpt = attr.bufferTxpt ;
	// This happens when it is called from PickTranscriptsByDP, the first subexon might be the end.
	subTxpt.seVector.Reset() ;	
	for ( i = 0 ; i < vcnt ; ++i )
		subTxpt.seVector.Set( visit[i] ) ;
	subTxpt.first = visit[0] ;
	subTxpt.last = visit[vcnt - 1] ;

	if ( subexons[ visit[vcnt - 1] ].canBeEnd )
	{
		subTxpt.partial = false ;
		double cover = 0 ;
		for ( i = tcStartInd ; i < size ; ++i )		
		{
			if ( tc[i].first > subTxpt.last )
				break ;

			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].normAbund <= attr.minAbundance )
				{
					cover = -2 ;
					break ;
				}

				if ( tc[i].abundance <= 0 )
					continue ;
				if ( attr.forAbundance )
				{
					if ( tc[i].normAbund < cover || cover == 0 )	
						cover = tc[i].normAbund ;
				}
				else
					++cover ;	
			}
		}

		visitdp.seVector.Assign( subTxpt.seVector ) ;		
		visitdp.cover = cover ;
		visitdp.first = subTxpt.first ;
		visitdp.last = subTxpt.last ;
		visitdp.strand = strand ;
	}
	
	// Now we extend.
	size = tc.size() ;
	int *extends = new int[tc.size() - tcStartInd + 1] ;
	int extendCnt = 0 ;
	subTxpt.partial = true ;
	for ( i = tcStartInd ; i < size ; ++i )
	{
		if ( tc[i].first > subTxpt.last )
			break ;
		if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 2 )
		{
			extends[extendCnt] = i ;
			++extendCnt ;
		}
	}

	// Sort the extend by the index of the last subexon.
	if ( extendCnt > 0 )
	{
		struct _pair32 *extendsPairs = new struct _pair32[extendCnt] ;

		for ( i = 0 ; i < extendCnt ; ++i )
		{
			extendsPairs[i].a = extends[i] ;
			extendsPairs[i].b = tc[ extends[i] ].last ;
		}
		qsort( extendsPairs, extendCnt, sizeof( struct _pair32 ), CompPairsByB ) ;

		for ( i = 0 ; i < extendCnt ; ++i )
			extends[i] = extendsPairs[i].a ;

		delete[] extendsPairs ;
	}

	size = subexons[ visit[vcnt - 1] ].nextCnt ;
	int nextvCnt = 1 ;
	if ( extendCnt > 0 && tc[ extends[ extendCnt - 1 ] ].last - visit[ vcnt - 1 ] > 1 )
		nextvCnt = tc[ extends[ extendCnt - 1 ] ].last - visit[ vcnt - 1 ] ;
	int *nextv = new int[ nextvCnt ] ;
	for ( i = 0 ; i < size ; ++i )
	{
		int a = visit[vcnt - 1] ;
		int b = subexons[a].next[i] ;
		if ( ( SubexonGraph::IsSameStrand( subexons[a].rightStrand, strand ) 
			&& SubexonGraph::IsSameStrand( subexons[b].leftStrand, strand ) )
			||
			subexons[b].start == subexons[a].end + 1 )		
		{
			int backupStrand = strand ;
			if ( subexons[b].start > subexons[a].end + 1 ) 
				strand = subexons[a].rightStrand ;
			SearchSubTranscript( subexons[ visit[vcnt - 1] ].next[i], strand, visit, vcnt, visitdp, nextv, 0, extends, extendCnt, tc, tcStartInd, attr ) ;		
			strand = backupStrand ;

		}
	}
	//printf( "%s %d(%d) %d %d %d: %lf\n", __func__, visit[0], subexons[ visit[vcnt - 1] ].canBeEnd, size, extendCnt, strand, visitdp.cover ) ;	
	delete[] nextv ;
	delete[] extends ;

	// store the result in the dp structure.
	// We return the structure stored in dp to simplify the memory access pattern.
	// In other words, we assume the structure returned from this function always uses the memory from attr.dp 
	if ( vcnt == 1 )
	{
		SetDpContent( attr.f1[ visit[0] ], visitdp, attr ) ;
		visitdp.seVector.Release() ;
		return attr.f1[ visit[0] ] ;
	}
	else if ( vcnt == 2 && attr.f2 )
	{
		SetDpContent( attr.f2[ visit[0] ][ visit[1] ], visitdp, attr ) ;
		visitdp.seVector.Release() ;
		return attr.f2[ visit[0] ][ visit[1] ] ;
	}
	else
	{
		int key = 0 ;	
		for ( i = 0 ; i < vcnt ; ++i )
			key = ( key * attr.seCnt + visit[i] ) % hashMax ;
		if ( key < 0 )
			key += hashMax ;

		//static int hashUsed = 0 ;
		//if (  attr.hash[key].cover == -1 )
		//	++hashUsed ;
		//printf( "%d/%d\n", hashUsed, HASH_MAX) ;
		//printf( "hash write: %d\n", key ) ;	
		SetDpContent( attr.hash[key], visitdp, attr ) ;	
		attr.hash[key].cnt = vcnt ;
		visitdp.seVector.Release() ;
		return attr.hash[key] ;
	}
}


void TranscriptDecider::PickTranscriptsByDP( struct _subexon *subexons, int seCnt, int iterBound, Constraints &constraints, SubexonCorrelation &correlation, struct _dpAttribute &attr, std::vector<struct _transcript> &alltranscripts )
{
	int i, j, k ;
	
	std::vector<struct _transcript> transcripts ;
	std::vector<struct _constraint> &tc = constraints.constraints ;
	int tcCnt = tc.size() ;
	int coalesceThreshold = 1024 ;

	//printf( "tcCnt=%d\n", tcCnt ) ;

	attr.timeStamp = 1 ;
	attr.bufferTxpt.seVector.Init( seCnt ) ;
	attr.subexons = subexons ;
	attr.seCnt = seCnt ;

	double maxAbundance = -1 ;
	// Initialize the dp data structure
	/*memset( attr.f1, -1, sizeof( struct _dp ) * seCnt ) ;
	for ( i = 0 ; i < seCnt ; ++i )
		memset( attr.f2[i], -1, sizeof( struct _dp ) * seCnt ) ;
	memset( attr.hash, -1, sizeof( struct _dp ) * HASH_MAX ) ;*/
	for ( i = 0 ; i < seCnt ; ++i )
		ResetDpContent( attr.f1[i] ) ;
	for ( i = 0 ; i < seCnt && attr.f2 ; ++i )
		for ( j = i ; j < seCnt ; ++j )
			ResetDpContent( attr.f2[i][j] ) ;
	for ( i = 0 ; i < hashMax ; ++i )
		ResetDpContent( attr.hash[i] ) ;

	// Set the uncovered pair
	attr.uncoveredPair.clear() ;
	BitTable bufferTable( seCnt ) ;
	k = 0 ;
	for ( i = 0 ; i < seCnt ; ++i )
	{	
		for ( ; k < tcCnt ; ++k )
		{
			if ( tc[k].last >= i )
				break ;
		}

		if ( k >= tcCnt || tc[k].first > i )
		{
			for ( j = 0 ; j < subexons[i].nextCnt ; ++j )
			{
				attr.uncoveredPair[i * seCnt + subexons[i].next[j] ]	= 1 ;
			}
			continue ;
		}
		
		for ( j = 0 ; j < subexons[i].nextCnt ; ++j )			
		{
			bool covered = false ;
			int l, n ;

			n = subexons[i].next[j] ;
			for ( l = k ; l < tcCnt ; ++l )
			{
				if ( tc[l].first > i )
					break ;
				if ( tc[l].vector.Test( i ) && tc[l].vector.Test( n ) )	
				{
					if ( n == i + 1 )
					{
						covered = true ;
						break ;
					}
					else
					{
						bufferTable.Assign( tc[l].vector ) ;
						bufferTable.MaskRegionOutside( i + 1, n - 1 ) ;
						if ( bufferTable.IsAllZero() )
						{
							covered = true ;
							break ;
						}
					}
				}
			}

			if ( !covered )
			{
				//printf( "set!: (%d: %d %d) (%d: %d %d)\n", i, subexons[i].start, subexons[i].end, n, subexons[n].start, subexons[n].end ) ;
				attr.uncoveredPair[ i * seCnt + n ] = 1 ;
			}
		}
	}
	bufferTable.Release() ;
		

	// Find the max abundance 
	attr.forAbundance = true ;
	attr.minAbundance = 0 ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].canBeStart )	
		{
			int visit[1] = {i} ;
			struct _dp tmp ;
			
			tmp = SolveSubTranscript( visit, 1, 0, tc, 0, attr ) ;
			
			if ( tmp.cover > maxAbundance )
				maxAbundance = tmp.cover ;
		}
	}
	//PrintLog( "maxAbundance=%lf", maxAbundance ) ;
	//exit( 1 ) ;

	// Pick the transcripts. Quantative Set-Cover
	// Notice that by the logic in SearchSubTxpt and SolveSubTxpt, we don't need to reinitialize the data structure.
	attr.forAbundance = false ;
	int *coveredTc = new int[tcCnt] ;
	int coveredTcCnt ;
	struct _dp maxCoverDp ;
	struct _dp bestDp ;
	std::map<double, struct _dp> cachedCoverResult ;

	maxCoverDp.seVector.Init( seCnt ) ;
	bestDp.seVector.Init( seCnt ) ;
	int iterCnt = 0 ;

	while ( 1 )
	{
		double bestScore ;
		
		// iterately assign constraints
		attr.minAbundance = 0 ;	
		
		// Find the best candidate transcript.
		bestDp.cover = -1 ;
		bestScore = -1 ;
		while ( 1 )
		{
			// iterate the change of minAbundance
			if ( cachedCoverResult.find( attr.minAbundance ) != cachedCoverResult.end() )
			{
				struct _dp tmp = cachedCoverResult[ attr.minAbundance ] ;
				SetDpContent( maxCoverDp, tmp, attr ) ;
			}
			else
			{
				maxCoverDp.cover = -1 ;
				++attr.timeStamp ;
				for ( i = 0 ; i < seCnt ; ++i )		
				{
					if ( subexons[i].canBeStart == false )
						continue ;
					int visit[1] = {i} ;
					struct _dp tmp ;
					tmp = SolveSubTranscript( visit, 1, 0, tc, 0, attr ) ;

					if ( tmp.cover > maxCoverDp.cover && tmp.cover > 0 )
					{
						SetDpContent( maxCoverDp, tmp, attr ) ;
					}
					//if ( subexons[i].start == 6870264 || subexons[i].start == 6872237 )
					//	printf( "%d: %lf\n", i, tmp.cover ) ;
				}

				if ( maxCoverDp.cover == -1 )
					break ;
				struct _dp ccr ;
				ccr.seVector.Init( seCnt ) ;
				SetDpContent( ccr, maxCoverDp, attr ) ;
				cachedCoverResult[ attr.minAbundance ] = ccr ; 
			}
			// the abundance for the max cover txpt.
			double min = -1 ;
			struct _transcript &subTxpt = attr.bufferTxpt ;
			subTxpt.seVector.Assign( maxCoverDp.seVector ) ;
			subTxpt.first = maxCoverDp.first ;
			subTxpt.last = maxCoverDp.last ;
			
			for ( i = 0 ; i < tcCnt ; ++i )
			{
				if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )	
				{
					if ( tc[i].normAbund < min || min == -1 )
						min = tc[i].normAbund ;
				}
			}

			if ( attr.minAbundance == 0 )
			{
				std::vector<int> subexonIdx ;
				maxCoverDp.seVector.GetOnesIndices( subexonIdx ) ;
				int size = subexonIdx.size() ;
				for ( i = 0 ; i < size - 1 ; ++i )
					if ( attr.uncoveredPair.find( subexonIdx[i] * seCnt + subexonIdx[i + 1] ) != attr.uncoveredPair.end() )
					{
						min = 1e-6 ;
						break ;
					}
			}
			
			double score = ComputeScore( maxCoverDp.cover, 1.0, min, maxAbundance, 0 ) ;
			if ( bestScore == -1 || score > bestScore )	
			{
				bestScore = score ;
				SetDpContent( bestDp, maxCoverDp, attr ) ;
			}
			else if ( score < bestScore )
			{
				if ( ComputeScore( maxCoverDp.cover, 1.0, maxAbundance, maxAbundance, 0 ) < bestScore )
					break ;
			}
			//PrintLog( "normAbund=%lf maxCoverDp.cover=%lf score=%lf timeStamp=%d", min, maxCoverDp.cover, score, attr.timeStamp ) ;
			attr.minAbundance = min ;
		} // end of iteration for minAbundance.

		if ( bestDp.cover == -1 )
			break ;
		// Assign the constraints.
		coveredTcCnt = 0 ;
		double update = -1 ;
		struct _transcript &subTxpt = attr.bufferTxpt ;
		subTxpt.seVector.Assign( bestDp.seVector ) ;
		subTxpt.first = bestDp.first ;
		subTxpt.last = bestDp.last ;
		subTxpt.partial = false ;
		for ( i = 0 ; i < tcCnt ; ++i )
		{
			if ( IsConstraintInTranscript( subTxpt, tc[i] ) == 1 )
			{
				if ( tc[i].abundance > 0 && 
					( tc[i].abundance < update || update == -1 ) )
				{
					update = tc[i].abundance ;		
				}
				coveredTc[ coveredTcCnt ] = i ;
				++coveredTcCnt ;
			}
			/*else
			{
				printf( "%d: ", i ) ;
				tc[i].vector.Print() ;
				if ( i == 127 )
				{
					printf( "begin debug:\n" ) ;
					IsConstraintInTranscriptDebug( subTxpt, tc[i] ) ;
				}
			}*/
		}
		update *= ( 1 + iterCnt / 50 ) ;//* ( 1 + iterCnt / 50 )  ; 

		//PrintLog( "%d: update=%lf %d %d. %d %d %d", iterCnt, update, coveredTcCnt, tcCnt, 
		//	bestDp.first, bestDp.last, subexons[ bestDp.first ].start ) ;
		//bestDp.seVector.Print() ;

		struct _transcript nt ;
		nt.seVector.Duplicate( bestDp.seVector ) ; 
		nt.first = bestDp.first ;
		nt.last = bestDp.last ;
		nt.partial = false ;
		nt.abundance = 0 ;
		for ( i = 0 ; i < coveredTcCnt ; ++i )
		{
			j = coveredTc[i] ;
			if ( tc[j].abundance > 0 )	
			{
				double tmp = ( tc[j].abundance > update ? update : tc[j].abundance ) ;
				tc[j].abundance -= tmp ;
				double factor = 1 ;

				nt.abundance += ( tc[j].support * update / tc[j].normAbund * factor ) ;

				if ( tc[j].abundance <= 0 )
				{
					std::vector<double> removeKey ;
					for ( std::map<double, struct _dp>::iterator it = cachedCoverResult.begin() ; it != cachedCoverResult.end() ; ++it )
					{
						subTxpt.seVector.Assign( it->second.seVector ) ;
						subTxpt.first = it->second.first ;
						subTxpt.last = it->second.last ;
						subTxpt.partial = false ;
						if ( IsConstraintInTranscript( subTxpt, tc[j] ) == 1 )
						{
							it->second.seVector.Release() ;
							removeKey.push_back( it->first ) ;
						}
					}
					int size = removeKey.size() ;
					int l ;
					for ( l = 0 ; l < size ; ++l )
						cachedCoverResult.erase( removeKey[l] ) ;
				}
			}

			if ( tc[j].abundance < 0 )
				tc[j].abundance = 0 ;
		}
		
		transcripts.push_back( nt ) ;
		if ( transcripts.size() >= transcripts.capacity() && (int)transcripts.size() >= coalesceThreshold )
		{
			CoalesceSameTranscripts( transcripts ) ;
			if ( transcripts.size() >= transcripts.capacity() / 2 )
				coalesceThreshold *= 2 ;
		}
		++iterCnt ;

		if ( iterCnt >= iterBound )
			break ;
	}
	CoalesceSameTranscripts( transcripts ) ;
	int size = transcripts.size() ;
	// Compute the correlation score
	for ( i = 0 ; i < size ; ++i )
	{
		std::vector<int> subexonInd ;
		transcripts[i].seVector.GetOnesIndices( subexonInd ) ;
		double cor = 2.0 ;
		int s = subexonInd.size() ;
		for ( j = 0 ; j < s ; ++j )
			for ( k = j + 1 ; k < s ; ++k )
			{
				double tmp = correlation.Query( subexonInd[j], subexonInd[k] ) ;
				if ( tmp < cor )					
					cor = tmp ;
			}
		if ( cor > 1 )
			cor = 0 ;
		transcripts[i].correlationScore = cor ;
	}

	// store the result
	for ( i = 0 ; i < size ; ++i )
		alltranscripts.push_back( transcripts[i] ) ;

	// Release the memory
	for ( std::map<double, struct _dp>::iterator it = cachedCoverResult.begin() ; it != cachedCoverResult.end() ; ++it )
	{
		it->second.seVector.Release() ;
	}
	attr.bufferTxpt.seVector.Release() ;

	delete[] coveredTc ;	
	maxCoverDp.seVector.Release() ;
	bestDp.seVector.Release() ;
}


// Add the preifx/suffix of transcripts to the list
void TranscriptDecider::AugmentTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &alltranscripts, int limit, bool extend )
{
	int i, j, k ;
	int size = alltranscripts.size() ;
	if ( size >= limit )
		return ;

	// Augment suffix, prefix transcripts
	for ( i = 0 ; i < size ; ++i )
	{
		std::vector<int> subexonIdx ;
		alltranscripts[i].seVector.GetOnesIndices( subexonIdx ) ;
		int seIdxCnt = subexonIdx.size() ;
		// suffix
		for ( j = 1 ; j < seIdxCnt ; ++j )
		{
			if ( subexons[ subexonIdx[j] ].canBeStart )
			{
				struct _transcript nt ;
				nt.first = subexonIdx[j] ;
				nt.last = alltranscripts[i].last ;
				nt.seVector.Duplicate( alltranscripts[i].seVector ) ;
				nt.seVector.MaskRegionOutside( nt.first, nt.last ) ;
				nt.partial = false ;
				nt.correlationScore = 0 ;
				nt.abundance = 0 ;
				nt.constraintsSupport = NULL ;
				
				alltranscripts.push_back( nt ) ;
				if ( alltranscripts.size() >= limit )
					return ;
			}
		}
		
		// prefix
		for ( j = 0 ; j < seIdxCnt - 1 ; ++j )
		{
			if ( subexons[ subexonIdx[j] ].canBeEnd )
			{
				struct _transcript nt ;
				nt.first = alltranscripts[i].first ;
				nt.last = subexonIdx[j] ;
				nt.seVector.Duplicate( alltranscripts[i].seVector ) ;
				nt.seVector.MaskRegionOutside( nt.first, nt.last ) ;
				nt.partial = false ;
				nt.correlationScore = 0 ;
				nt.abundance = 0 ;
				nt.constraintsSupport = NULL ;

				alltranscripts.push_back( nt ) ;
				if ( alltranscripts.size() >= limit )
					return ;
			}
		}

		if ( extend )
		{
			//Extentions right.
			for ( j = 0 ; j < seIdxCnt ; ++j )
			{
				if ( subexons[ subexonIdx[j] ].nextCnt > 1 )
				{
					for ( k = 0 ; k < subexons[ subexonIdx[j] ].nextCnt ; ++k )
					{
						int idx = subexons[ subexonIdx[j] ].next[k] ;

						if ( alltranscripts[i].seVector.Test( idx ) )
							continue ;
						int l ;
						std::vector<int> visited ;
						while ( 1 )
						{
							if ( subexons[idx].nextCnt > 1 || subexons[idx].prevCnt > 1 )
							{
								break ;
							}

							visited.push_back( idx ) ;
							if ( subexons[idx].canBeEnd && subexons[idx].nextCnt == 0 )
							{
								struct _transcript nt ;
								nt.first = alltranscripts[i].first ;
								nt.last = idx ;
								nt.seVector.Duplicate( alltranscripts[i].seVector ) ;
								nt.seVector.MaskRegionOutside( nt.first, subexonIdx[j] ) ;
								int visitedSize = visited.size() ;
								for ( l = 0 ; l < visitedSize ; ++l )
									nt.seVector.Set( visited[l] ) ;
								nt.partial = false ;
								nt.correlationScore = 0 ;
								nt.abundance = 0 ;
								nt.constraintsSupport = NULL ;

								alltranscripts.push_back( nt ) ;
								if ( alltranscripts.size() >= limit )
									return ;
							}

							if ( subexons[idx].nextCnt == 1 )
								idx = subexons[idx].next[0] ;
							else
								break ;
						}
					}
				}
			}

			// Extension towards left
			for ( j = 0 ; j < seIdxCnt ; ++j )
			{
				if ( subexons[ subexonIdx[j] ].prevCnt > 1 )
				{
					for ( k = 0 ; k < subexons[ subexonIdx[j] ].prevCnt ; ++k )
					{
						int idx = subexons[ subexonIdx[j] ].prev[k] ;

						if ( alltranscripts[i].seVector.Test( idx ) )
							continue ;
						int l ;
						std::vector<int> visited ;
						while ( 1 )
						{
							if ( subexons[idx].nextCnt > 1 || subexons[idx].prevCnt > 1 ) 
							{
								break ;
							}

							visited.push_back( idx ) ;
							if ( subexons[idx].canBeStart && subexons[idx].prevCnt == 0 )
							{
								struct _transcript nt ;
								nt.first = idx ;
								nt.last = alltranscripts[i].last ;
								nt.seVector.Duplicate( alltranscripts[i].seVector ) ;
								nt.seVector.MaskRegionOutside( subexonIdx[j], nt.last ) ;
								int visitedSize = visited.size() ;
								for ( l = 0 ; l < visitedSize ; ++l )
									nt.seVector.Set( visited[l] ) ;
								nt.partial = false ;
								nt.correlationScore = 0 ;
								nt.abundance = 0 ;
								nt.constraintsSupport = NULL ;

								alltranscripts.push_back( nt ) ;
								if ( alltranscripts.size() >= limit )
									return ;
							}

							if ( subexons[idx].prevCnt == 1 )
								idx = subexons[idx].prev[0] ;
							else
								break ;
						}
					}
				}
			}
		} // for if-extend
	}

	CoalesceSameTranscripts( alltranscripts ) ;
}

// Pick the transcripts from given transcripts.
void TranscriptDecider::PickTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &alltranscripts, Constraints &constraints, 
		SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) 
{
	int i, j, k ;
	std::vector<int> chosen ;
	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	int atcnt = alltranscripts.size() ;
	int tcCnt = tc.size() ; // transcript constraints
	int seCnt = 0 ;
	
	if ( tcCnt == 0 )
		return ;
	if ( atcnt > 0 )
		seCnt = alltranscripts[0].seVector.GetSize() ;
	else
		return ;

	double inf = -1 ; // infinity
	int coalesceThreshold = 1024 ;
	int *transcriptSeCnt = new int[ atcnt ] ;
	int *transcriptLength = new int[atcnt] ;
	double *transcriptAbundance = new double[atcnt] ; // the roughly estimated abundance based on constraints.
	double *avgTranscriptAbundance = new double[atcnt] ; // the average normAbund from the compatible constraints.

	BitTable *btable = new BitTable[ atcnt ] ; 
	//BitTable lowCovSubexon ; // force the abundance to 0 for the transcript contains the subexon.
	double *coveredPortion = new double[atcnt] ;

	memset( avgTranscriptAbundance, 0 ,sizeof( double ) * atcnt ) ;
	for ( i = 0 ; i < atcnt ; ++i )
		btable[i].Init( tcCnt ) ;
	for ( j = 0 ; j < tcCnt ; ++j )
	{
		int a = constraints.matePairs[j].i ;
		int b = constraints.matePairs[j].j ;

		if ( constraints.constraints[a].support > inf )
			inf = constraints.constraints[a].support ;
		if ( constraints.constraints[b].support > inf )
			inf = constraints.constraints[b].support ;

		if ( tc[j].normAbund > inf )
			inf = tc[j].normAbund ;

		tc[j].abundance = tc[j].normAbund ;
	}
	++inf ;
	bool btableSet = false ;
	for ( i = 0 ; i < atcnt ; ++i )
	{
		//printf( "correlation %d: %lf\n", i, alltranscripts[i].correlationScore ) ;
		/*for ( int l = 0 ; l < subexonInd.size() ; ++l )
		{
			for ( int m = l ; m < subexonInd.size() ; ++m )
				printf( "%lf ", seCorrelation.Query( l, m ) ) ;
			printf( "\n" ) ;
		}*/

		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int a = tc[j].i ;
			int b = tc[j].j ;

			//printf( "try set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
			//alltranscripts[i].seVector.Print() ;
			//constraints.constraints[a].vector.Print() ;
			//constraints.constraints[b].vector.Print() ;
			if ( IsConstraintInTranscript( alltranscripts[i], constraints.constraints[a] ) == 1 
					&& IsConstraintInTranscript( alltranscripts[i], constraints.constraints[b] ) == 1 )
			{
				//printf( "set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
				btable[i].Set( j ) ;
				btableSet = true ;
			}
		}
		transcriptSeCnt[i] = alltranscripts[i].seVector.Count() ;
	}
	if ( btableSet == false )
	{
		for ( i = 0 ; i < atcnt ; ++i )
			btable[i].Release() ;
		delete[] btable ;
		return ;
	}

	double maxAbundance = -1 ; // The abundance of the most-abundant transcript
	double *adjustScore = new double[atcnt] ;
	memset( adjustScore, 0, sizeof( double ) * atcnt ) ;
	if ( atcnt > 0 /*&& alltranscripts[0].abundance == -1*/ )
	{
		struct _pair32 *chain = new struct _pair32[seCnt] ;
		bool *covered = new bool[seCnt] ;
		bool *usedConstraints = new bool[constraints.constraints.size() ] ;
		std::vector<BitTable> togetherChain ; // those subexons is more likely to show up in the same transcript, like an IR with overhang, should be together to represent a 3'/5'-end

		/*lowCovSubexon.Init( seCnt ) ;
		double *avgDepth = new double[seCnt ] ;

		memset( avgDepth, 0, sizeof( double ) * seCnt ) ;
		int size = constraints.constraints.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			std::vector<int> subexonIdx ;
			constraints.constraints[i].GetOnesIndices( subexonIdx ) ;
			int seIdxCnt = subexonidx.size() ;
			for ( j = 0 ; j < seIdxCnt ; ++j )
				avgDepth[ subexonidx[j] ] += constraints.constraints[i].support ;
		}
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( avgDepth[i] * alignments.readLen / (double)( subexons[i].end - subexons[i].start + 1 ) < 1 )
		}*/

		struct _pair32 firstRegion, lastRegion ;
		
		for ( i = 0 ; i < seCnt ;  )
		{
			for ( j = i + 1 ; j < seCnt ; ++j ) 
			{
				if ( subexons[j].start > subexons[j - 1].end + 1 )
					break ;
			}

						
			int cnt = 0 ;
			for ( k = i ; k < j ; ++k )	
			{
				if ( ( subexons[k].leftType == 2 && subexons[k].rightType == 1 )
					|| ( subexons[k].leftType == 0 && subexons[k].rightType == 1 ) 
					|| ( subexons[k].leftType == 2 && subexons[k].rightType == 0 ) )
					++cnt ;
			}
			
			if ( cnt <= 1 )
			{
				i = j ;
				continue ;
			}

			BitTable tmpTable( seCnt ) ;
			for ( k = i ; k < j ; ++k )	
			{
				if ( ( subexons[k].leftType == 2 && subexons[k].rightType == 1 )
					|| ( subexons[k].leftType == 0 && subexons[k].rightType == 1 ) 
					|| ( subexons[k].leftType == 2 && subexons[k].rightType == 0 ) )
					tmpTable.Set( k ) ;
			}
			togetherChain.push_back( tmpTable ) ;
			i = j ;
		}

		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = inf ;
			int tag = -1 ;

			alltranscripts[i].abundance = 0 ;
			alltranscripts[i].constraintsSupport = new double[tcCnt] ;

			std::vector<int> subexonIdx ;
			alltranscripts[i].seVector.GetOnesIndices( subexonIdx ) ;
			int seIdxCnt = subexonIdx.size() ;
			transcriptLength[i] = 0 ;
			
			firstRegion.a = subexonIdx[0] ;
			for ( j = 1 ; j < seIdxCnt ; ++j )
			{
				if ( subexons[ subexonIdx[j] ].start > subexons[ subexonIdx[j - 1] ].end + 1 )
					break ;
			}
			firstRegion.b = subexonIdx[j - 1] ;
			lastRegion.b = subexonIdx[ seIdxCnt - 1 ] ;
			for ( j = seIdxCnt - 2 ; j >= 0 ; --j )
			{
				if ( subexons[ subexonIdx[j] ].end < subexons[ subexonIdx[j + 1] ].start - 1 )
					break ;
			}
			lastRegion.a = subexonIdx[j + 1] ;

			for ( j = 0 ; j < seIdxCnt ; ++j )
				transcriptLength[i] += subexons[ subexonIdx[j] ].end - subexons[ subexonIdx[j] ].start + 1 ;
			
			//for ( j = firstRegion.b ; j < lastRegion.a ; ++j )
			for ( j = 0 ; j < seIdxCnt - 1 ; ++j )
			{
				chain[j].a = subexonIdx[j] ;
				chain[j].b = subexonIdx[j + 1] ;
				covered[j] = false ;
			}
			memset( usedConstraints, false, sizeof( bool ) * constraints.constraints.size() ) ;
			int compatibleCnt = 0 ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				alltranscripts[i].constraintsSupport[j] = 0 ;
				if ( btable[i].Test(j) && tc[j].abundance > 0 )
				{	
					++compatibleCnt ;
					double adjustAbundance = tc[j].abundance ;
					if ( seIdxCnt > 1 )
					{
						if ( tc[j].i == tc[j].j 
								&& ( constraints.constraints[ tc[j].i ].first + 
									constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].first 
									||  constraints.constraints[ tc[j].i ].first + 
									constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].last ) )
						{
							adjustAbundance = inf ;
						}
						else if ( tc[j].i != tc[j].j 
								&& ( constraints.constraints[ tc[j].i ].first + 
									constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].first 
									||  constraints.constraints[ tc[j].i ].first + 
									constraints.constraints[ tc[j].i ].last == 2 * alltranscripts[i].last ) )
						{
							adjustAbundance = constraints.constraints[ tc[j].j ].normAbund ;
						}
						else if ( tc[j].i != tc[j].j 
								&& ( constraints.constraints[ tc[j].j ].first + 
									constraints.constraints[ tc[j].j ].last == 2 * alltranscripts[i].first 
									||  constraints.constraints[ tc[j].j ].first + 
									constraints.constraints[ tc[j].j ].last == 2 * alltranscripts[i].last ) )
						{
							adjustAbundance = constraints.constraints[ tc[j].i ].normAbund ;
						}
					}
					if ( adjustAbundance < value )
						/*!( seIdxCnt > 1 
							&& ( ( ( constraints.constraints[ tc[j].i ].first >= firstRegion.a && constraints.constraints[ tc[j].i ].last <= firstRegion.b ) 
									&& ( constraints.constraints[ tc[j].j ].first >= firstRegion.a && constraints.constraints[ tc[j].j ].last <= firstRegion.b ) )
						          || ( ( constraints.constraints[ tc[j].i ].first >= lastRegion.a && constraints.constraints[ tc[j].i ].last <= lastRegion.b )
									&& ( constraints.constraints[ tc[j].j ].first >= lastRegion.a && constraints.constraints[ tc[j].j ].last <= lastRegion.b ) ) ) )
					   )*/
					{
						// Not use the constraints totally within the 3'/5'-end in the transcript
						value = adjustAbundance ;
						tag = j ;
					}
					avgTranscriptAbundance[i] += tc[j].abundance ;

					if ( !usedConstraints[ tc[j].i ] )
					{
						struct _constraint &c = constraints.constraints[ tc[j].i ] ;
						for ( k = 0 ; k < seIdxCnt - 1 ; ++k )
						{
							// Note that since the constraint is already compatible with the txpt,
							//   chain[k].a/b must be also adjacent in this constraint.
							if ( c.vector.Test( chain[k].a ) && c.vector.Test( chain[k].b ) ) 
								covered[k] = true ;
						}
						usedConstraints[ tc[j].i ] = true ;
					}

					if ( !usedConstraints[ tc[j].j ] )
					{
						struct _constraint &c = constraints.constraints[ tc[j].j ] ;
						for ( k = 0 ; k < seIdxCnt - 1 ; ++k )
						{
							if ( c.vector.Test( chain[k].a ) && c.vector.Test( chain[k].b ) )
								covered[k] = true ;
						}
						usedConstraints[ tc[j].j ] = true ;
					}
				}
			}

			// Get some penalty if something should together did not show up together
			int size = togetherChain.size() ;
			if ( size > 0 )
			{
				BitTable bufferTable( seCnt ) ;
				for ( j = 0 ; j < size ; ++j )
				{
					bufferTable.Assign( togetherChain[j] ) ;
					bufferTable.And( alltranscripts[i].seVector ) ;
					//if ( !bufferTable.IsAllZero() && !bufferTable.IsEqual( togetherChain[j] ) ) 
					//	value /= 2 ;

					if ( !bufferTable.IsAllZero() )
					{
						if ( bufferTable.IsEqual( togetherChain[j] ) ) 
							//printf( "nice together!\n" ) ;
							;
						else
							value /= 2 ;
							//printf( "bad together!\n" ) ;
					}
				}
				bufferTable.Release() ;
			}


			// Every two-subexon chain should be covered by some reads if a transcript is expressed highly enough
			int cnt = 0 ;
			for ( j = 0 ; j < seIdxCnt - 1 ; ++j )	
				if ( covered[j] == false ) // && j >= firstRegion.b && j <= lastRegion.a - 1 )
				{
					value = 0 ;
				}
				else
					++cnt ;
			if ( seIdxCnt > 1 )
				coveredPortion[i] = (double)cnt / (double)( seIdxCnt - 1 ) ;
			else
				coveredPortion[i] = 1 ;
			if ( coveredPortion[i] == 0 )
				coveredPortion[i] = (double)0.5 / ( seIdxCnt ) ;

			// For short subexon (readLength-subexon_length-1>30), we further require a constraint cover three conseuctive subexon
			/*memset( usedConstraints, false, sizeof( bool ) * constraints.constraints.size() ) ;
			for ( j = 1 ; j < seIdxCnt - 1 ; ++j )
			{
				int k = subexonIdx[j] ;	
				if ( alignments.readLen - ( subexons[k].end - subexons[k].start + 1 ) - 1 <= 30 )
					continue ;
				// We need at least one of the side subexons are adjacent to the center one.
				if ( subexons[ subexonIdx[j - 1] ].end + 1 < subexons[k].start && subexons[k].end + 1 < subexons[ subexonIdx[j + 1] ].start )
					continue ;
				
				int l = 0 ; 
				for ( l = 0 ; l < tcCnt ; ++l )
				{
					if ( btable[i].Test(l) && tc[l].abundance > 0 )
					{
						if ( !usedConstraints[ tc[l].i ] )
						{
							struct _constraint &c = constraints.constraints[ tc[l].i ] ;
							if ( c.vector.Test( subexonIdx[j - 1] ) && c.vector.Test( subexonIdx[j] ) &&
								c.vector.Test( subexonIdx[j + 1] ) ) 
								break ;
							usedConstraints[ tc[l].i ] = true ;
						}

						if ( !usedConstraints[ tc[l].j ] )
						{
							struct _constraint &c = constraints.constraints[ tc[l].j ] ;
							if ( c.vector.Test( subexonIdx[j - 1] ) && c.vector.Test( subexonIdx[j] ) &&
								c.vector.Test( subexonIdx[j + 1] ) ) 
								break ;
							usedConstraints[ tc[l].j ] = true ;
						}
					}	
				}
				// It is not covered
				if ( l >= tcCnt )
				{
					int residual = alignments.readLen - ( subexons[k].end - subexons[k].start + 1 ) - 1 ;
					//printf( "residual: %d %d %lf\n", k, residual, value ) ;
					if ( value * residual > 2 )
					{
						value = 1 / (double)residual ;
					}
				}
			}*/

			if ( tag == -1 ) 
				value = 0 ;
			if ( value > maxAbundance )
				maxAbundance = value ;
			transcriptAbundance[i] = value ;
			if ( tag != -1 )
				avgTranscriptAbundance[i] /= compatibleCnt ;

			//printf( "abundance %d: %lf %lf ", i, value, avgTranscriptAbundance[i] ) ;
			//alltranscripts[i].seVector.Print() ;
		}
		if ( maxAbundance == 0 )
		{
			for ( i = 0 ; i < atcnt ; ++i )
			{
				transcriptAbundance[i] = coveredPortion[i] ;
			}
			maxAbundance = 1 ;
		}
		//printf( "%s: %lf\n", __func__, maxAbundance ) ;
		int size = togetherChain.size() ;
		for ( j = 0 ; j < size ; ++j )
			togetherChain[j].Release() ;
		delete[] usedConstraints ;
		delete[] covered ;
		delete[] chain ;
	}
	else 
	{
		for ( i = 0 ; i < atcnt ; ++i )
		{
			transcriptAbundance[i] = alltranscripts[i].abundance ;
			if ( transcriptAbundance[i] > maxAbundance )
				maxAbundance = transcriptAbundance[i] ;
			coveredPortion[i] = 1 ;
		}
		if ( maxAbundance == 0 )
			maxAbundance = 1 ;
	}

	// Obtain the prefix, suffix information of the transcripts.
	int *nextSuffix, *nextPrefix ;
	struct _pair32 *txptRank ;
	nextSuffix = new int[atcnt] ;
	nextPrefix = new int[atcnt] ;
	txptRank = new struct _pair32[atcnt] ;
	memset( nextSuffix, -1, sizeof( int ) * atcnt ) ;
	memset( nextPrefix, -1, sizeof( int ) * atcnt ) ;
	/*for ( i = 0 ; i < atcnt ; ++i )
	{
		std::vector<int> subexonIdx ;
		txptRank[i].a = i ;
		alltranscripts[i].seVector.GetOnesIndices( subexonIdx ) ;
		txptRank[i].b = subexonIdx.size() ; 
	}
	qsort( txptRank, atcnt, sizeof( struct _pair32 ), CompPairsByB) ;
	BitTable bufferTable( seCnt ) ;
	for ( i = atcnt - 1 ; i >= 0 ; --i )
	{
		int a = txptRank[i].a ;
		for ( j = i - 1 ; j >= 0 ; --j )
		{
			if ( txptRank[i].b == txptRank[j].b )
				continue ;

			int b = txptRank[j].a ;

			if ( alltranscripts[b].last != alltranscripts[a].last )
				continue ;

			bufferTable.Assign( alltranscripts[a].seVector ) ; 
			bufferTable.MaskRegionOutside( alltranscripts[b].first, alltranscripts[b].last ) ;
			if ( bufferTable.IsEqual( alltranscripts[b].seVector ) )
			{
				nextSuffix[a] = b ;
				break ;
			}
		}
	}
	for ( i = atcnt - 1 ; i >= 0 ; --i )
	{
		int a = txptRank[i].a ;
		for ( j = i - 1 ; j >= 0 ; --j )
		{
			if ( txptRank[i].b == txptRank[j].b )
				continue ;

			int b = txptRank[j].a ;

			if ( alltranscripts[b].first != alltranscripts[a].first )
				continue ;

			bufferTable.Assign( alltranscripts[a].seVector ) ; 
			bufferTable.MaskRegionOutside( alltranscripts[b].first, alltranscripts[b].last ) ;
			if ( bufferTable.IsEqual( alltranscripts[b].seVector ) )
			{
				nextPrefix[a] = b ;
				break ;
			}
		}
	}

	bufferTable.Release() ;*/
	delete[] txptRank ;

	// Quantative Set-Cover
	int iterCnt = -1 ;
	double *coverCnt = new double[atcnt] ;
	for ( i = 0 ; i < atcnt ; ++i )
		coverCnt[i] = -1 ;
	int *list = new int[atcnt] ;
	int listCnt ;

	while ( 1 )
	{
		double max = -1 ;
		int maxtag = -1 ;
		double maxcnt = -1 ;
		++iterCnt ;

		// Find the optimal candidate.
		for ( i = 0 ; i < atcnt ; ++i )
		{
			double value = inf ;
			double cnt = 0 ;

			if ( coverCnt[i] == -1 )
			{
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( tc[j].abundance > 0 && btable[i].Test( j ) )
					{
						cnt += tc[j].effectiveCount ; 
					}
				}
				/*else
				  {
				  std::vector<int> tcIdx ;
				  btable[i].GetOnesIndices( tcIdx ) ;
				  int size = tcIdx.size() ;
				  for ( j = 0 ; j < size ; ++j )
				  {
				  if ( tc[ tcIdx[j] ].abundance > 0 )	
				  {
				  cnt += tc[ tcIdx[j] ].effectiveCount ;
				  }
				  }
				  }*/
				coverCnt[i] = cnt ;
			}
			else
			{
				cnt = coverCnt[i] ;
			}

			value = transcriptAbundance[i] ;
			if ( cnt < 1 ) // This transcript does not satisfy any undepleted constraints.
				continue ;
			cnt *= coveredPortion[i] ;

			double weight = 1 ; //* seCnt / transcriptSeCnt[i] ;
			//if ( maxAbundance >= 1 && value / maxAbundance >= 0.2 )
			//	seCntAdjust = sqrt( (double)( transcriptSeCnt[i] ) / seCnt ) ;//< 0.5 ? 0.5 : (double)( transcriptSeCnt[i] ) / seCnt ;
			
			if ( alltranscripts[i].FPKM > 0 && sampleCnt > 1 )
				weight = ( 1 + alltranscripts[i].FPKM / sampleCnt ) ;

			double score = ComputeScore( cnt, weight, value, maxAbundance, alltranscripts[i].correlationScore ) ;
			if ( cnt > maxcnt )
				maxcnt = cnt ;
			score += adjustScore[i] ;
			if ( score > max )
			{
				max = score ;
				maxtag = i ;
			}
			else if ( score == max )
			{
				if ( avgTranscriptAbundance[maxtag] < avgTranscriptAbundance[i] )	
				{
					max = score ;
					maxtag = i ;
				}
			}
			//printf( "score: %d %lf -> %lf\n", i, cnt, score ) ;
		}

		if ( maxcnt == 0 || maxtag == -1 )
			break ;

		// Find the constraint that should be depleted.
		double update = inf ;
		int updateTag = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[ maxtag ].Test( j ) && tc[j].abundance > 0 && 
					tc[j].abundance <= update )
			{
				update = tc[j].abundance ;	
				updateTag = j ;
			}
		}

		// Search suffix and prefix to see whether these fit better.
		int p = nextSuffix[ maxtag] ;
		while ( p != -1 )
		{
			if ( transcriptAbundance[p] >= 10.0 * transcriptAbundance[maxtag] 
				&& btable[p].Test( updateTag ) )
			{
				//printf( "%d\n", p ) ;
				maxtag = p ;
				break ;
			}
			p = nextSuffix[p] ;
		}
		p = nextPrefix[maxtag] ;
		while ( p != -1 )
		{
			if ( transcriptAbundance[p] >= 10.0 * transcriptAbundance[maxtag] 
				&& btable[p].Test( updateTag ) )
			{
				maxtag = p ;
				break ;
			}
			p = nextPrefix[p] ;
		}

		
		// Update the abundance.
		int supportCnt = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( btable[maxtag].Test( j ) )
			{
				if ( tc[j].abundance > 0 )
				{
					tc[j].abundance -= 1 * update ;
					double factor = tc[j].effectiveCount ;
					double tmp = ( tc[j].support * update / tc[j].normAbund * factor ) ;
					alltranscripts[maxtag].constraintsSupport[j] += tmp ;
					alltranscripts[maxtag].abundance += tmp ;

					if ( tc[j].abundance <= 0 )
					{
						int l ;
						for ( l = 0 ; l < atcnt ; ++l )
						{
							if ( btable[l].Test(j) )
								coverCnt[l] -= tc[j].effectiveCount ;
						}
					}
					++supportCnt ;
				}
				else if ( alltranscripts[maxtag].constraintsSupport[j] == 0 ) 
				{
					double sum = 0 ;
					double takeOut = 0 ;
					double factor = tc[j].effectiveCount ;
					listCnt = 0 ;
					for ( i = 0 ; i < atcnt ; ++i )
					{
						if ( i == maxtag )
							continue ;

						if ( alltranscripts[i].abundance > 0 && btable[i].Test(j) )
						{
							sum += alltranscripts[i].constraintsSupport[j] ;

							double tmp =  ( alltranscripts[i].constraintsSupport[j] + alltranscripts[maxtag].constraintsSupport[j] ) * 
								transcriptAbundance[maxtag] / ( transcriptAbundance[maxtag] + transcriptAbundance[i] ) 
								- alltranscripts[maxtag].constraintsSupport[j] ;
							if ( tmp > 0 )
							{
								list[ listCnt ] = i ;
								++listCnt ;
								takeOut += tmp ; //alltranscripts[i].constraintsSupport[j] * transcriptAbundance[maxtag] / ( transcriptAbundance[maxtag] + transcriptAbundance[i] ) ;
							}
						}
					}

					double ratio = 1 ;
					double takeOutFactor = 0.5 ;
					if ( update < tc[j].normAbund )
					{
						if ( takeOut > ( tc[j].support * update / tc[j].normAbund * factor ) * takeOutFactor )
							ratio = ( tc[j].support * update / tc[j].normAbund * factor  ) * takeOutFactor / takeOut ;
					}
					else
					{
						if ( takeOut > ( tc[j].support * factor ) * takeOutFactor )
							ratio = ( tc[j].support * factor ) * takeOutFactor / takeOut ;
					}

					if ( 1 ) //update < tc[j].normAbund )
					{
						for ( i = 0 ; i < listCnt ; ++i )
						{
							//double tmp = ( tc[j].support * update / tc[j].normAbund * factor ) * 
							//	( alltranscripts[ list[i] ].constraintsSupport[j] / sum  ) ;
							//if ( alltranscripts[ list[i] ].constraintsSupport[j] < tmp )
							//	printf( "WARNING! %lf %lf, %lf\n", alltranscripts[ list[i] ].constraintsSupport[j], sum, tmp ) ;

							//double tmp = alltranscripts[ list[i] ].constraintsSupport[j] * transcriptAbundance[maxtag] / ( transcriptAbundance[maxtag] + transcriptAbundance[ list[i] ] ) * ratio ; 
							double tmp =  ( ( alltranscripts[ list[i] ].constraintsSupport[j] + alltranscripts[maxtag].constraintsSupport[j] ) * 
									transcriptAbundance[maxtag] / ( transcriptAbundance[maxtag] + transcriptAbundance[ list[i] ] ) 
									- alltranscripts[maxtag].constraintsSupport[j] ) * ratio ;


							alltranscripts[ list[i] ].constraintsSupport[j] -= tmp ;
							alltranscripts[ list[i] ].abundance -= tmp ;
						}
						//double tmp = ( tc[j].support * update / tc[j].normAbund * factor ) ;
						//printf( "%lf %lf. %lf %lf\n", takeOut, ratio, update, tc[j].normAbund ) ;
						double tmp = takeOut * ratio ;
						alltranscripts[maxtag].constraintsSupport[j] += tmp ;
						alltranscripts[maxtag].abundance += tmp ;
					}
					/*else
					{
					  	double tmp = ( tc[j].support / (double)( listCnt + 1 ) ) * factor ;
						for ( i = 0 ; i < listCnt ; ++i )	
						{
							alltranscripts[ list[i] ].abundance -= alltranscripts[ list[i] ].constraintsSupport[j] ;
	
							alltranscripts[ list[i] ].constraintsSupport[j] = tmp ;
							alltranscripts[ list[i] ].abundance += tmp ;
						}
						alltranscripts[maxtag].constraintsSupport[j] += tmp ;
						alltranscripts[maxtag].abundance += tmp ;
					}*/

				}
			}

			if ( tc[j].abundance < 0 )
			{
				tc[j].abundance = 0 ;

			}
		}
		tc[ updateTag ].abundance = 0 ;
		if ( supportCnt == 0 )
			break ;
		//adjustScore[maxtag] += 1 / (double)tcCnt ;
		//printf( "maxtag=%d %lf %d\n", maxtag, update, updateTag ) ;
	}

	for ( i = 0 ; i < atcnt ; ++i )
	{
		if ( alltranscripts[i].abundance > 0 )
		{
			struct _transcript nt = alltranscripts[i] ;
			nt.seVector.Nullify() ;
			nt.seVector.Duplicate( alltranscripts[i].seVector ) ;
			nt.constraintsSupport = NULL ;
			if ( transcriptAbundance[i] == 0 )
				nt.correlationScore = -1 ;
			else
				nt.correlationScore = 0 ;
			nt.id = i ;
			transcripts.push_back( nt ) ;
		}
	}
	
	// Release the memory of btable.
	for ( i = 0 ; i < atcnt ; ++i )
	{
		delete[] alltranscripts[i].constraintsSupport ;
		btable[i].Release() ;
	}
	delete[] btable ;
	
	delete[] list ;
	delete[] transcriptSeCnt ;
	delete[] transcriptLength ;
	delete[] transcriptAbundance ;
	delete[] avgTranscriptAbundance ;
	delete[] coveredPortion ;
	delete[] adjustScore ;
	delete[] coverCnt ;
	
	delete[] nextPrefix ;
	delete[] nextSuffix ;

	// Redistribute weight if there is some constraints that are unbalanced.
	/*tcnt = transcripts.size() ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int maxRatio = -1 ;
		for ( j = 0 ; j < tcCnt ; ++j )
			if ( transcripts[i].constraintsSupport[j] > 0 )
			{
				double factor = tc[j].effectiveCount ;
				if ( transcripts[])
			}
	}*/
}

void TranscriptDecider::AbundanceEstimation( struct _subexon *subexons, int seCnt, Constraints &constraints, std::vector<struct _transcript> &transcripts )
{
	int tcnt = transcripts.size() ;
	int size ;
	int i, j ;
	if ( tcnt <= 0 )
		return ;
	
	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	int tcCnt = tc.size() ; // transcript constraints

	BitTable *btable = new BitTable[ tcnt ] ; 
	int *transcriptLength = new int[tcnt] ;
	int *compatibleList = new int[tcnt] ;
	double *rho = new double[tcnt] ; // the abundance.
	int iterCnt = 0 ;

	for ( i = 0 ; i < tcnt ; ++i )
		transcripts[i].constraintsSupport = new double[ tcCnt ] ;
	
	for ( i = 0 ; i < tcnt ; ++i )
	{
		btable[i].Init( tcCnt ) ;
		double min = -1 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int a = tc[j].i ;
			int b = tc[j].j ;

			if ( IsConstraintInTranscript( transcripts[i], constraints.constraints[a] ) == 1 
					&& IsConstraintInTranscript( transcripts[i], constraints.constraints[b] ) == 1 )
			{
				//printf( "set btble[ %d ].Set( %d ): %d %d\n", i, j, a, b ) ;
				btable[i].Set( j ) ;

				if ( min == -1 || tc[j].normAbund < min )
					min = tc[j].normAbund ;
			}
		}

		std::vector<int> subexonIdx ;
		transcripts[i].seVector.GetOnesIndices( subexonIdx ) ;
		int subexonIdxCnt = subexonIdx.size() ;
		int len = 0 ;
		for ( j = 0 ; j < subexonIdxCnt ; ++j )
			len += subexons[ subexonIdx[j] ].end - subexons[ subexonIdx[j] ].start + 1 ;
		transcriptLength[i] = len - alignments.fragLen + 2 * alignments.fragStdev ;
		if ( transcriptLength[i] < 1 )
			transcriptLength[i] = 1 ;
		rho[i] = transcripts[i].abundance / transcriptLength[i]  ; // use the rough estimation generated before.
		if ( transcripts[i].correlationScore == -1 && rho[i] > 0.1 / (double)alignments.readLen )
			rho[i] = 0.1 / (double)alignments.readLen  ;
	}
	
	while ( 1 )
	{
		for ( i = 0 ; i < tcnt ; ++i )
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				transcripts[i].constraintsSupport[j] = 0 ;
			}
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			int clCnt = 0 ;
			double sum = 0 ;
			for ( i = 0 ; i < tcnt ; ++i )
			{
				if ( btable[i].Test(j) )
				{
					compatibleList[ clCnt ] = i ;
					++clCnt ;
					sum += rho[i] ;
				}
			}

			for ( i = 0 ; i < clCnt ; ++i )
			{
				double factor = tc[j].effectiveCount ;
				transcripts[ compatibleList[i] ].constraintsSupport[j] = ( rho[ compatibleList[i] ] / sum ) * tc[j].support * factor ;
			}
		}

		double diff = 0 ;
		for ( i = 0 ; i < tcnt ; ++i )
		{
			double newAbund = 0 ;
			for ( j = 0 ; j < tcCnt ; ++j )
				newAbund += transcripts[i].constraintsSupport[j] ;
			double old = rho[i] ;
			rho[i] = newAbund / transcriptLength[i] ;
			//printf( "rho[%d]=%lf\n", i, rho[i] ) ;
			if ( transcripts[i].correlationScore == -1 && rho[i] > 0.1 / (double)alignments.readLen )
				rho[i] = 0.1 / (double)alignments.readLen  ;
			
			double tmp = ( old - rho[i] ) ;
			diff += tmp < 0 ? -tmp : tmp ;
		}
		//printf( "%lf\n", diff ) ;
		if ( diff < 1e-3)
			break ;

		++iterCnt ;
		if ( iterCnt >= 1000 )
			break ;
	}

	for ( i = 0 ; i < tcnt ; ++i )
	{
		//printf( "%lf=>", transcripts[i].abundance ) ;
		transcripts[i].abundance = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			transcripts[i].abundance += transcripts[i].constraintsSupport[j] ;
		}
		//printf( "%lf. (%lf)\n", transcripts[i].abundance, transcripts[i].correlationScore ) ;
		//transcripts[i].seVector.Print() ;
	}

	for ( i = 0 ; i < tcnt ; ++i )
		delete[] transcripts[i].constraintsSupport ;

	// Release the memory of btable.
	for ( i = 0 ; i < tcnt ; ++i )
	{
		btable[i].Release() ;
	}
	delete[] compatibleList ;
	delete[] btable ;
	delete[] transcriptLength ;
	delete[] rho ;
}

int TranscriptDecider::RefineTranscripts( struct _subexon *subexons, int seCnt, bool aggressive,
	std::map<int, int> *subexonChainSupport, int *txptSampleSupport, std::vector<struct _transcript> &transcripts, Constraints &constraints ) 
{
	int i, j, k ;
	int tcnt = transcripts.size() ;
	if ( tcnt == 0 )
		return 0 ;
	int tcCnt = constraints.matePairs.size() ;

	std::vector<struct _matePairConstraint> &tc = constraints.matePairs ;
	std::vector<struct _constraint> &scc = constraints.constraints ; //single-end constraints.constraints
	
	// Remove transcripts whose FPKM are too small.
	//printf( "%d %d\n", usedGeneId, baseGeneId ) ;
	double *geneMaxFPKM = new double[usedGeneId - baseGeneId ] ;
	int *geneMaxFPKMTag = new int[usedGeneId - baseGeneId ] ;
	double *nonOverlapMaxFPKM = new double[ usedGeneId - baseGeneId ] ; // the max FPKM among all the transcripts not overlapping with maxFPKMTag transcripts.
	memset( geneMaxFPKM, 0, sizeof( double ) * ( usedGeneId - baseGeneId ) ) ;
	memset( geneMaxFPKMTag, 0, sizeof( int ) * ( usedGeneId - baseGeneId ) ) ;
	memset( nonOverlapMaxFPKM, 0, sizeof( double ) * ( usedGeneId - baseGeneId ) ) ;
	
	double *geneMaxCov = new double[ usedGeneId - baseGeneId ] ;
	memset( geneMaxCov, 0, sizeof( double ) * ( usedGeneId - baseGeneId ) ) ;
	int *txptGid = new int[tcnt] ;

	/*for ( i = 0 ; i < tcnt ; ++i )
	{
		printf( "%d: %lf ", i, transcripts[i].FPKM ) ;
		transcripts[i].seVector.Print() ;
	}*/

	/*==================================================================
	  Remove transcripts that has too few relative FPKM. (-f)
	  ====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int gid = GetTranscriptGeneId( transcripts[i], subexons ) ;
		int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
		//printf( "gid=%d\n", gid ) ;
		//printf( "%lf %lf %d\n", transcripts[i].abundance, transcripts[i].FPKM, len ) ;
		if ( transcripts[i].FPKM > geneMaxFPKM[gid - baseGeneId ] )
		{
			geneMaxFPKM[ gid - baseGeneId ] = transcripts[i].FPKM ;
			geneMaxFPKMTag[ gid - baseGeneId ] = i ;
		}
		if ( transcripts[i].abundance * alignments.readLen / len > geneMaxCov[gid - baseGeneId ] )
			geneMaxCov[gid - baseGeneId] = ( transcripts[i].abundance * alignments.readLen ) / len ;
		txptGid[i] = gid ;
	}

	for ( i = 0 ; i < tcnt ; ++i )
	{
		int tag = txptGid[i] - baseGeneId ;
		if ( ( transcripts[i].last < transcripts[ geneMaxFPKMTag[ tag ] ].first 
			|| transcripts[i].first > transcripts[ geneMaxFPKMTag[tag] ].last ) && transcripts[i].FPKM > nonOverlapMaxFPKM[tag] )
			nonOverlapMaxFPKM[tag] = transcripts[i].FPKM ;
	}
	BitTable bufferTable ;
	bufferTable.Duplicate( transcripts[0].seVector ) ;

	if ( !aggressive )
	{
		// Rescue the transcripts covering unique constraints.
		int cnt = 0 ;
		int tag = 0 ;
		int *uniqCount = new int[tcnt] ;
		memset( uniqCount, 0, sizeof( int ) * tcnt ) ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			cnt = 0 ;
			if ( tc[j].uniqSupport <= 5 )
				continue ;
			for ( i = 0 ; i < tcnt ; ++i )
			{
				if ( IsConstraintInTranscript( transcripts[i], scc[ tc[j].i ] ) && 
					IsConstraintInTranscript( transcripts[i], scc[ tc[j].j] ) )
				{
					tag = i ;
					++cnt ;
				}
				if ( cnt >= 2 )
					break ;
			}
			if ( cnt == 1 )
			{
				++uniqCount[tag] ;
			}
		}
		for ( i = 0 ; i < tcnt ; ++i )
		{
			if ( uniqCount[i] >= 2 )
			{
				transcripts[i].abundance *= 4 ;
				transcripts[i].FPKM *= 4 ;
			}
		}

		delete[] uniqCount ;
	}

	int sccCnt = scc.size() ;
	double filterFactor = 1.0 ;
	
	for ( i = 0 ; i < tcnt ; ++i )
	{
		//printf( "%d: %lf %lf\n", txptGid[i], transcripts[i].abundance, geneMaxFPKM[ txptGid[i] - baseGeneId ] ) ;

		if ( transcripts[i].FPKM < filterFactor * FPKMFraction * geneMaxFPKM[ txptGid[i] - baseGeneId ] )
		{
			/*int cnt = 0 ;
			int coverCnt = 0 ;
			for ( j = 0 ; j < tcCnt ; ++j )
			{
				if ( transcripts[i].constraintsSupport[j] > 0 )
					++coverCnt ;
				double factor = tc[j].effectiveCount ;
				if ( transcripts[i].constraintsSupport[j] >= factor * tc[j].support - 1e-3
					&& tc[j].support >= 10 
					&& tc[j].uniqSupport >= 0.95 * tc[j].support ) 
				{
					++cnt ;
				}
			}
			//cnt = 0 ;	
			if ( cnt >= 2 )
			{
				;
			}
			else*/
				transcripts[i].abundance = -transcripts[i].abundance ;
		}
		//if ( transcripts[i].FPKM >= 0.8 * geneMaxFPKM[ txptGid[i] - baseGeneId ] && geneMaxCov[ txptGid[i] - baseGeneId ] >= txptMinReadDepth )
		//	continue ;
	}

	if ( nonOverlapMaxFPKM != 0 )
	{
		// Go two iterations to rescue, the first iteration should be just for marking.
		std::vector<int> rescueList ;
		for ( i = 0 ; i < tcnt ; ++i )
		{
			if ( transcripts[i].abundance >= 0 )
				continue ;

			for ( j = 0 ; j < tcnt ; ++j )
			{
				if ( transcripts[j].abundance < 0 || txptGid[i] != txptGid[j] )
					continue ;
				if ( transcripts[i].first <= transcripts[j].last && transcripts[i].last >= transcripts[j].first )
				/*bufferTable.Assign( transcripts[i].seVector ) ;
				bufferTable.And( transcripts[j].seVector ) ;
				
				if ( !bufferTable.IsAllZero() )*/
					break ;
			}
			if ( j >= tcnt && transcripts[i].FPKM >= FPKMFraction * nonOverlapMaxFPKM[ txptGid[i] - baseGeneId ] )
			{
				//transcripts[i].abundance = -transcripts[i].abundance ;	
				rescueList.push_back( i ) ;
			}
		}
		
		int size = rescueList.size() ;
		for ( i = 0 ; i < size ; ++i )
			transcripts[ rescueList[i] ].abundance *= -1 ;
	}

	/*==================================================================
	  Remove transcripts that has too few read coverage (-d)
	  ====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].abundance >= 0 )
		{
			int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
			double cov = ( transcripts[i].abundance * alignments.readLen ) / len ;
			//printf( "%d: %d %d %lf %lf\n", i, len, transcripts[i].seVector.Count(), cov, geneMaxCov[ txptGid[i] - baseGeneId ]  ) ;
			
			if ( ( tcnt > 1 || len <= 1000 || transcripts[i].seVector.Count() <= 3 ) && cov < txptMinReadDepth )
			{
				//if ( usedGeneId == baseGeneId + 1 && /*transcripts[i].seVector.Count() > 3 
				//	&& len > 1000  &&*/ geneMaxCov[ txptGid[i] - baseGeneId ] == cov ) 
				if ( geneMaxCov[ txptGid[i] - baseGeneId ] == cov ) 
					continue ;

				// Test whether it has some very abundant constraints.
				/*int cnt = 0 ;
				for ( j = 0 ; j < tcCnt ; ++j )
				{
					if ( transcripts[i].constraintsSupport[j] >= tc[j].support / 2.0 
							&& tc[j].support >= 10 
							&& tc[j].uniqSupport >= 0.95 * tc[j].support 
							&& tc[j].normAbund >= 1  )	
					{
						++cnt ;
					}
				}

				if ( cnt >= 1 )
				{
					continue ;
				}*/

				// Test whether this transcript is fully covered. If so ,we can filter it.
				
				if ( geneMaxCov[ txptGid[i] - baseGeneId ] <= 5 )
				{
					bufferTable.Reset() ;
					for ( j = 0 ; j < sccCnt ; ++j )	
					{
						if ( !IsConstraintInTranscript( transcripts[i], scc[j] ) ) 
							continue ;
						bufferTable.Or( scc[j].vector ) ;
					}
					if ( bufferTable.IsEqual( transcripts[i].seVector ) )
						transcripts[i].abundance = -transcripts[i].abundance ;
				}
				else
					transcripts[i].abundance = -transcripts[i].abundance ;
					
				/*else
				  {
				  transcripts[i].seVector.Print() ;
				  bufferTable.Print() ;
				  OutputTranscript( stderr, subexons, transcripts[i] ) ;
				  }*/
			}
		}
	}
	
	/*==================================================================
	  Remove transcripts that is too short
	  ====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].abundance <= 0 )
			continue ;

		int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[i].abundance, transcripts[i].FPKM ) ;
		if ( len < 200 )
		{
			transcripts[i].abundance = -transcripts[i].abundance ;
		}
	}

	// Rescue transcripts that showed up in many samples.
	/*for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( transcripts[i].abundance > 0 )
			continue ;
		if ( txptSampleSupport[ transcripts[i].id ] >= 3 &&
			txptSampleSupport[transcripts[i].id ] >= (int)( sampleCnt / 2 ) ) 
			transcripts[i].abundance = -transcripts[i].abundance ;
	}*/

	// Rescue some transcripts covering subexon chains showed up in many samples, but missing after filtration.
	struct _constraint tmpC ;
	tmpC.vector.Init( seCnt ) ;
	
	std::vector< struct _pair32 > missingChain ;
	std::vector<int> recoverCandidate ;
	bool *used = new bool[tcnt] ;
	memset( used, false, sizeof( bool ) * tcnt ) ;

	// Obtain the list of transcripts that should be recovered.
	for ( i = 0 ; i < seCnt && sampleCnt > 1  ; ++i )
	{

		double maxFPKM = -1 ;
		for ( std::map<int, int>::iterator it = subexonChainSupport[i].begin() ; 
			it != subexonChainSupport[i].end() ; ++it )
		{
			if ( sampleCnt >= 0 && ( it->second < 3 || it->second < (int)( 0.5 * sampleCnt ) ) && it->second <= sampleCnt / 2 )	
				continue ;
			
			bool recover = true ;
			tmpC.vector.Reset() ;
			tmpC.vector.Set( i ) ;
			tmpC.vector.Set( it->first ) ;
			tmpC.first = i ; 
			tmpC.last = it->first ;

	
			for ( j = 0 ; j < tcnt ; ++j )
			{
				if ( transcripts[j].abundance < 0 )
					continue ;

				if ( IsConstraintInTranscript( transcripts[j], tmpC ) ) 
				{
					recover = false ;
					break ;
				}

				if ( recover )
				{
					for ( j = 0 ; j < tcnt ; ++j )
					{
						if ( transcripts[j].abundance > 0 )
							continue ;
						//printf( "%d %lf\n", IsConstraintInTranscript( transcripts[j], tmpC ), transcripts[j].FPKM ) ;
						if ( IsConstraintInTranscript( transcripts[j], tmpC ) )
						{
							/*if ( maxTag == -1 )
								maxTag = j ;
							else 
							{
								if ( txptSampleSupport[ transcripts[j].id ] > txptSampleSupport[ transcripts[maxTag ].id ] )
									maxTag = j ;
								else if ( txptSampleSupport[ transcripts[j].id ] == txptSampleSupport[ transcripts[maxTag ].id ])
								{
									if ( transcripts[j].FPKM > transcripts[maxTag].FPKM )
										maxTag = j ;
								}
							}*/

							struct _pair32 np ;
							np.a = i ; np.b = it->first ;
							missingChain.push_back( np ) ;

							if ( !used[j] )
							{
								recoverCandidate.push_back( j ) ;
								used[j] = true ;
							}
						}
					}

					/*if ( maxTag != -1 && txptSampleSupport[ transcripts[maxTag].id ] > 1 )
					{
						//printf( "recover %d %d\n", maxTag, txptSampleSupport[ transcripts[maxTag].id ] ) ;
						transcripts[maxTag].abundance *= -1 ;
					}*/
				}
			}

		}
	}

	int size = recoverCandidate.size() ;
	memset( used, false, sizeof( bool ) * tcnt ) ;
	// Recover the candidates in the order of reliability
	int *geneRecoverCnt = new int[ usedGeneId - baseGeneId ] ;
	memset( geneRecoverCnt, 0, sizeof( int ) * ( usedGeneId - baseGeneId ) ) ;
	int round = 1 ;
	if ( aggressive && size > 1 )
		round = 1 ;
		
	for ( i = 0 ; i < size ; ++i )
	{
		int maxTag = -1 ;
		int maxCover = -1 ;
		for ( j = 0 ; j < size ; ++j )
		{
			if ( !used[ recoverCandidate[j] ] )
			{
				/*int cover = 0 ; 

				k = missingChain.size() ;
				int l ;
				for ( l = 0 ; l < k ; ++l )
				{
					if ( missingChain[l].a == -1 )
						continue ;

					tmpC.vector.Reset() ;
					tmpC.vector.Set( missingChain[l].a ) ;
					tmpC.vector.Set( missingChain[l].b ) ;
					tmpC.first = missingChain[l].a ; 
					tmpC.last = missingChain[l].b ;

					if ( IsConstraintInTranscript( transcripts[ recoverCandidate[j] ], tmpC ) )
					{
						++cover ;
					}
				}*/

				if ( maxTag == -1 )
				{
					maxTag = recoverCandidate[j] ;
					//maxCover = cover ;
					continue ;
				}
								
				/*if ( cover > maxCover )
				{
					maxTag = recoverCandidate[j] ;
					maxCover = cover ;
				}
				else if ( cover == maxCover )
				{*/
				if ( txptSampleSupport[ transcripts[ recoverCandidate[j] ].id ] > 
						txptSampleSupport[ 
							transcripts[ maxTag ].id 
							] )
					maxTag = recoverCandidate[j] ;
				else if ( txptSampleSupport[ transcripts[ recoverCandidate[j] ].id ] ==
						txptSampleSupport[ transcripts[ maxTag ].id ] )
				{
					if ( transcripts[ recoverCandidate[j] ].FPKM > transcripts[ maxTag ].FPKM ) 
						maxTag = recoverCandidate[j] ;
				}

				/*else if ( transcripts[ recoverCandidate[j] ].FPKM > transcripts[ maxTag ].FPKM ) 
					maxTag = recoverCandidate[j] ;
				else if ( transcripts[ recoverCandidate[j] ].FPKM == transcripts[ maxTag ].FPKM ) 
				{
					if ( txptSampleSupport[ transcripts[ recoverCandidate[j] ].id ] > 
						txptSampleSupport[ transcripts[ maxTag ].id ] )
					maxTag = recoverCandidate[j] ;
				}*/
				//}
			}
		}

		if ( maxTag == -1 || txptSampleSupport[ transcripts[ maxTag ].id ] <= 2 
			|| txptSampleSupport[ transcripts[maxTag].id ] < 0.5 * sampleCnt )
			break ;
		
		used[maxTag] = true ;
		if ( geneRecoverCnt[ txptGid[maxTag] - baseGeneId ] >= round )
			continue ;
		++geneRecoverCnt[ txptGid[maxTag] - baseGeneId ] ;

		k = missingChain.size() ;
		int cnt = 0 ;
		for ( j = 0 ; j < k ; ++j )
		{
			if ( missingChain[j].a == -1 )
				continue ;

			tmpC.vector.Reset() ;
			tmpC.vector.Set( missingChain[j].a ) ;
			tmpC.vector.Set( missingChain[j].b ) ;
			tmpC.first = missingChain[j].a ; 
			tmpC.last = missingChain[j].b ;

			if ( IsConstraintInTranscript( transcripts[maxTag], tmpC ) )
			{
				missingChain[j].a = -1 ; 
				++cnt ;
			}
		}

		int len = GetTranscriptLengthFromAbundanceAndFPKM( transcripts[maxTag].abundance, transcripts[maxTag].FPKM ) ;
		double cov = ( transcripts[maxTag].abundance * alignments.readLen ) / len ;
		if ( cnt >= 1 &&  cov > 1.0 ) 
		{
			transcripts[maxTag].abundance *= -1 ;
		}
	}
	delete[] used ;
	delete[] geneRecoverCnt ;
	tmpC.vector.Release() ;


	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;



	delete []geneMaxCov ;
	bufferTable.Release() ;
	delete []geneMaxFPKM ;
	delete []geneMaxFPKMTag ;
	delete []nonOverlapMaxFPKM ;
	delete []txptGid ;

	/*==================================================================
	  Remove transcripts that seems duplicated
	  ====================================================================*/
	for ( i = 0 ; i < tcnt ; ++i )
	{
		int support = 0 ;
		int uniqSupport = 0 ;

		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( !IsConstraintInTranscript( transcripts[i], scc[ tc[j].i ] ) || !IsConstraintInTranscript( transcripts[i], scc[ tc[j].j ] ) )
				continue ;
			//support += scc[ tc[j].i ].support + scc[ tc[j].j ].support ;
			//uniqSupport += scc[ tc[j].i ].uniqSupport + scc[ tc[j].j ].uniqSupport ; 
			support += tc[j].support ;
			uniqSupport += tc[j].uniqSupport ;

			//printf( "constraint uniqness: %d: %d %d\n", i, tc[j].uniqSupport, tc[j].support ) ;
		}
		//printf( "%d: %d %d\n", i, uniqSupport, support ) ;
		if ( (double)uniqSupport < 0.03 * support )
			transcripts[i].abundance = -1 ;
	}
	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;

	
	/*==================================================================
	  Remove shadow transcripts, the abnormal 2-exon txpt whose intron is very close to the true one or one of the anchor exon is shorter than 25bp....
	  ====================================================================*/
	int minusCnt = 0, plusCnt = 0 ;
	int mainStrand ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].rightStrand == 1 )	
			++plusCnt ;
		else if ( subexons[i].rightStrand == -1 )
			++minusCnt ;
	}
	if ( plusCnt > minusCnt )
		mainStrand = 1 ;
	else
		mainStrand = -1 ;
	
	for ( i = 0 ; i < tcnt ; ++i )
	{ 
		std::vector<int> subexonIdx ;
		transcripts[i].seVector.GetOnesIndices( subexonIdx ) ;
		int size = subexonIdx.size() ;
		int intronCnt = 0 ;
		int anchorIdx = 0 ; // the subexon adjacent to the only intron.
		
		for ( j = 0 ; j < size - 1 ; ++j )
		{
			if ( subexons[ subexonIdx[j] ].end + 1 < subexons[ subexonIdx[j + 1] ].start )	
			{
				++intronCnt ;
				anchorIdx = j ;
			}
		}
		if ( intronCnt != 1 )
			continue ;
		
		int anchorExonLength[2] = {0, 0};
		int tag = 0 ;
		for ( j = 0 ; j < size ; ++j )
		{
			anchorExonLength[tag] += subexons[ subexonIdx[j] ].end - subexons[ subexonIdx[j] ].start + 1 ;
			if ( tag == 0 && subexons[ subexonIdx[j] ].end + 1 < subexons[ subexonIdx[j + 1] ].start )	
				++tag ;
		}
		
		int flag = 0 ;
		if ( subexons[ subexonIdx[anchorIdx] ].rightStrand == mainStrand )
		{
			j = subexonIdx[ anchorIdx ] ;
			if ( subexons[j].end - subexons[j].start + 1 <= 20 || 
					( subexons[j+ 1].start == subexons[j].end + 1 && subexons[j + 1].end - subexons[j + 1].start + 1 <= 20 
					  && subexons[j + 1].rightStrand == mainStrand ) )
				++flag ;
			j = subexonIdx[ anchorIdx + 1 ] ;
			if ( subexons[j].end - subexons[j].start + 1 <= 20 || 
					( subexons[j].start == subexons[j - 1].end + 1 && subexons[j - 1].end - subexons[j - 1].start + 1 <= 20 
					  && subexons[j - 1].leftStrand == mainStrand ) )
				++flag ;
		}

		if ( anchorExonLength[0] <= 25 || anchorExonLength[1] <= 25 )
			flag = 2 ;
		
		// the alignment support the intron must be unique and has enough support.
		int support = 0 ;
		int uniqSupport = 0 ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			if ( !IsConstraintInTranscript( transcripts[i], scc[ tc[j].i ] ) || !IsConstraintInTranscript( transcripts[i], scc[ tc[j].j ] ) )
				continue ;
			if ( ( scc[ tc[j].i ].vector.Test( subexonIdx[ anchorIdx ] ) && scc[ tc[j].i ].vector.Test( subexonIdx[ anchorIdx + 1 ] ) ) 
			 	||  ( scc[ tc[j].j ].vector.Test( subexonIdx[ anchorIdx ] ) && scc[ tc[j].j ].vector.Test( subexonIdx[ anchorIdx + 1 ] ) ) )
			{
				support += tc[j].support ;
				uniqSupport += tc[j].uniqSupport ;
			}

		}

		if ( (double)uniqSupport < 0.3 * support || support < txptMinReadDepth )
		{
			flag = 2 ;
		}

		if ( flag == 2 )
			transcripts[i].abundance = -1 ;		
	
	}
	tcnt = RemoveNegativeAbundTranscripts( transcripts )  ;
	
	return transcripts.size() ;
}

void TranscriptDecider::ComputeTranscriptsScore( struct _subexon *subexons, int seCnt, std::map<int, int> *subexonChainSupport, std::vector<struct _transcript> &transcripts ) 
{
	int i, j ;
	int tcnt = transcripts.size() ;
	struct _constraint tmpC ;
	tmpC.vector.Init( seCnt ) ;

	for ( i = 0 ; i < tcnt ; ++i )
		transcripts[i].correlationScore = 0 ;
	
	for ( i = 0 ; i < seCnt ; ++i )
	{
		for ( std::map<int, int>::iterator it = subexonChainSupport[i].begin() ; 
				it != subexonChainSupport[i].end() ; ++it )
		{
			if ( sampleCnt >= 0 && ( it->second < 3 || it->second < (int)( 0.1 * sampleCnt ) ) && it->second <= sampleCnt / 2 )	
				continue ;

			tmpC.vector.Reset() ;
			tmpC.vector.Set( i ) ;
			tmpC.vector.Set( it->first ) ;
			tmpC.first = i ; 
			tmpC.last = it->first ;
			
			for ( j = 0 ; j < tcnt ; ++j )
			{
				if ( IsConstraintInTranscript( transcripts[j], tmpC ) )
					++transcripts[j].correlationScore ;
			}
		}
	}
	
	tmpC.vector.Release() ;
}

int TranscriptDecider::Solve( struct _subexon *subexons, int seCnt, std::vector<Constraints> &constraints, SubexonCorrelation &subexonCorrelation )
{
	int i, j, k ;
	int cnt = 0 ;
	int *f = new int[seCnt] ; // this is a general buffer for a type of usage.	
	bool useDP = false ;

	compatibleTestVectorT.Init( seCnt ) ; // this is the bittable used in compatible test function.	
	compatibleTestVectorC.Init( seCnt ) ;

	for ( i = 0 ; i < seCnt ; ++i )
	{
		subexons[i].canBeStart = subexons[i].canBeEnd = false ;

		if ( subexons[i].prevCnt == 0 )
			subexons[i].canBeStart = true ;
		else if ( subexons[i].leftClassifier < canBeSoftBoundaryThreshold && subexons[i].leftClassifier != -1 
			&& subexons[i].leftStrand != 0 ) // The case of overhang.
		{
			// We then look into whether there is a left-side end already showed up before this subexon in this region of subexons.
			bool flag = true ;
			for ( j = i - 1 ; j >= 0 ; --j )
			{
				if ( subexons[j].end + 1 != subexons[j + 1].start )	
					break ;
				if ( subexons[i].canBeStart == true )
				{
					flag = false ;
					break ;
				}	
			}
			subexons[i].canBeStart = flag ;
		}

		if ( subexons[i].nextCnt == 0 )
			subexons[i].canBeEnd = true ;
		else if ( subexons[i].rightClassifier < canBeSoftBoundaryThreshold && subexons[i].rightClassifier != -1 
			&& subexons[i].rightStrand != 0 )
		{
			subexons[i].canBeEnd = true ;
		}
		// Remove other soft end already showed up in this region of subexons.
		if ( subexons[i].canBeEnd == true )
		{
			for ( j = i - 1 ; j >= 0 ; --j )
			{
				if ( subexons[j].end + 1 != subexons[j + 1].start )
					break ;
				if ( subexons[j].canBeEnd == true )
				{
					subexons[j].canBeEnd = false ;
					break ;
				}
			}
		}
		//printf( "%d: %d %lf\n", subexons[i].canBeStart, subexons[i].prevCnt, subexons[i].leftClassifier ) ;
	}

	// Go through the cases of mixture region to set canBeStart/End.
	// e.g: +[...]+_____+[....]-...]+____+[..)_____-[...]-
	//                   ^ then we need to force a start point here.
	// Do we need to associate a strand information with canBeStart, canBeEnd?
	for ( i = 0 ; i < seCnt ; )
	{
		// [i, j) is a region.
		for ( j = i + 1 ; j < seCnt ; ++j )
			if ( subexons[j].start > subexons[j - 1].end + 1 )
				break ;
		if ( subexons[i].canBeStart == false ) // then subexons[i] must has a hard left boundary.
		{
			int leftStrandCnt[2] = {0, 0} ;
			for ( k = i ; k < j ; ++k )
			{
				if ( !SubexonGraph::IsSameStrand( subexons[k].rightStrand, subexons[i].leftStrand ) )
					break ;
				if ( subexons[k].leftStrand != 0 )
					++leftStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] ;
			}
			if ( k < j && leftStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] == 0 )
				subexons[i].canBeStart = true ;
		}

		if ( subexons[j - 1].canBeEnd == false )
		{
			int rightStrandCnt[2] = {0, 0} ;
			for ( k = j - 1 ; k >= i ; --k )	
			{
				if ( !SubexonGraph::IsSameStrand( subexons[k].leftStrand, subexons[j - 1].rightStrand ) )
					break ;
				if ( subexons[k].rightStrand != 0 )
					++rightStrandCnt[ ( subexons[k].rightStrand + 1 ) / 2 ] ;
			}
			if ( k >= i && rightStrandCnt[ ( subexons[k].leftStrand + 1 ) / 2 ] == 0 )
				subexons[j - 1].canBeEnd = true ;
		}

		//if ( subexons[i].start == 6870264)
		//	printf( "hi %d %d\n",i , subexons[i].canBeStart ) ;
		i = j ;
	}
	/*for ( i = 0 ; i < seCnt ; ++i )
	{
		printf( "%d %d: %d %d\n", subexons[i].start, subexons[i].end, subexons[i].canBeStart, subexons[i].canBeEnd ) ;
	}*/

	// Find the gene ids.
	baseGeneId = subexons[0].lcCnt ;
	usedGeneId = subexons[0].rcCnt ;
	defaultGeneId[0] = defaultGeneId[1] = -1 ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].geneId < 0 )
			continue ;

		//if ( baseGeneId == -1 || subexons[i].geneId < baseGeneId )	
		//	baseGeneId = subexons[i].geneId ;
		//if ( subexons[i].geneId > usedGeneId )
		//	usedGeneId = subexons[i].geneId ;
	
		if ( ( subexons[i].rightStrand == -1 || subexons[i].leftStrand == -1 ) &&
			( defaultGeneId[0] == -1 || subexons[i].geneId < defaultGeneId[0] ) )
			defaultGeneId[0] = subexons[i].geneId ;
		if ( ( subexons[i].rightStrand == 1 || subexons[i].leftStrand == 1 ) &&
			( defaultGeneId[1] == -1 || subexons[i].geneId < defaultGeneId[1] ) )
			defaultGeneId[1] = subexons[i].geneId ;
	}
	if ( defaultGeneId[0] == -1 )
		defaultGeneId[0] = baseGeneId ;
	if ( defaultGeneId[1] == -1 )
		defaultGeneId[1] = usedGeneId - 1 ;

	// Go through the constraints to find the chain of subexons that should be kept.
	std::map<int, int> *subexonChainSupport = new std::map<int, int>[ seCnt ] ; 
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		std::vector<int> subexonIdx ;
		std::vector<struct _pair32> chain ;
		
		int tcCnt = constraints[i].constraints.size() ;
		int size ;
		for ( j = 0 ; j < tcCnt ; ++j )
		{
			struct _constraint c = constraints[i].constraints[j] ;
			if ( c.uniqSupport < 0.95 * c.support || c.support < 3 )
				continue ;
			
			subexonIdx.clear() ;
			c.vector.GetOnesIndices( subexonIdx ) ;
			size = subexonIdx.size() ;

			for ( k = 0 ; k < size - 1 ; ++k )
			{
				struct _pair32 p ;
				
				p.a = subexonIdx[k] ;
				p.b = subexonIdx[k + 1] ;
				//if ( subexons[p.a].end + 1 == 113235898 && subexons[ p.b ].start + 1 == 113236121 )
				//	printf( "bad bad %d %d %d\n", i, c.uniqSupport, c.support ) ;
				
				if ( subexons[ p.a ].end + 1 < subexons[ p.b ].start )	
					chain.push_back( p ) ;
			}
		}
		// Remove redundancy.
		sort( chain.begin(), chain.end(), CompSortPairs ) ;
		size = chain.size() ;
		k = 0 ;
		for ( j = 1 ; j < size ; ++j )
		{
			if ( chain[j].a == chain[k].a && chain[j].b == chain[k].b )			
				continue ;
			else
			{
				++k ;
				chain[k] = chain[j] ;
			}
		}
		chain.resize( k + 1 ) ;

		// Add those to sample count
		size = k + 1 ;
		for ( j = 0 ; j < size ; ++j )
		{
			if ( subexonChainSupport[ chain[j].a ].count( chain[j].b ) )
			{
				++subexonChainSupport[ chain[j].a ][ chain[j].b ] ;
			}
			else
				subexonChainSupport[ chain[j].a ][ chain[j].b ] = 1 ;
		}
	}

	/*for ( i = 0 ; i < seCnt ; ++i )
	{
		printf( "%d:", i ) ;
		for ( std::map<int, int>::iterator it = subexonChainSupport[i].begin() ; it != subexonChainSupport[i].end() ; ++it )
			printf( " (%d %d) ", it->first, it->second ) ;
		printf( "\n" ) ;
	}*/
	
	//printf( "%d %d %d\n", defaultGeneId[0], baseGeneId, usedGeneId ) ;
	cnt = 0 ;
	memset( f, -1, sizeof( int ) * seCnt ) ;
	for ( i = 0 ; i < seCnt ; ++i )
	{
		if ( subexons[i].canBeStart )	
		{
			cnt += SubTranscriptCount( i, subexons, f ) ;
		}
	}
	if ( cnt <= USE_DP )
	{
		for ( i = 0 ; i < seCnt ; ++i )
			if ( f[i] > USE_DP )
			{
				useDP = true ;
				break ;
			}
	}
	else
		useDP = true ;
	if ( !useDP )
	{
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			double msize = constraints[i].matePairs.size() ;
			double csize = constraints[i].constraints.size() ;
			if ( cnt > ( csize / msize ) * ( csize / msize ) * seCnt 
				&& cnt > USE_DP / ( msize * msize ) && cnt > 50 )
			{
				useDP = true ;
				break ;
			}
		}
	}
	
	int atCnt = cnt ;
	printf( "%d: atCnt=%d seCnt=%d %d %d %d\n", subexons[0].start + 1, atCnt, seCnt, useDP, (int)constraints[0].constraints.size(), (int)constraints[0].matePairs.size() ) ;
	fflush( stdout ) ;
	std::vector<struct _transcript> alltranscripts ;
	
	if ( !useDP )
	{
		int origSize = atCnt ;
		alltranscripts.resize( atCnt ) ;
		for ( i = 0 ; i < atCnt ; ++i )
		{
			alltranscripts[i].seVector.Init( seCnt ) ; 
			alltranscripts[i].correlationScore = 1 ;
		}

		atCnt = 0 ;
		for ( i = 0 ; i < seCnt ; ++i )
		{
			if ( subexons[i].canBeStart )
				EnumerateTranscript( i, 0, f, 0, subexons, subexonCorrelation, 1, alltranscripts, atCnt ) ;
		}

		for ( i = atCnt ; i < origSize ; ++i )
			alltranscripts[i].seVector.Release() ;

		alltranscripts.resize( atCnt ) ;
		//printf( "transcript cnt: %d\n", atCnt ) ;
		//printf( "%d %d\n", alltranscripts[0].seVector.Test( 1 ), constraints[0].matePairs.size() ) ;
	}
	else // Use dynamic programming to pick a set of candidate transcript.
	{
		std::vector<struct _transcript> sampleTranscripts ;

		// pre allocate the memory.
		struct _dpAttribute attr ;
		attr.f1 = new struct _dp[seCnt] ;
		if ( seCnt <= 10000 )
		{
			attr.f2 = new struct _dp*[seCnt] ;
			for ( i = 0 ; i < seCnt ; ++i )
				attr.f2[i] = new struct _dp[seCnt] ;
		}
		else 
			attr.f2 = NULL ;
		
		hashMax = HASH_MAX ;
		if (seCnt > 500)
			hashMax = 1000003 ;
		else if (seCnt > 1000)
			hashMax = 10000019 ; 
		else if (seCnt > 1500)
			hashMax = 20000003 ;

		attr.hash = dpHash ; 
		if ( hashMax != HASH_MAX )
			attr.hash = new struct _dp[hashMax] ;

		for ( i = 0 ; i < seCnt ; ++i )
		{
			attr.f1[i].seVector.Nullify() ;
			attr.f1[i].seVector.Init( seCnt ) ;
			if ( attr.f2 )
			{
				for ( j = i ; j < seCnt ; ++j )
				{
					attr.f2[i][j].seVector.Nullify() ;
					attr.f2[i][j].seVector.Init( seCnt ) ;
				}
			}
		}
		for ( i = 0 ; i < hashMax ; ++i )
		{
			attr.hash[i].seVector.Nullify() ;
			attr.hash[i].seVector.Init( seCnt ) ;
		}

		// select candidate transcripts from each sample.
		struct _pair32 *sampleComplexity = new struct _pair32[ sampleCnt ] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			sampleComplexity[i].a = i ;
			sampleComplexity[i].b = constraints[i].constraints.size() ;
		}
		qsort( sampleComplexity, sampleCnt, sizeof( sampleComplexity[0] ), CompPairsByB ) ;
		int downsampleCnt = -1 ;

		for ( i = sampleCnt - 1 ; i >= 0 ; --i )
		{
			sampleTranscripts.clear() ;
			int iterBound = constraints[ sampleComplexity[i].a ].constraints.size() ;
			if ( i < sampleCnt - 1 )
				iterBound = 100 ;

			if ( i < sampleCnt - 10 && alltranscripts.size() > 1000 )
				iterBound = 10 ;
			//printf( "%d %d: %d %d %d %d\n", subexons[0].start + 1, sampleComplexity[i].a, constraints[ sampleComplexity[i].a ].constraints.size(), constraints[ sampleComplexity[i].a ].matePairs.size(),
			//		alltranscripts.size(), iterBound ) ; fflush( stdout ) ;	
			if ( maxDpConstraintSize > 0 )
			{
				Constraints truncatedConstraints ;
				truncatedConstraints.TruncateConstraintsCoverFrom( constraints[ sampleComplexity[i].a ], seCnt, maxDpConstraintSize ) ;
				PickTranscriptsByDP( subexons, seCnt, iterBound, truncatedConstraints, 
						subexonCorrelation, attr, sampleTranscripts ) ;		
			}
			else if ( ( constraints[ sampleComplexity[i].a ].constraints.size() > 1000 
				&& constraints[ sampleComplexity[i].a ].constraints.size() * 10 < constraints[ sampleComplexity[i].a ].matePairs.size() ) 
				|| ( downsampleCnt > 0 && (int)constraints[ sampleComplexity[i].a ].constraints.size() >= downsampleCnt ) 
				|| seCnt >= 1500 )
			{
				Constraints downsampledConstraints ;
				int stride = (int)constraints[ sampleComplexity[i].a ].matePairs.size() / (int)constraints[ sampleComplexity[i].a ].constraints.size() ;
				if ( downsampleCnt > 0 )
					stride = (int)constraints[ sampleComplexity[i].a ].constraints.size() / downsampleCnt ;
				if ( stride < 1 )
					stride = 1 ;
				downsampledConstraints.DownsampleConstraintsFrom( constraints[ sampleComplexity[i].a ], stride ) ; 
				if ( downsampleCnt <= 0 )
					downsampleCnt = downsampledConstraints.constraints.size() ;
				if ( iterBound <= 10 )
					continue ;
				PickTranscriptsByDP( subexons, seCnt, iterBound, downsampledConstraints, subexonCorrelation, attr, sampleTranscripts ) ;		
			}
			else
			{
				PickTranscriptsByDP( subexons, seCnt, iterBound, constraints[ sampleComplexity[i].a ], subexonCorrelation, attr, sampleTranscripts ) ;		
			}
			int size = sampleTranscripts.size() ;
			for ( j = 0 ; j < size ; ++j )
				alltranscripts.push_back( sampleTranscripts[j] ) ;
			
			// we can further pick a smaller subsets of transcripts here if the number is still to big.
			CoalesceSameTranscripts( alltranscripts ) ;
		
			AugmentTranscripts( subexons, alltranscripts, 1000, false ) ;
		}

		// release the memory.
		delete[] sampleComplexity ;
		for ( i = 0 ; i < seCnt ; ++i )	
		{
			attr.f1[i].seVector.Release() ;
			for ( j = i ; j < seCnt && attr.f2 ; ++j )
				attr.f2[i][j].seVector.Release() ;
		}
		for ( i = 0 ; i < hashMax ; ++i )
			attr.hash[i].seVector.Release() ;

		delete[] attr.f1 ;
		for ( i = 0 ; i < seCnt && attr.f2 ; ++i )
			delete[] attr.f2[i] ;
		if ( attr.f2 )
			delete[] attr.f2 ;
		if (hashMax != HASH_MAX)
			delete[] attr.hash ;

	}
	
	transcriptId = new int[usedGeneId - baseGeneId] ;
	std::vector<struct _transcript> *predTranscripts = new std::vector<struct _transcript>[sampleCnt] ;
	
	atCnt = alltranscripts.size() ;
	for ( i = 0 ; i < atCnt ; ++i )
		alltranscripts[i].FPKM = 0 ;
	
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = alltranscripts.size() ;
		for ( j = 0 ; j < size ; ++j )
			alltranscripts[j].abundance = -1 ;
		//printf( "pick: %d: %d %d\n", i, constraints[i].matePairs.size(), alltranscripts.size() ) ;
		PickTranscripts( subexons, alltranscripts, constraints[i], subexonCorrelation, predTranscripts[i] ) ;
		
		/*double tmp = FPKMFraction ;
		FPKMFraction = 0 ;
		size = predTranscripts.size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			ConvertTranscriptAbundanceToFPKM( subexons, predTranscripts[j] ) ;
		}
		RefineTranscripts( subexons, seCnt, predTranscripts, constraints[i] ) ;
		FPKMFraction = tmp ;*/
		
	}
	
	atCnt = alltranscripts.size() ;
	int *txptSampleSupport = new int[atCnt] ;
	memset( txptSampleSupport, 0, sizeof( int ) * atCnt ) ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = predTranscripts[i].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			++txptSampleSupport[ predTranscripts[i][j].id ] ;
			++alltranscripts[ predTranscripts[i][j].id ].FPKM ;
		}
	}

	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = alltranscripts.size() ;
		for ( j = 0 ; j < size ; ++j )
			alltranscripts[j].abundance = -1 ;
		//printf( "pick: %d: %d %d\n", i, constraints[i].matePairs.size(), alltranscripts.size() ) ;
		
		size = predTranscripts[i].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			predTranscripts[i][j].seVector.Release() ;
		}
		predTranscripts[i].clear() ;
		PickTranscripts( subexons, alltranscripts, constraints[i], subexonCorrelation, predTranscripts[i] ) ;
	}
	
	std::vector<int> *rawPredTranscriptIds = new std::vector<int>[sampleCnt] ;
	std::vector<double> *rawPredTranscriptAbundance = new std::vector<double>[sampleCnt] ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = predTranscripts[i].size() ;

		for ( j = 0 ; j < size ; ++j )
		{
			rawPredTranscriptIds[i].push_back( predTranscripts[i][j].id ) ;
			rawPredTranscriptAbundance[i].push_back( predTranscripts[i][j].abundance ) ;
		}
	}

	// Do the filtration.
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = predTranscripts[i].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			ConvertTranscriptAbundanceToFPKM( subexons, predTranscripts[i][j] ) ;
		}
		size = RefineTranscripts( subexons, seCnt, false, subexonChainSupport, txptSampleSupport, predTranscripts[i], constraints[i] ) ;
		
		// Recompute the abundance.
		AbundanceEstimation( subexons, seCnt, constraints[i], predTranscripts[i] ) ;
		for ( j = 0 ; j < size ; ++j )
			ConvertTranscriptAbundanceToFPKM( subexons, predTranscripts[i][j] ) ;
		size = RefineTranscripts( subexons, seCnt, true, subexonChainSupport, txptSampleSupport, predTranscripts[i], constraints[i] ) ;
		
		//ComputeTranscriptsScore( subexons, seCnt, subexonChainSupport, predTranscripts[i] ) ;
	}
	
	// Rescue some filtered transcripts
	memset( txptSampleSupport, 0, sizeof( int ) * atCnt ) ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		int size = predTranscripts[i].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			++txptSampleSupport[ predTranscripts[i][j].id ] ;
		}
	}

	bool *predicted = new bool[atCnt] ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		memset( predicted, false, sizeof( bool ) * atCnt ) ;
		if ( predTranscripts[i].size() != rawPredTranscriptIds[i].size() )	
		{
			int psize = predTranscripts[i].size() ;
			int rsize = rawPredTranscriptIds[i].size() ;
			int tcCnt = constraints[i].matePairs.size() ;

			for ( j = 0 ; j < psize ; ++j )
				predicted[ predTranscripts[i][j].id ] = true ;
			
			for ( j = 0 ; j < rsize ; ++j )
			{
				int id = rawPredTranscriptIds[i][j] ;
				if ( predicted[ id ] == false &&
					( txptSampleSupport[ id ] >= 3 && txptSampleSupport[id] >= 0.25 * sampleCnt ) )
				{
					struct _transcript nt = alltranscripts[id] ;
					nt.seVector.Nullify() ;
					nt.seVector.Duplicate( alltranscripts[id].seVector ) ;
					nt.constraintsSupport = NULL ;
					nt.correlationScore = -1 ;
					nt.abundance = rawPredTranscriptAbundance[i][j] ;
					nt.id = id ;
					predTranscripts[i].push_back( nt ) ;
				}
			}
			if ( psize != predTranscripts[i].size() )
				AbundanceEstimation( subexons, seCnt, constraints[i], predTranscripts[i] ) ;
		}

		int size = predTranscripts[i].size() ;

		if ( 0 ) //size == 1 )
		{
			//AugmentTranscripts( subexons, predTranscripts[i], false ) ;	
			
			int l = predTranscripts[i].size() ;
			int tcCnt = constraints[i].matePairs.size() ;
			for ( j = 0 ; j < l ; ++j )
			{
				predTranscripts[i][j].abundance = 1.0 / alignments.readLen ;
			}
			AbundanceEstimation( subexons, seCnt, constraints[i], predTranscripts[i] ) ;
			
			std::vector<int> subexonIdx ;
			for ( j = 0 ; j < l ; ++j )
			{
				subexonIdx.clear() ;
				predTranscripts[i][j].seVector.GetOnesIndices( subexonIdx ) ;
				int subexonIdxCnt = subexonIdx.size() ;
				int len = 0 ;
				for ( k = 0 ; k < subexonIdxCnt ; ++k )
					len += subexons[ subexonIdx[k] ].end - subexons[ subexonIdx[k] ].start + 1 ;

				if ( predTranscripts[i][j].abundance * alignments.readLen / len < 2.0 )
					predTranscripts[i][j].abundance = -1 ;
				else
					ConvertTranscriptAbundanceToFPKM( subexons, predTranscripts[i][j] ) ;

			}
			RemoveNegativeAbundTranscripts( predTranscripts[i] )  ;
		}
		
		// Output
		size = predTranscripts[i].size() ;
		InitTranscriptId() ;
		for ( j = 0 ; j < size ; ++j )
		{
			OutputTranscript( i, subexons, predTranscripts[i][j] ) ;
		}
		for ( j = 0 ; j < size ; ++j )
		{
			predTranscripts[i][j].seVector.Release() ;
		}
	}

	delete []predicted ;
	delete []transcriptId ;
	delete []predTranscripts ;
	delete []rawPredTranscriptIds ;
	delete []rawPredTranscriptAbundance ;
	delete []txptSampleSupport ;

	atCnt = alltranscripts.size() ;
	for ( i = 0 ; i < atCnt ; ++i )
		alltranscripts[i].seVector.Release() ;
	compatibleTestVectorT.Release() ;
	compatibleTestVectorC.Release() ;
	delete[] f ;
	delete[] subexonChainSupport ;
	return 0 ;	
}

void *TranscriptDeciderSolve_Wrapper( void *a ) 
{
	int i ;
	
	struct _transcriptDeciderThreadArg &arg = *( (struct _transcriptDeciderThreadArg *)a ) ;
	TranscriptDecider transcriptDecider( arg.FPKMFraction, arg.classifierThreshold, arg.txptMinReadDepth, arg.sampleCnt, *( arg.alignments ) ) ;
	transcriptDecider.SetNumThreads( arg.numThreads + 1 ) ;
	transcriptDecider.SetMultiThreadOutputHandler( arg.outputHandler ) ;
	transcriptDecider.SetMaxDpConstraintSize( arg.maxDpConstraintSize ) ;
	transcriptDecider.Solve( arg.subexons, arg.seCnt, arg.constraints, arg.subexonCorrelation ) ;
	
	int start = arg.subexons[0].start ;
	int end = arg.subexons[ arg.seCnt - 1 ].end ;
	int chrId = arg.subexons[0].chrId ;
	// Release memory
	for ( i = 0 ; i < arg.seCnt ; ++i )
	{
		delete[] arg.subexons[i].prev ;
		delete[] arg.subexons[i].next ;
	}
	delete[] arg.subexons ;

	// Put the work id back to the free threads queue.
	pthread_mutex_lock( arg.ftLock ) ;
	arg.freeThreads[ *( arg.ftCnt ) ] = arg.tid ;
	++*( arg.ftCnt ) ;
	if ( *( arg.ftCnt ) == 1 )
		pthread_cond_signal( arg.fullWorkCond ) ;
	pthread_mutex_unlock( arg.ftLock) ;
	printf( "Thread %d: %s %d %d finished.\n", arg.tid, arg.alignments->GetChromName(chrId), start + 1, end + 1 ) ;
	fflush( stdout ) ;


	pthread_exit( NULL ) ;
}
