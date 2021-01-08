#ifndef _MOURISL_CLASSES_TRANSCRIPTDECIDER_HEADER
#define _MOURISL_CLASSES_TRANSCRIPTDECIDER_HEADER

#include <pthread.h>
#include <map>
#include <time.h>
#include <stdarg.h>

#include "alignments.hpp"
#include "SubexonGraph.hpp"
#include "SubexonCorrelation.hpp"
#include "BitTable.hpp"
#include "Constraints.hpp"

#define HASH_MAX 100003 // default HASH_MAX
#define USE_DP 200000

struct _transcript
{
	BitTable seVector ;
	double abundance ;
	double correlationScore ; 
	double FPKM ;
	double *constraintsSupport ; // record the assign ment of constraints.

	int first, last ; // indicate the index of the first and last subexons.
	bool partial ; // wehther this is a partial transcript.
	int id ; // the id for various usage: i.e transcript index in the alltranscripts.
} ;

struct _outputTranscript
{
	int chrId ;
	int geneId, transcriptId ;
	struct _pair32 *exons ;
	int ecnt ;
	char strand ;
	int sampleId ;

	double FPKM ;
	double TPM ;
	double cov ;
} ;

struct _dp
{
	BitTable seVector ; 
	int first, last ;
	// The "cnt" is for the hash structure.
	// the first cnt set bits represent the subexons that are the key of the hash
	// the remaining set bits are the optimal subtranscript follow the key.
	int cnt ; 
	double cover ;

	double minAbundance ;
	int timeStamp ;
	int strand ;
} ;


struct _dpAttribute
{
	struct _dp *f1, **f2 ;
	struct _dp *hash ;

	struct _transcript bufferTxpt ;

	bool forAbundance ;

	struct _subexon *subexons ;
	int seCnt ;

	std::map<uint64_t, int> uncoveredPair ;
	
	double minAbundance ;
	int timeStamp ;
} ;

class MultiThreadOutputTranscript ;

struct _transcriptDeciderThreadArg
{
	int tid ;
	struct _subexon *subexons ;
	int seCnt ;
	int sampleCnt ;
	int numThreads ;

	int maxDpConstraintSize ;
	double FPKMFraction, classifierThreshold, txptMinReadDepth ;
	Alignments *alignments ;
	std::vector<Constraints> constraints ;
	SubexonCorrelation subexonCorrelation ;
	MultiThreadOutputTranscript *outputHandler ;

	int *freeThreads ; // the stack for free threads
	int *ftCnt ;
	pthread_mutex_t *ftLock ;
	pthread_cond_t *fullWorkCond ;
} ;

class MultiThreadOutputTranscript
{
private:
	std::vector<struct _outputTranscript> outputQueue ;
	pthread_t *threads ;
	pthread_mutex_t outputLock ;
	int sampleCnt ;
	int numThreads ;
	std::vector<FILE *> outputFPs ;
	Alignments &alignments ;

public:
	static int CompTranscripts( const struct _outputTranscript &a, const struct _outputTranscript &b )
	{
		int i ;
		if ( a.geneId != b.geneId )
			return a.geneId - b.geneId ;
		if ( a.ecnt != b.ecnt )
			return a.ecnt - b.ecnt ;

		for ( i = 0 ; i < a.ecnt ; ++i )
		{
			if ( a.exons[i].a != b.exons[i].a )					
				return a.exons[i].a - b.exons[i].a ;

			if ( a.exons[i].b != b.exons[i].b )					
				return a.exons[i].b - b.exons[i].b ;
		}
		return 0 ;
	}

	static bool CompSortTranscripts( const struct _outputTranscript &a, const struct _outputTranscript &b )
	{
		int tmp = CompTranscripts( a, b ) ;
		if ( tmp < 0 )
			return true ;
		else if ( tmp > 0 )
			return false ;
		else
			return a.sampleId < b.sampleId ;
	}

	MultiThreadOutputTranscript( int cnt, Alignments &a ): alignments( a )
	{
		sampleCnt = cnt ;
		pthread_mutex_init( &outputLock, NULL ) ;
	}
	~MultiThreadOutputTranscript()
	{
		pthread_mutex_destroy( &outputLock ) ;
		int i ;
		for ( i = 0 ; i < sampleCnt ; ++i )
			fclose( outputFPs[i] ) ;
	}

	void SetThreadsPointer( pthread_t *t, int n )
	{
		threads = t ;
		numThreads = n ;
	}

	void SetOutputFPs( char *outputPrefix ) 
	{
		int i ;
		char buffer[1024] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			if ( outputPrefix[0] )
				sprintf( buffer, "%s_sample_%d.gtf", outputPrefix, i ) ;
			else
				sprintf( buffer, "sample_%d.gtf", i ) ;
			FILE *fp = fopen( buffer, "w" ) ;
			outputFPs.push_back( fp ) ;
		}
	}

	void Add( struct _outputTranscript &t ) 
	{
		pthread_mutex_lock( &outputLock ) ;
		outputQueue.push_back( t ) ;
		pthread_mutex_unlock( &outputLock ) ;
	}
	
	void Add_SingleThread( struct _outputTranscript &t ) 
	{
		outputQueue.push_back( t ) ;
	}
	
	void ComputeFPKMTPM( std::vector<Alignments> &alignmentFiles )
	{
		int i ;
		int qsize = outputQueue.size() ;
		double *totalFPK = new double[ sampleCnt ] ;
		memset( totalFPK, 0, sizeof( double ) * sampleCnt ) ;
		for ( i = 0 ; i < qsize ; ++i )
		{
			totalFPK[ outputQueue[i].sampleId ] += outputQueue[i].FPKM ;
		}

		for ( i = 0 ; i < qsize ; ++i )
		{
			outputQueue[i].TPM = outputQueue[i].FPKM / ( totalFPK[ outputQueue[i].sampleId ] / 1000000.0 ) ;
			outputQueue[i].FPKM /= ( alignmentFiles[ outputQueue[i].sampleId ].totalReadCnt / 1000000.0 ) ;
		}

		delete[] totalFPK ;
	}

	void OutputCommandInfo( int argc, char *argv[] )
	{
		int i ;
		int j ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			fprintf( outputFPs[i], "#PsiCLASS_v1.0.2\n#" ) ;	
			for ( j = 0 ; j < argc - 1 ; ++j )
			{
				fprintf( outputFPs[i], "%s ", argv[j] ) ;
			}
			fprintf( outputFPs[i], "%s\n", argv[j] ) ;
		}
	}

	void OutputCommentToSampleGTF( int sampleId, char *s )
	{
		fprintf( outputFPs[ sampleId ], "#%s\n", s ) ;
	}

	void Flush()
	{
		std::sort( outputQueue.begin(), outputQueue.end(), CompSortTranscripts ) ;	
		int i, j ;
		int qsize = outputQueue.size() ;
		char prefix[10] = "" ;

		// Recompute the transcript id
		int gid = -1 ;
		int tid = 0 ;
		for ( i = 0 ; i < qsize ; )
		{
			for ( j = i + 1 ; j < qsize ; ++j )
			{
				if ( CompTranscripts( outputQueue[i], outputQueue[j] ) )
					break ;
			}
			int l ;	
			if ( outputQueue[i].geneId != gid )
			{
				gid = outputQueue[i].geneId ;
				tid = 0 ;
			}
			else
				++tid ;

			for ( l = i ; l < j ; ++l )
				outputQueue[l].transcriptId = tid ;

			i = j ;
		}

		// output
		for ( i = 0 ; i < qsize ; ++i )
		{
			struct _outputTranscript &t = outputQueue[i] ;
			char *chrom = alignments.GetChromName( t.chrId ) ;

			fprintf( outputFPs[t.sampleId], "%s\tPsiCLASS\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; FPKM \"%.6lf\"; TPM \"%.6lf\"; cov \"%.6lf\";\n",
					chrom, t.exons[0].a, t.exons[t.ecnt - 1].b, t.strand,
					prefix, chrom, t.geneId,
					prefix, chrom, t.geneId, t.transcriptId, t.FPKM, t.TPM, t.cov ) ;
			for ( j = 0 ; j < t.ecnt ; ++j )
			{
				fprintf( outputFPs[ t.sampleId ], "%s\tPsiCLASS\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; "
						"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; FPKM \"%.6lf\"; TPM \"%.6lf\"; cov \"\%.6lf\";\n",
						chrom, t.exons[j].a, t.exons[j].b, t.strand,
						prefix, chrom, t.geneId,
						prefix, chrom, t.geneId, t.transcriptId,
						j + 1, t.FPKM, t.TPM, t.cov ) ;
			}
			delete []t.exons ;
		}
	}
} ;

class TranscriptDecider
{
private:
	int sampleCnt ;
	int numThreads ;
	double FPKMFraction ;
	double txptMinReadDepth ;
	int hashMax ;
	int maxDpConstraintSize ;

	Constraints *constraints ;
	//struct _subexon *subexons ;
	//int seCnt ;

	int usedGeneId ;
	int baseGeneId, defaultGeneId[2] ;

	int *transcriptId ; // the next transcript id for each gene id (we shift the gene id to 0 in this array.)
	Alignments &alignments ; // for obtain the chromosome names.

	std::vector<FILE *> outputFPs ;

	BitTable compatibleTestVectorT, compatibleTestVectorC ;
	double canBeSoftBoundaryThreshold ;

	MultiThreadOutputTranscript *outputHandler ;

	// Test whether subexon tag is a start subexon in a mixture region that corresponds to the start of a gene on another strand.
	bool IsStartOfMixtureStrandRegion( int tag, struct _subexon *subexons, int seCnt ) ;

	// The functions to pick transcripts through dynamic programming
	struct _dp *dpHash ;
	void SearchSubTranscript( int tag, int strand, int parents[], int pcnt, struct _dp &pdp, int visit[], int vcnt, int extends[], int extendCnt, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) ;
	struct _dp SolveSubTranscript( int visit[], int vcnt, int strand, std::vector<struct _constraint> &tc, int tcStartInd, struct _dpAttribute &attr ) ;
	void PickTranscriptsByDP( struct _subexon *subexons, int seCnt, int iterBound, Constraints &constraints, SubexonCorrelation &correlation, struct _dpAttribute &attr, std::vector<struct _transcript> &allTranscripts ) ;

	void SetDpContent( struct _dp &a, struct _dp &b, const struct _dpAttribute &attr )
	{
		a.seVector.Assign( b.seVector ) ;
		a.first = b.first ;
		a.last = b.last ;
		a.cnt = b.cnt ;
		a.cover = b.cover ;
		
		a.strand = b.strand ;
		a.minAbundance = attr.minAbundance ;
		a.timeStamp = attr.timeStamp ;
	} 

	void ResetDpContent( struct _dp &d )
	{
		d.seVector.Reset() ;
		d.first = -1 ;
		d.last = -1 ;
		d.cnt = -1 ;
		d.cover = -1 ;
		d.minAbundance = -1 ;
		d.timeStamp = -1 ;
	}

	void AugmentTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &alltranscripts, int limit, bool extend ) ;
	// Test whether a constraints is compatible with the transcript.
	// Return 0 - uncompatible or does not overlap at all. 1 - fully compatible. 2 - Head of the constraints compatible with the tail of the transcript
	int IsConstraintInTranscript( struct _transcript transcript, struct _constraint &c ) ;
	int IsConstraintInTranscriptDebug( struct _transcript transcript, struct _constraint &c ) ;
	
	// Count how many transcripts are possible starting from subexons[tag].
	int SubTranscriptCount( int tag, struct _subexon *subexons, int f[] ) ;

	// The methods when there is no need for DP
	void EnumerateTranscript( int tag, int strand, int visit[], int vcnt, struct _subexon *subexons, SubexonCorrelation &correlation, double correlationScore, std::vector<struct _transcript> &alltranscripts, int &atcnt ) ;
	// For the simpler case, we can pick sample by sample.
	void PickTranscripts( struct _subexon *subexons, std::vector<struct _transcript> &alltranscripts, Constraints &constraints, SubexonCorrelation &seCorrelation, std::vector<struct _transcript> &transcripts ) ; 
	
	static bool CompSortTranscripts( const struct _transcript &a, const struct _transcript &b )
	{
		if ( a.first < b.first )
			return true ;
		else if ( a.first > b.first )
			return false ;

		int diffPos = a.seVector.GetFirstDifference( b.seVector ) ;
		if ( diffPos == -1 )
			return false ;
		if ( a.seVector.Test( diffPos ) )
			return true ;
		else
			return false ;
	} 
	
	static bool CompSortPairs( const struct _pair32 &x, const struct _pair32 &y )
	{
		if ( x.a != y.a )
			return x.a < y.a ;
		else 
			return x.b < y.b ;
	}

	static bool CompSortPairsByB( const struct _pair32 &x, const struct _pair32 &y )
	{
		return x.b < y.b ;
	}

	static int CompPairsByB( const void *p1, const void *p2 )
	{
		return ((struct _pair32 *)p1)->b - ((struct _pair32 *)p2)->b ;
	}

	double ComputeScore( double cnt, double weight, double a, double A, double correlation )
	{
		if ( a > A * 0.1 )
			return ( cnt * weight ) * ( 1 + pow( a / A, 0.25 ) ) + correlation ;
		else 
			return ( cnt * weight ) * ( 1 + a / A ) + correlation ;
		//return ( cnt ) * ( exp( 1 + a / A ) ) + correlation ; 
	}

	int GetFather( int f, int *father ) ;
	
	void ConvertTranscriptAbundanceToFPKM( struct _subexon *subexons, struct _transcript &t, int readCnt = 1000000 )
	{
		int txptLen = 0 ;
		int i, size ;

		std::vector<int> subexonInd ;
		t.seVector.GetOnesIndices( subexonInd ) ;
		size = subexonInd.size() ;
		for ( i = 0 ; i < size ; ++i )
			txptLen += ( subexons[ subexonInd[i] ].end - subexons[ subexonInd[i] ].start + 1 ) ;
		double factor = 1 ;
		if ( alignments.matePaired )
			factor = 0.5 ;
		t.FPKM = t.abundance * factor / ( ( readCnt / 1000000.0 ) * ( txptLen / 1000.0 ) ) ;
	}

	int GetTranscriptLengthFromAbundanceAndFPKM( double abundance, double FPKM, int readCnt = 1000000 )
	{
		double factor = 1 ;
		if ( alignments.matePaired )
			factor = 0.5 ;
		return int( abundance * factor / ( FPKM / 1000.0 ) / ( readCnt / 1000000.0  ) + 0.5 ) ;
	}

	void CoalesceSameTranscripts( std::vector<struct _transcript> &t ) ;


	// Initialize the structure to store transcript id 
	void InitTranscriptId() ; 
	
	int GetTranscriptGeneId( std::vector<int> &subexonInd, struct _subexon *subexons ) ;
	int GetTranscriptGeneId( struct _transcript &t, struct _subexon *subexons ) ;
	int RemoveNegativeAbundTranscripts( std::vector<struct _transcript> &transcripts )
	{
		int i, j ;
		int tcnt = transcripts.size() ;
		j = 0 ;
		for ( i = 0 ; i < tcnt ; ++i )
		{
			if ( transcripts[i].abundance < 0 )
			{
				transcripts[i].seVector.Release() ; // Don't forget release the memory.
				continue ;
			}
			transcripts[j] = transcripts[i] ;
			++j ;
		}
		transcripts.resize( j ) ;
		return j ;
	}
		
	void AbundanceEstimation( struct _subexon *subexons, int seCnt, Constraints &constraints, std::vector<struct _transcript> &transcripts ) ;

	int RefineTranscripts( struct _subexon *subexons, int seCnt, bool aggressive, std::map<int, int> *subexonChainSupport, int *txptSampleSupport, std::vector<struct _transcript> &transcripts, Constraints &constraints ) ;

	void ComputeTranscriptsScore( struct _subexon *subexons, int seCnt, std::map<int, int> *subexonChainSupport, std::vector<struct _transcript> &transcripts ) ;

	void OutputTranscript( int sampleId, struct _subexon *subexons, struct _transcript &transcript ) ;

	void PrintLog( const char *fmt, ... )
	{
		char buffer[10021] ;
		va_list args ;
		va_start( args, fmt ) ;
		vsprintf( buffer, fmt, args ) ;

		time_t mytime = time(NULL) ;
		struct tm *localT = localtime( &mytime ) ;
		char stime[500] ;
		strftime( stime, sizeof( stime ), "%c", localT ) ;
		fprintf( stderr, "[%s] %s\n", stime, buffer ) ;
	}


public:
	TranscriptDecider( double f, double c, double d, int sampleCnt, Alignments &a ): alignments( a )  
	{
		FPKMFraction = f ;
		canBeSoftBoundaryThreshold = c ;
		txptMinReadDepth = d ;
		usedGeneId = 0 ;
		defaultGeneId[0] = -1 ;
		defaultGeneId[1] = -1 ;
		maxDpConstraintSize = -1 ;
		numThreads = 1 ;
		this->sampleCnt = sampleCnt ;
		dpHash = new struct _dp[ HASH_MAX ] ; // pre-allocated buffer to hold dp information.
	}
	~TranscriptDecider() 
	{
		int i ;
		if ( numThreads == 1 )
		{
			int size = outputFPs.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				fclose( outputFPs[i] ) ;
			}
		}
		delete[] dpHash ;
	}


	// @return: the number of assembled transcript 
	int Solve( struct _subexon *subexons, int seCnt, std::vector<Constraints> &constraints, SubexonCorrelation &subexonCorrelation ) ;

	void SetOutputFPs( char *outputPrefix )
	{
		int i ;
		char buffer[1024] ;
		for ( i = 0 ; i < sampleCnt ; ++i )
		{
			if ( outputPrefix[0] )
				sprintf( buffer, "%s_sample_%d.gtf", outputPrefix, i ) ;
			else
				sprintf( buffer, "sample_%d.gtf", i ) ;
			FILE *fp = fopen( buffer, "w" ) ;
			outputFPs.push_back( fp ) ;
		}
	}

	void SetMultiThreadOutputHandler( MultiThreadOutputTranscript *h ) 
	{
		outputHandler = h ;
	}

	void SetNumThreads( int t )
	{
		numThreads = t ;
	}

	void SetMaxDpConstraintSize(int size)
	{
		maxDpConstraintSize = size ;
	}
} ;

void *TranscriptDeciderSolve_Wrapper( void *arg ) ;

#endif
