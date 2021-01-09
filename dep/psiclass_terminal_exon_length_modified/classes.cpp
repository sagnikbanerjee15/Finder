#include <stdio.h>
#include <getopt.h>
#include <vector>
#include <pthread.h>

#include "alignments.hpp"
#include "SubexonGraph.hpp"
#include "SubexonCorrelation.hpp"
#include "Constraints.hpp"
#include "TranscriptDecider.hpp"

char usage[] = "./classes [OPTIONS]:\n"
	"Required:\n"
	"\t-s STRING: path to the subexon file.\n"
	"\t-b STRING: path to the BAM file.\n"
	"\t\tor\n"
	"\t--lb STRING: path to the file of the list of BAM files.\n"
	"Optional:\n"
	"\t-p INT: number of threads. (default: 1)\n"
	"\t-o STRING: the prefix of the output file. (default: not used)\n"
	"\t-c FLOAT: only use the subexons with classifier score <= than the given number. (default: 0.05)\n" 
	"\t-f FLOAT: filter the transcript from the gene if its abundance is lower than the given number percent of the most abundant one. (default: 0.05)\n"
	"\t-d FLOAT: filter the transcript whose average read depth is less than the given number. (default: 2.5)\n"
	"\t--ls STRING: path to the file of the list of single-sample subexon files. (default: not used)\n"
	"\t--hasMateIdSuffix: the read id has suffix such as .1, .2 for a mate pair. (default: false)\n"
	"\t--maxDpConstraintSize: the maximum number of subexons a constraint can cover in dynamic programming. (default: 7; -1 for inf)\n"
	"\t--primaryParalog: use primary alignment to retain paralog genes instead of unique alignments. (default: not used)\n"
	;

static const char *short_options = "s:b:f:o:d:p:c:h" ;
static struct option long_options[] =
	{
		{ "ls", required_argument, 0, 10000 },
		{ "lb", required_argument, 0, 10001 },
		{ "hasMateIdSuffix", no_argument, 0, 10002 },
		{ "primaryParalog", no_argument, 0, 10003 },
		{ "maxDpConstraintSize", required_argument, 0, 10004 },
		{ (char *)0, 0, 0, 0} 
	} ;

struct _getAlignmentsInfoThreadArg
{
	std::vector<Alignments> *pAlignmentFiles ;
	int numThreads ;
	int tid ;
} ;

struct _getConstraintsThreadArg
{
	std::vector<Constraints> *pMultiSampleConstraints ;	
	int numThreads ;
	int tid ;

	struct _subexon *subexons ;
	int seCnt ;
	int start, end ;
} ;

void *GetAlignmentsInfo_Thread( void *pArg )
{
	int i ;
	std::vector <Alignments> &alignmentFiles = *( ( (struct _getAlignmentsInfoThreadArg *)pArg )->pAlignmentFiles ) ;
	int tid = ( (struct _getAlignmentsInfoThreadArg *)pArg )->tid ;
	int numThreads = ( (struct _getAlignmentsInfoThreadArg *)pArg )->numThreads ;
	int size = alignmentFiles.size() ;
	for ( i = 0 ; i < size ; ++i )
	{
		if ( i % numThreads == tid )
		{
			alignmentFiles[i].GetGeneralInfo( true ) ;
			alignmentFiles[i].Rewind() ;
		}
	}
	
	pthread_exit( NULL ) ;
}


void *GetConstraints_Thread( void *pArg )
{
	int i ;
	struct _getConstraintsThreadArg &arg = *( (struct _getConstraintsThreadArg *)pArg ) ;
	std::vector<Constraints> &multiSampleConstraints = *( arg.pMultiSampleConstraints ) ;
	int tid = arg.tid ;
	int numThreads = arg.numThreads ;
	int size = multiSampleConstraints.size() ;
	for ( i = 0 ; i < size ; ++i )
	{
		if ( i % numThreads == tid )
			multiSampleConstraints[i].BuildConstraints( arg.subexons, arg.seCnt, arg.start, arg.end ) ;
	}
	pthread_exit( NULL ) ;
}

int main( int argc, char *argv[] )
{
	int i, j ;
	int size ;

	if ( argc <= 1 )
	{
		printf( "%s", usage ) ;
		return 0 ;
	}

	int c, option_index ; // For getopt
	option_index = 0 ;
	FILE *fpSubexon = NULL ;
	double FPKMFraction = 0.05 ; 
	double classifierThreshold ;
	double txptMinReadDepth = 2.5 ;
	char outputPrefix[1024] = "" ;
	int numThreads = 1 ;
	bool hasMateReadIdSuffix = false ;
	bool usePrimaryAsUnique = false ;
	int maxDpConstraintSize = 7 ;
	
	std::vector<Alignments> alignmentFiles ;
	SubexonCorrelation subexonCorrelation ;
	
	classifierThreshold = 0.05 ;
	while ( 1 )
	{
		c = getopt_long( argc, argv, short_options, long_options, &option_index ) ;
		if ( c == -1 )
			break ;

		if ( c == 's' )
		{
			fpSubexon = fopen( optarg, "r" ) ;
		}
		else if ( c == 'b' )
		{
			Alignments a ;
			a.Open( optarg ) ;
			//a.SetAllowClip( false ) ;
			alignmentFiles.push_back( a ) ;
		}
		else if ( c == 'c' )
		{
			classifierThreshold = atof( optarg ) ;
		}
		else if ( c == 'f' )
		{
			FPKMFraction = atof( optarg ) ; 		
		}
		else if ( c == 'o' )
		{
			strcpy( outputPrefix, optarg ) ;	
		}
		else if ( c == 'd' )
		{
			txptMinReadDepth = atof( optarg ) ;
		}
		else if ( c == 'p' )
		{
			numThreads = atoi( optarg ) ;
		}
		else if ( c == 10000 ) // the list of subexon files.
		{
			subexonCorrelation.Initialize( optarg ) ;
		}
		else if ( c == 10001 ) // the list of bam files.
		{
			FILE *fp = fopen( optarg, "r" ) ;
			char buffer[1024] ;
			while ( fgets( buffer, sizeof( buffer ), fp ) != NULL )
			{
				int len = strlen( buffer ) ;
				if ( buffer[len - 1] == '\n' )
				{
					buffer[len - 1] = '\0' ;
					--len ;

				}
				Alignments a ;
				a.Open( buffer ) ;
				alignmentFiles.push_back( a ) ;
			}
			fclose( fp ) ;
		}
		else if ( c == 10002 ) // the mate pair read id has suffix.
		{
			hasMateReadIdSuffix = true ;
		}
		else if ( c == 10003 ) // do not check unique alignment
		{
			usePrimaryAsUnique = true ;
		}
		else if ( c == 10004 ) // maxDpConstraintSize
		{
			maxDpConstraintSize = atoi(optarg) ;
		}
		else
		{
			printf( "%s", usage ) ;
			exit( 1 ) ;
		}
	}
	if ( fpSubexon == NULL )			
	{
		printf( "Cannot find combined subexon file.\n" ) ;
		exit( 1 ) ;
	}
	if ( alignmentFiles.size() < 1 )
	{
		printf( "Must use -b option to specify BAM files.\n" ) ;
		exit( 1 ) ;
	}


	if ( alignmentFiles.size() < 50 )
	{
		size = alignmentFiles.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			alignmentFiles[i].GetGeneralInfo( true ) ;
			alignmentFiles[i].Rewind() ;
		}
	}
	else
	{
		struct _getAlignmentsInfoThreadArg *args = new struct _getAlignmentsInfoThreadArg[ numThreads ] ;
		pthread_attr_t pthreadAttr ;
		pthread_t *threads ;
		
		pthread_attr_init( &pthreadAttr ) ;
		pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;

		threads = new pthread_t[ numThreads ] ;

		for ( i = 0 ; i < numThreads ; ++i )
		{
			args[i].pAlignmentFiles = &alignmentFiles ;
			args[i].tid = i ;
			args[i].numThreads = numThreads ;

			pthread_create( &threads[i], &pthreadAttr, GetAlignmentsInfo_Thread, &args[i] ) ;
		}

		for ( i = 0 ; i < numThreads ; ++i )
		{
			pthread_join( threads[i], NULL ) ;
		}

		pthread_attr_destroy( &pthreadAttr ) ;
		delete[] args ;
		delete[] threads ;
	}

	// Build the subexon graph
	SubexonGraph subexonGraph( classifierThreshold, alignmentFiles[0], fpSubexon ) ;
	subexonGraph.ComputeGeneIntervals() ;
	
	// Solve gene by gene
	int sampleCnt = alignmentFiles.size() ;
	std::vector<Constraints> multiSampleConstraints ;
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		Constraints constraints( &alignmentFiles[i] ) ;
		constraints.SetHasMateReadIdSuffix( hasMateReadIdSuffix ) ;
		constraints.SetUsePrimaryAsUnique( usePrimaryAsUnique ) ;
		multiSampleConstraints.push_back( constraints ) ;
	}
	MultiThreadOutputTranscript outputHandler( sampleCnt, alignmentFiles[0] ) ;
	outputHandler.SetOutputFPs( outputPrefix ) ;

	if ( numThreads <= 1 )
	{
		TranscriptDecider transcriptDecider( FPKMFraction, classifierThreshold, txptMinReadDepth, sampleCnt, alignmentFiles[0] ) ;
		
		transcriptDecider.SetMultiThreadOutputHandler( &outputHandler ) ;
		transcriptDecider.SetNumThreads( numThreads ) ;
		transcriptDecider.SetMaxDpConstraintSize( maxDpConstraintSize ) ;

		int giCnt = subexonGraph.geneIntervals.size() ;
		for ( i = 0 ; i < giCnt ; ++i )
		{
			struct _geneInterval gi = subexonGraph.geneIntervals[i] ;
			struct _subexon *intervalSubexons = new struct _subexon[ gi.endIdx - gi.startIdx + 1 ] ;
			subexonGraph.ExtractSubexons( gi.startIdx, gi.endIdx, intervalSubexons ) ;
			printf( "%d: %d %s %d %d\n", i, gi.endIdx - gi.startIdx + 1, 
					alignmentFiles[0].GetChromName( intervalSubexons[0].chrId ), 
					gi.start + 1, gi.end + 1 ) ;	
			fflush( stdout ) ;

			subexonCorrelation.ComputeCorrelation( intervalSubexons, gi.endIdx - gi.startIdx + 1, alignmentFiles[0] ) ;
			for ( j = 0 ; j < sampleCnt ; ++j )
				multiSampleConstraints[j].BuildConstraints( intervalSubexons, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	

			transcriptDecider.Solve( intervalSubexons, gi.endIdx - gi.startIdx + 1, multiSampleConstraints, subexonCorrelation ) ;

			for ( j = 0 ; j < gi.endIdx - gi.startIdx + 1 ; ++j )
			{
				delete[] intervalSubexons[j].prev ;
				delete[] intervalSubexons[j].next ;
			}
			delete[] intervalSubexons ;
		}
	}
	else // multi-thread case.
	{
		--numThreads ; // one thread is used for read in the data.
		
		// Allocate memory
		struct _transcriptDeciderThreadArg *pArgs = new struct _transcriptDeciderThreadArg[ numThreads ] ;
		pthread_mutex_t ftLock ;
		int *freeThreads ;
		int ftCnt ;
		pthread_cond_t fullWorkCond ;
		pthread_attr_t pthreadAttr ;
		pthread_t *threads ;
		pthread_t *getConstraintsThreads ;
		bool *initThreads ;

		pthread_mutex_init( &ftLock, NULL ) ;
		pthread_cond_init( &fullWorkCond, NULL ) ;
		pthread_attr_init( &pthreadAttr ) ;
		pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;

		threads = new pthread_t[ numThreads ] ;
		getConstraintsThreads = new pthread_t[ numThreads ] ;
		initThreads = new bool[numThreads] ;
		freeThreads = new int[ numThreads ] ;
		ftCnt = numThreads ;
		for ( i = 0 ; i < numThreads ; ++i )
		{
			pArgs[i].tid = i ;
			pArgs[i].sampleCnt = sampleCnt ;
			pArgs[i].numThreads = numThreads ;
			pArgs[i].maxDpConstraintSize = maxDpConstraintSize ;
			pArgs[i].FPKMFraction = FPKMFraction ;
			pArgs[i].classifierThreshold = classifierThreshold ;
			pArgs[i].txptMinReadDepth = txptMinReadDepth ;
			pArgs[i].alignments = &alignmentFiles[0] ;
			//pArgs[i].constraints = new std::vector<Constraints> ;
			pArgs[i].outputHandler = &outputHandler ;

			freeThreads[i] = i ;
			pArgs[i].freeThreads = freeThreads ;
			pArgs[i].ftCnt = &ftCnt ;
			pArgs[i].ftLock = &ftLock ;
			pArgs[i].fullWorkCond = &fullWorkCond ;
			
			for ( j = 0 ; j < sampleCnt ; ++j )
			{
				Constraints constraints( &alignmentFiles[j] ) ;
				pArgs[i].constraints.push_back( constraints ) ;
			}

			initThreads[i] = false ;
		}


		// Read in and distribute the work
		int giCnt = subexonGraph.geneIntervals.size() ;
		for ( i = 0 ; i < giCnt ; ++i )
		{
			struct _geneInterval gi = subexonGraph.geneIntervals[i] ;
			struct _subexon *intervalSubexons = new struct _subexon[ gi.endIdx - gi.startIdx + 1 ] ;
			subexonGraph.ExtractSubexons( gi.startIdx, gi.endIdx, intervalSubexons ) ;
			subexonCorrelation.ComputeCorrelation( intervalSubexons, gi.endIdx - gi.startIdx + 1, alignmentFiles[0] ) ;
			printf( "%d: %d %s %d %d. Free threads: %d/%d\n", i, gi.endIdx - gi.startIdx + 1, 
					alignmentFiles[0].GetChromName( intervalSubexons[0].chrId ), 
					gi.start + 1, gi.end + 1, ftCnt, numThreads + 1 ) ;	
			fflush( stdout ) ;
			
			int gctCnt = ftCnt ;
			if ( gctCnt > 1 && sampleCnt > 1 )
			{
				gctCnt = ( gctCnt < sampleCnt ? gctCnt : sampleCnt ) ;	
				struct _getConstraintsThreadArg *args = new struct _getConstraintsThreadArg[ gctCnt ] ;
				for ( j = 0 ; j < gctCnt ; ++j )
				{
					args[j].pMultiSampleConstraints = &multiSampleConstraints ;
					args[j].numThreads = gctCnt ;
					args[j].tid = j ;
					args[j].subexons = intervalSubexons ;
					args[j].seCnt = gi.endIdx - gi.startIdx + 1 ;
					args[j].start = gi.start ; args[j].end = gi.end ;
					pthread_create( &getConstraintsThreads[j], &pthreadAttr, GetConstraints_Thread, &args[j] ) ;
				}


				
				for ( j = 0 ; j < gctCnt ; ++j )
					pthread_join( getConstraintsThreads[j], NULL ) ;		

				delete[] args ;
			}
			else
			{
				for ( j = 0 ; j < sampleCnt ; ++j )
					multiSampleConstraints[j].BuildConstraints( intervalSubexons, gi.endIdx - gi.startIdx + 1, gi.start, gi.end ) ;	
			}
			
			// Search for the free queue.
			int tag = -1 ; // get the working thread.
			pthread_mutex_lock( &ftLock ) ;
			if ( ftCnt == 0 )
			{
				pthread_cond_wait( &fullWorkCond, &ftLock ) ;	
			}
			tag = freeThreads[ ftCnt - 1 ] ;
			--ftCnt ;
			pthread_mutex_unlock( &ftLock ) ;
			
			if ( initThreads[tag] )
				pthread_join( threads[tag], NULL ) ; // Make sure the chosen thread exits.
			
			// Assign the subexons, the constraints and correlation content.
			pArgs[tag].subexons = new struct _subexon[gi.endIdx - gi.startIdx + 1] ;
			pArgs[tag].seCnt = gi.endIdx - gi.startIdx + 1 ;
			for ( j = 0 ; j < pArgs[tag].seCnt ; ++j )
			{
				pArgs[tag].subexons[j] = intervalSubexons[j] ;
				int cnt = intervalSubexons[j].prevCnt ;
				pArgs[tag].subexons[j].prev = new int[cnt] ;
				memcpy( pArgs[tag].subexons[j].prev, intervalSubexons[j].prev, sizeof( int ) * cnt ) ;
				cnt = intervalSubexons[j].nextCnt ;
				pArgs[tag].subexons[j].next = new int[cnt] ;
				memcpy( pArgs[tag].subexons[j].next, intervalSubexons[j].next, sizeof( int ) * cnt ) ;
			}

			for ( j = 0 ; j < sampleCnt ; ++j )
			{
				pArgs[tag].constraints[j].Assign( multiSampleConstraints[ j ] ) ;
			}
			pArgs[tag].subexonCorrelation.Assign( subexonCorrelation ) ;
			pthread_create( &threads[tag], &pthreadAttr, TranscriptDeciderSolve_Wrapper, &pArgs[tag] ) ;
			initThreads[tag] = true ;
			for ( j = 0 ; j < gi.endIdx - gi.startIdx + 1 ; ++j )
			{
				delete[] intervalSubexons[j].prev ;
				delete[] intervalSubexons[j].next ;
			}
			delete[] intervalSubexons ;
		}

		for ( i = 0 ; i < numThreads ; ++i )
		{
			if ( initThreads[i] )
				pthread_join( threads[i], NULL ) ;
		}

		// Release memory
		for ( i = 0 ; i < numThreads ; ++i )
		{
			std::vector<Constraints>().swap( pArgs[i].constraints ) ;
		}
		delete []pArgs ;
		pthread_attr_destroy( &pthreadAttr ) ;
		pthread_mutex_destroy( &ftLock ) ;
		pthread_cond_destroy( &fullWorkCond ) ;

		delete[] threads ;
		delete[] getConstraintsThreads ;
		delete[] initThreads ;
		delete[] freeThreads ;
	} // end of else for multi-thread.

	outputHandler.OutputCommandInfo( argc, argv ) ; 
	for ( i = 0 ; i < sampleCnt ; ++i )
	{
		char buffer[1024] ;
		alignmentFiles[i].GetFileName( buffer ) ;
		int separator = -1 ;
		
		// deprive the directory.
		for ( j = 0 ; buffer[j] ; ++j )
		{
			if ( buffer[j] == '/' )
				separator = j ;
		}
		if ( separator != -1 )
		{
			++separator ;
			for ( j = separator ; buffer[j] ; ++j )
				buffer[j - separator] = buffer[j] ;
			buffer[j - separator] = '\0' ;
		}
		outputHandler.OutputCommentToSampleGTF( i, buffer ) ;
	}
	outputHandler.ComputeFPKMTPM( alignmentFiles ) ;
	outputHandler.Flush() ;
	
	for ( i = 0 ; i < sampleCnt ; ++i )
		alignmentFiles[i].Close() ;
	fclose( fpSubexon ) ;
	return 0 ;
}
