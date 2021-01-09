// Compare the reference and the predictions
// Format: ./a.out ref.gtf prediction.gtf
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_TXPT 2000000

#define DEBUG 0

bool flagMultiExonOnly = false ;
bool flagValidChr = true ;
int relaxWidth = 20 ;

struct _pair
{
	int a, b ;
} ;

struct _transcript
{
	char tid[50] ;
	char chrom[50] ;
	char strand ;
	struct _pair *exons ;
	int ecnt ;
} ;

struct _info
{
	double coverage ;
	int byId ; // It gets this coverage by comparing with byId.
} ;

struct _intron
{
	char chrom[50] ;
	int start, end  ;
	int tindex ; // The index to the corresponding transcripts
} ;

struct _exon
{
	char chrom[50] ;
	int start, end ;
	//int tindex ; 
	int soft ; // left-bit: left side is soft, right-bit: right side is soft
	bool matched ;
} ;

int TranscriptComp( const void *p1, const void *p2 )
{
	const struct _transcript *a = (struct _transcript *)p1 ;
	const struct _transcript *b = ( struct _transcript *)p2 ;
	int tmp = strcmp( a->chrom, b->chrom ) ;
	if ( tmp != 0 )
		return tmp ;
	return a->exons[0].a - b->exons[0].a ;
}

int TranscriptComp_ByIntron( const void *p1, const void *p2 )
{
	const struct _transcript *a = (struct _transcript *)p1 ;
	const struct _transcript *b = ( struct _transcript *)p2 ;
	if ( a->strand != b->strand )
		return a->strand - b->strand ;
	int tmp = strcmp( a->chrom, b->chrom ) ;
	if ( tmp != 0 )
		return tmp ;
	if ( a->ecnt != b->ecnt ) 
		return a->ecnt - b->ecnt ;
	int i ;
	for ( i = 0 ; i < a->ecnt - 1 ; ++i )
	{
		if ( a->exons[i].b != b->exons[i].b )
			return a->exons[i].b - b->exons[i].b ;
		if ( a->exons[i + 1].a != b->exons[i + 1].a )
			return a->exons[i + 1].a - b->exons[i + 1].a ;
	}

	return 0 ; //strcmp( a->tid, b->tid ) ;
}

bool validChrom( char *chrom )
{
	if ( !flagValidChr )
		return true ;
	if ( chrom[0] != 'c' || chrom[1] != 'h' || chrom[2] != 'r' )
		return false ;
	// only consider chr1-22,x,y,z
	if ( chrom[3]=='x' || chrom[3]=='X' 
		|| chrom[3] == 'y' || chrom[3] == 'Y' 
		|| chrom[3] == 'm' || chrom[3] == 'M'
		|| ( chrom[3] >= '3' && chrom[3] <= '9' ) )
	{
		if ( chrom[4] == '\0' )
			return true ;
		else
			return false ;
	}
	if ( chrom[3] == '1' )
	{
		if ( chrom[4] == '\0' )
			return true ;
		if ( chrom[4] >= '0' && chrom[4] <= '9' && chrom[5] == '\0' )
			return true ;
		return false ;
	}
	if ( chrom[3] == '2' )
	{
		if ( chrom[4] == '\0' )
			return true ;
		if ( chrom[4] >= '0' && chrom[4] <= '2' && chrom[5] == '\0' )
			return true ;
		return false ;
	}
	return false ;
}

void ReverseExonList( struct _pair *exons, int ecnt )
{
	if ( ecnt < 2 || ( exons[0].a < exons[1].a ) )
		return ;
	int i, j ;
	struct _pair tmp ;

	for ( i = 0, j = ecnt - 1 ; i < j ; ++i, --j )
	{
		tmp = exons[i] ;
		exons[i] = exons[j] ;
		exons[j] = tmp ;
	}
}

/**
  Merge exons that are next to each other. Showed up in scripture.
*/
int CleanExonList( struct _pair *exons, int ecnt )
{
	int i, k ;
	for ( i = 0 ; i < ecnt - 1 ; ++i )
		if ( exons[i + 1].a - exons[i].b - 1 < 20 )
		{
			exons[i + 1].a = exons[i].a ;
			exons[i].a = -1 ;
		}
	k = 0 ;
	for ( i = 0 ; i < ecnt ; ++i )
	{
		if ( exons[i].a == -1 )
			continue ;
		exons[k] = exons[i] ;
		++k ;
	}
	return k ;
}

/**
  Remove the duplicated intron chain in the list
*/
int RemoveDuplicateIntronChain( struct _transcript *t, int tcnt )
{
	int i, j, k ;

	qsort( t, tcnt, sizeof( *t ), TranscriptComp_ByIntron ) ;
	for ( i = 1 ; i < tcnt ; ++i )
	{
		k = i - 1 ;
		while ( t[k].ecnt == -1 )
			--k ;

		if ( t[i].ecnt != t[k].ecnt || strcmp( t[i].chrom, t[k].chrom ) || t[i].strand != t[k].strand )
			continue ;
		for ( j = 0 ; j < t[i].ecnt - 1 ; ++j )
			if ( t[i].exons[j].b != t[k].exons[j].b || 
				t[i].exons[j + 1].a != t[k].exons[j + 1].a )
				break ;
		if ( j >= t[i].ecnt - 1 && t[i].ecnt != 1 )
			t[i].ecnt = -1 ;
		if ( t[i].ecnt == 1 && t[i].exons[0].a == t[k].exons[0].a &&
			t[i].exons[0].b == t[k].exons[0].b )
			t[i].ecnt = -1 ;
		/*if ( t[i].ecnt == -1 )
		{
			printf( "%s <- %s\n", t[i].tid, t[k].tid ) ;
		}*/
	}

	k = 0 ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		if ( t[i].ecnt == -1 )
		{
			free( t[i].exons ) ;
			continue ;
		}
		t[k] = t[i] ;
		++k ; 
	}
	tcnt = k ;
	return tcnt ;
}


double CompareTranscripts( const struct _transcript &ref, const struct _transcript &pred )
{
	int i, j ;
	struct _pair *refExons = ref.exons ;
	int refECnt = ref.ecnt ;
	struct _pair *predExons = pred.exons  ;
	int predECnt = pred.ecnt ;

	// Prediction must be a "subset" of the reference
	if ( refECnt < predECnt )
		return -1 ;

	// single exon case
	if ( refECnt == 1 )
	{
		if ( refExons[0].b < predExons[0].a || refExons[0].a > predExons[0].b )
			return -1 ;
		else
			return 1 ;
	}
	
	if ( predECnt == 1 )
	{
		for ( i = 0 ; i < refECnt ; ++i )
		{
			if ( refExons[i].b < predExons[0].a || refExons[i].a > predExons[0].b )
				continue ;
			else
			{
				/*if ( i == 0 && predExons[0].a >= refExons[0].a - relaxWidth )
					return 0 ;
				else if ( i == refECnt - 1 && predExons[0].b <= refExons[i].b + relaxWidth )
					return 0 ;
				else if ( i > 0 && i < refECnt - 1 && predExons[0].a >= refExons[i].a - relaxWidth && 
					predExons[0].b <= refExons[i].b + relaxWidth )
					return 0 ;

				return -1 ;*/
				return 0 ;
			}
		}
		return -1 ;
	}

	if ( predECnt > 0 && ref.strand != pred.strand )
		return -1 ;

	// Test the intron chain
	for ( i = 0 ; i < refECnt - 1 ; ++i )
		if ( refExons[i].b == predExons[0].b )
			break ;
	if ( i >= refECnt - 1 )
		return -1 ;
	//if ( i != 0 && predExons[0].a < refExons[i].a - relaxWidth )
	//	return -1 ;

	for ( j = 0 ; i < refECnt - 1 && j < predECnt - 1 ; ++i, ++j )
	{
		if ( refExons[i].b != predExons[j].b || refExons[i + 1].a != predExons[j + 1].a )	
			break ;
	}
	if ( j >= predECnt - 1 )
	{
		/*if ( i == refECnt - 1 ||
			predExons[ predECnt - 1 ].b <= refExons[i].b + relaxWidth )
			return (double)( predECnt - 1 ) / (double)( refECnt - 1 ) ;
		else
			return -1 ;*/
		return (double)( predECnt - 1 ) / (double)( refECnt - 1 ) ;
	}
	else
		return -1 ;
}

bool WithinRange( int a, int b, int w = 10 )
{
	if ( a - b >= -w && a - b <= w )
		return true ;
	return false ;
}

int GetIntrons( struct _transcript *transcripts, int tcnt, struct _intron *introns )
{
	int i, j, k, tag ;
	int cnt = 0 ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		for ( j = 0 ; j < transcripts[i].ecnt - 1 ; ++j )
		{
			bool flag = false ;

			tag = 0 ;
			for ( k = cnt - 1 ; k >= 0 ; --k )
			{
				if ( strcmp( introns[k].chrom, transcripts[i].chrom ) )
				{
					//printf( "hi\n" ) ;
					tag = k + 1 ;
					break ;
				}
				if ( introns[k].start < transcripts[i].exons[j].b )
				{
					tag = k + 1 ;
					break ;
				}
				else if ( introns[k].start == transcripts[i].exons[j].b &&
					introns[k].end == transcripts[i].exons[j + 1].a )
				{
					flag = true ;
					break ;
				}
			}

			if ( !flag )
			{
				for ( k = cnt ; k > tag ; --k )
					introns[k] = introns[k - 1] ;
				strcpy( introns[k].chrom, transcripts[i].chrom ) ;
				introns[k].start = transcripts[i].exons[j].b ;
				introns[k].end = transcripts[i].exons[j + 1].a ;
				introns[k].tindex = i ;
				++cnt ;
			}
		}
	}
	return cnt ;
}

int GetExons( struct _transcript *transcripts, int tcnt, struct _exon *exons  )
{
	int i, j, k, tag ;
	int cnt = 0 ;
	for ( i = 0 ; i < tcnt ; ++i )
	{
		for ( j = 0 ; j < transcripts[i].ecnt ; ++j )
		{
			bool flag = false ;
			int soft = 0 ;
			if ( j == 0 )
				soft = soft | 2 ;
			if ( j == transcripts[i].ecnt - 1 )
				soft = soft | 1 ;

			tag = 0 ;
			for ( k = cnt - 1 ; k >= 0 ; --k )
			{
				if ( strcmp( exons[k].chrom, transcripts[i].chrom ) )
				{
					//printf( "hi\n" ) ;
					tag = k + 1 ;
					break ;
				}

				if ( exons[k].start < transcripts[i].exons[j].a )
				{
					tag = k + 1 ;
					break ;
				}
				else if ( exons[k].start == transcripts[i].exons[j].a &&
					exons[k].end == transcripts[i].exons[j].b &&
					exons[k].soft == soft )
				{
					flag = true ;
					break ;
				}
			}

			if ( !flag )
			{
				for ( k = cnt ; k > tag ; --k )
					exons[k] = exons[k - 1] ;
				strcpy( exons[k].chrom, transcripts[i].chrom ) ;
				exons[k].start = transcripts[i].exons[j].a ;
				exons[k].end = transcripts[i].exons[j].b ;
				//exons[k].tindex = i ;
				exons[k].soft = soft ;
				exons[k].matched = false ;
				++cnt ;
			}
		}
	}
	return cnt ;
}

//chr1	HAVANA	transcript	320162	324461	.	+	.	gene_id "ENSG00000237094.6"; transcript_id "ENST00000423728.1"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP4-669L17.10"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "RP4-669L17.10-002"; level 2; havana_gene "OTTHUMG00000156968.4"; havana_transcript "OTTHUMT00000346879.1";
int main( int argc, char *argv[] )
{
	FILE *fp ;
	struct _transcript *ref, *pred ;		
	struct _info *refInfo, *predInfo ;

	int refCnt, predCnt ;
	char line[10000], filename[10000] ;
	struct _pair tmpExons[10000] ;
	int tmpECnt = 0 ;

	char chrom[50], tool[20], type[40], strand[3] ;
	char tid[50] ;
	char buffer[50] ;
	int start, end ;

	int i, j, k, tag ;

	if ( argc == 1 )
	{
		printf( "Compare the reference transcripts and predicted transcripts.\n"
			"Format: ./grader ref.gtf prediction.gtf\n" ) ;
		exit( 1 ) ;
	}
	// Parse the arguments
	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-M" ) )
			flagMultiExonOnly = true ;
		else if ( !strcmp( argv[i], "-ac" ) )
			flagValidChr = false ;
		else
		{
			printf( "Unknown argument: %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	/*========================================================================================
	Stage 1: Read the transcripts from the reference and prediction's GTF files.
	========================================================================================*/
	fp = NULL ;
	fp = fopen( argv[1], "r" ) ;
	if ( fp == NULL )
	{
		printf( "Could not open file %s.\n", argv[1] ) ;
		exit( 1 ) ;
	}
	refCnt = 0 ;
	ref = ( struct _transcript *)malloc( MAX_TXPT * sizeof( struct _transcript ) ) ;
	
	strcpy( ref[0].tid, "tmp_-1" ) ;
	while ( fgets( line, sizeof( line ), fp ) != NULL )
	{
		if ( line[0] == '#' )
			continue ;

		sscanf( line, "%s %s %s %d %d %s %s", chrom, tool, type, &start, &end, buffer, strand ) ;
		
		if ( strcmp( type, "exon" ) || !validChrom( chrom ) )
			continue ;

		char *p = strstr( line, "transcript_id" ) ;
		for ( ; *p != ' ' ; ++p )
			;
		p += 2 ;
		sscanf( p, "%s", tid ) ;
		//printf( "+%s %d\n", tid, strlen( tid )  ) ;
		p = tid + strlen( tid ) ;
		while ( *p != '\"' )
			--p ;
		*p = '\0' ;
		//printf( "%s\n", tid ) ;
		if ( strcmp( tid, ref[ refCnt ].tid ) && tmpECnt )
		{
			ReverseExonList( tmpExons, tmpECnt ) ;
			tmpECnt = CleanExonList( tmpExons, tmpECnt ) ;
			if ( tmpECnt > 1 || !flagMultiExonOnly )
			{
				ref[ refCnt ].ecnt = tmpECnt ;
				ref[ refCnt ].exons = ( struct _pair * )malloc( sizeof( struct _pair ) * tmpECnt ) ; 
				memcpy( ref[ refCnt ].exons, tmpExons, sizeof( tmpExons[0] ) * tmpECnt ) ;
				++refCnt ;
			}
			tmpECnt = 0 ;
		}

		tmpExons[ tmpECnt ].a = start ;
		tmpExons[ tmpECnt ].b = end ;
		++tmpECnt ;
		
		ref[ refCnt ].strand = strand[0] ;
		strcpy( ref[ refCnt ].chrom, chrom ) ;
		strcpy( ref[ refCnt ].tid, tid ) ;
	}
	if ( tmpECnt != 0 )
	{
		ReverseExonList( tmpExons, tmpECnt ) ;
		tmpECnt = CleanExonList( tmpExons, tmpECnt ) ;
		if ( tmpECnt > 1 || !flagMultiExonOnly )
		{
			ref[ refCnt ].ecnt = tmpECnt ;
			ref[ refCnt ].exons = ( struct _pair * )malloc( sizeof( struct _pair ) * tmpECnt ) ; 
			memcpy( ref[ refCnt ].exons, tmpExons, sizeof( tmpExons[0] ) * tmpECnt ) ;
			++refCnt ;
		}
		tmpECnt = 0 ;
	}
	fclose( fp ) ;

	fp = NULL ;
	fp = fopen( argv[2], "r" ) ;
	if ( fp == NULL )
	{
		printf( "Could not open file %s.\n", argv[2] ) ;
		exit( 1 ) ;
	}
	predCnt = 0 ;
	pred = ( struct _transcript *)malloc( MAX_TXPT * sizeof( struct _transcript ) ) ;
	
	strcpy( pred[0].tid, "tmp_-1" ) ;
	tmpECnt = 0 ;
	while ( fgets( line, sizeof( line ), fp ) != NULL )
	{
		if ( line[0] == '#' )
			continue ;
		sscanf( line, "%s %s %s %d %d %s %s", chrom, tool, type, &start, &end, buffer, strand ) ;
		
		if ( strcmp( type, "exon" ) || !validChrom( chrom ) )
			continue ;

		char *p = strstr( line, "transcript_id" ) ;
		//char *p = strstr( line, "gene_name" ) ;
		for ( ; *p != ' ' ; ++p )
			;
		p += 2 ;
		sscanf( p, "%s", tid ) ;
		//tid[ strlen( tid ) - 2 ] = '\0' ;
		p = tid + strlen( tid ) ;
		while ( *p != '\"' )
			--p ;
		*p = '\0' ;

		if ( strcmp( tid, pred[ predCnt ].tid ) && tmpECnt )
		{
			ReverseExonList( tmpExons, tmpECnt ) ;
			tmpECnt = CleanExonList( tmpExons, tmpECnt ) ;
			if ( tmpECnt > 1 || !flagMultiExonOnly )
			{
				pred[ predCnt ].ecnt = tmpECnt ;
				pred[ predCnt ].exons = ( struct _pair * )malloc( sizeof( struct _pair ) * tmpECnt ) ; 
				memcpy( pred[ predCnt ].exons, tmpExons, sizeof( tmpExons[0] ) * tmpECnt ) ;
				++predCnt ;
			}
			tmpECnt = 0 ;
		}
		
		tmpExons[ tmpECnt ].a = start ;
		tmpExons[ tmpECnt ].b = end ;
		++tmpECnt ;
		
		pred[ predCnt ].strand = strand[0] ;
		strcpy( pred[ predCnt ].chrom, chrom ) ;
		strcpy( pred[ predCnt ].tid, tid ) ;
	}
	if ( tmpECnt > 0 )
	{
		ReverseExonList( tmpExons, tmpECnt ) ;
		tmpECnt = CleanExonList( tmpExons, tmpECnt ) ;
		if ( tmpECnt > 1 || !flagMultiExonOnly )
		{
			pred[ predCnt ].ecnt = tmpECnt ;
			pred[ predCnt ].exons = ( struct _pair * )malloc( sizeof( struct _pair ) * tmpECnt ) ; 
			memcpy( pred[ predCnt ].exons, tmpExons, sizeof( tmpExons[0] ) * tmpECnt ) ;
			++predCnt ;
		}
		tmpECnt = 0 ;
	}
	
	refCnt = RemoveDuplicateIntronChain( ref, refCnt ) ;
	//printf( "predCnt = %d\n", predCnt ) ;	
	predCnt = RemoveDuplicateIntronChain( pred, predCnt ) ;
	//printf( "predCnt = %d\n", predCnt ) ;	
	qsort( ref, refCnt, sizeof( ref[0] ), TranscriptComp ) ;
	qsort( pred, predCnt, sizeof( pred[0] ), TranscriptComp ) ;

	
	/*========================================================================================
	  Stage 2: Compute the recall and precision.
	  ========================================================================================*/
	refInfo = ( struct _info * )malloc( refCnt * sizeof( *refInfo ) ) ;
	predInfo = ( struct _info * )malloc( predCnt * sizeof( *predInfo ) ) ;
	for ( i = 0 ; i < refCnt ; ++i )
	{
		refInfo[i].coverage = -1 ;
		refInfo[i].byId = -1 ;
#if DEBUG	
		printf( "=%d %s %s %d\n", i, ref[i].chrom, ref[i].tid, ref[i].ecnt ) ;
		for ( j = 0 ; j < ref[i].ecnt ; ++j )
			printf( "%d %d\n", ref[i].exons[j].a, ref[i].exons[j].b ) ;
#endif

	}
	for ( i = 0 ; i < predCnt ; ++i )
	{
		predInfo[i].coverage = -1 ;
		predInfo[i].byId = -1 ;
#if DEBUG	
		printf( "#%d %s %s %d\n", i, pred[i].chrom, pred[i].tid, pred[i].ecnt ) ;
		for ( j = 0 ; j < pred[i].ecnt ; ++j )
			printf( "%d %d\n", pred[i].exons[j].a, pred[i].exons[j].b ) ;
#endif	
	}
	tag = 0 ;
	for ( i = 0 ; i < predCnt ; ++i )
	{
		while ( tag < refCnt && 
				( ( strcmp( ref[tag].chrom, pred[i].chrom ) < 0 )  || 
				  ( !strcmp( ref[tag].chrom, pred[i].chrom  ) && 
				    ref[tag].exons[ ref[tag].ecnt - 1 ].b < pred[i].exons[0].a - relaxWidth ) ) )
			++tag ;

		//	if ( i == 16474 )
		//		printf( "%d %d %d", i, tag, refCnt ) ;
		for ( j = tag ; j < refCnt && ref[j].exons[0].a <= pred[i].exons[ pred[i].ecnt - 1 ].b + relaxWidth && !strcmp( ref[j].chrom, pred[i].chrom ); ++j )
		{
			/*if ( !strcmp( pred[i].tid, "CUFF.4256.1" ) )
			{
				printf( "%s ? %s\n", pred[i].tid, ref[j].tid ) ;
			}*/
			//		if ( i == 16474 )
			//		{
			//			printf( "%d %d\n", i, j ) ;
			//		}
			double coverage = CompareTranscripts( ref[j], pred[i] ) ;
			if ( coverage > refInfo[j].coverage )
			{
				refInfo[j].coverage = coverage ;
				refInfo[j].byId = i ;
			}
			if ( coverage > predInfo[i].coverage )
			{
				predInfo[i].coverage = coverage ;
				predInfo[i].byId = j ;
			}
		}
	}

	/*========================================================================================
	  Stage 3: Dump the information into output files.
	  ========================================================================================*/
	//sprintf( filename, "grader.%s.recall", argv[1] ) ;
	sprintf( filename, "grader.recall" ) ;
	fp = fopen( filename, "w" ) ;
	for ( i = 0 ; i < refCnt ; ++i )
	{
		fprintf( fp, "%s\t%d\t%lf\t%s\n", ref[i].tid, ref[i].ecnt, refInfo[i].coverage, 
			refInfo[i].byId == -1 ? "-" :  pred[ refInfo[i].byId ].tid  ) ;		
	}
	fclose( fp ) ;

	//sprintf( filename, "grader.%s.precision", argv[2] ) ;
	sprintf( filename, "grader.precision" ) ;
	fp = fopen( filename, "w" ) ;
	for ( i = 0 ; i < predCnt ; ++i )
	{
		fprintf( fp, "%s\t%d\t%lf\t%s\n", pred[i].tid, pred[i].ecnt, predInfo[i].coverage,
			predInfo[i].byId == -1 ? "-" : ref[predInfo[i].byId].tid ) ;		
	}
	fclose( fp ) ;

	// print the summary to the stdout
	const int binCnt = 20 ;
	int bins[binCnt + 1] ;
	memset( bins, 0, sizeof( bins ) ) ;
	k = 0 ;
	for ( i = 0 ; i < refCnt ; ++i )
	{
		if ( refInfo[i].coverage == -1 )
			continue ;
		++bins[ (int)( refInfo[i].coverage * binCnt ) ] ;
		++k ;
	}
	for ( i = 0 ; i <= binCnt ; ++i )
	{
		printf( "recall %lf %lf %d/%d\n", (double)i / binCnt, (double)k / refCnt, k, refCnt ) ;
		k -= bins[i] ;
	}

	memset( bins, 0, sizeof( bins ) ) ;
	k = 0 ;
	for ( i = 0 ; i < predCnt ; ++i )
	{
		if ( predInfo[i].coverage == -1 )
			continue ;
		++bins[ (int)( predInfo[i].coverage * binCnt ) ] ;
		++k ;
	}
	for ( i = 0 ; i <= binCnt ; ++i )
	{
		printf( "precision %lf %lf %d/%d\n", (double)i / binCnt, (double)k / predCnt, k, predCnt ) ;
		k -= bins[i] ;
	}

	/*========================================================================================
	Stage 4: Evaluations for introns.
	========================================================================================*/
	struct _intron *refIntrons = ( struct _intron * )malloc( sizeof( struct _intron ) * MAX_TXPT ) ;
	int riCnt = GetIntrons( ref, refCnt, refIntrons ) ;
	struct _intron *predIntrons = ( struct _intron * )malloc( sizeof( struct _intron ) * MAX_TXPT ) ;
	int piCnt = GetIntrons( pred, predCnt, predIntrons )  ;
	int matchedIntron = 0 ;
	/*for ( i = 0 ; i < riCnt ; ++i )
	{
		printf( "%d %s %d %d\n", i, refIntrons[i].chrom, refIntrons[i].start, refIntrons[i].end ) ;
	}*/

	fp = fopen( "grader.intron", "w" ) ;
	tag = 0 ;
	for ( i = 0 ; i < piCnt ; ++i )
	{
		bool flag = false ;
		while ( 1 )
		{
			if ( tag >= riCnt )
				break ;
			int tmp = strcmp( refIntrons[tag].chrom, predIntrons[i].chrom ) ;			
			if ( tmp < 0 || ( tmp == 0 && refIntrons[tag].start < predIntrons[i].start ) )
				++tag ;
			else
				break ;
		}
		//printf( "%d (%s %d %d) %d=>", i, predIntrons[i].chrom, predIntrons[i].start, predIntrons[i].end, tag ) ;
		/*if ( predIntrons[i].start == 1613457 )
			printf( "%d %d\n", tag, riCnt ) ;*/
		for ( j = tag ; j < riCnt ; ++j )
		{
			if ( strcmp( refIntrons[j].chrom, predIntrons[i].chrom ) > 0 || 
				refIntrons[j].start > predIntrons[i].start /*+ 10*/ )	
				break ;
			if ( refIntrons[j].start == predIntrons[i].start &&
				refIntrons[j].end == predIntrons[i].end )
			//if ( WithinRange( refIntrons[j].start, predIntrons[i].start ) &&
			//	WithinRange( refIntrons[j].end, predIntrons[i].end ) )
			{
				++matchedIntron ;
				flag = true ;
				break ;
			}
		}
		//printf( "%d\n", j ) ;
		if ( flag )
			fprintf( fp, "Y\t%s\t%d\t%d\t%s\n", predIntrons[i].chrom, predIntrons[i].start, predIntrons[i].end, pred[ predIntrons[i].tindex ].tid ) ;
		else
			fprintf( fp, "N\t%s\t%d\t%d\t%s\n", predIntrons[i].chrom, predIntrons[i].start, predIntrons[i].end, pred[ predIntrons[i].tindex ].tid  ) ;
	}
	fclose( fp ) ;
	printf( "\n" ) ;
	printf( "Intron evaluation\n") ;
	//printf( "\tmatched\t%d\n", matchedIntron ) ;
	printf( "\trecall\t%lf\t%d/%d\n", matchedIntron / (double)riCnt, matchedIntron, riCnt ) ;
	printf( "\tprecision\t%lf\t%d/%d\n", matchedIntron / (double)piCnt, matchedIntron, piCnt ) ;

	/*========================================================================================
	Stage 4: Evaluations for exons.
	========================================================================================*/
	struct _exon *refExons = ( struct _exon * )malloc( sizeof( struct _exon ) * MAX_TXPT ) ;
	int reCnt = GetExons( ref, refCnt, refExons ) ;
	struct _exon *predExons = ( struct _exon * )malloc( sizeof( struct _exon ) * MAX_TXPT ) ;
	int peCnt = GetExons( pred, predCnt, predExons )  ;
	int matchedExons = 0, matchedInternalExons = 0 ;
	int refInternalExons = 0, predInternalExons = 0 ; 

	for ( i = 0 ; i < reCnt ; ++i )
	{
		//printf( "%d %d\n", refExons[i].start, refExons[i].end ) ;
		if ( refExons[i].soft == 0 )
			++refInternalExons ;
		//if ( refExons[i].soft == 3 )
		//	printf( "hi\n" ) ;
	}
	for ( i = 0 ; i < peCnt ; ++i )
		if ( predExons[i].soft == 0 )
			++predInternalExons ;
	tag = 0 ;
	for ( i = 0 ; i < peCnt ; ++i )
	{
		bool flag = false ;
		while ( 1 )
		{
			if ( tag >= reCnt )
				break ;
			int tmp = strcmp( refExons[tag].chrom, predExons[i].chrom ) ;			
			if ( tmp < 0 || ( tmp == 0 && refExons[tag].end < predExons[i].start ) )
				++tag ;
			else
				break ;
		}
		for ( j = tag ; j < reCnt ; ++j )
		{
			if ( strcmp( refExons[j].chrom, predExons[i].chrom ) > 0 || 
				refExons[j].start > predExons[i].end /*+ 10*/ )	
				break ;
			if ( refExons[j].soft != predExons[i].soft 
				|| refExons[j].end < predExons[i].start )
				continue ;

			if ( refExons[j].start == predExons[i].start &&
				refExons[j].end == predExons[i].end && predExons[i].soft == 0 )
			//if ( WithinRange( refIntrons[j].start, predIntrons[i].start ) &&
			//	WithinRange( refIntrons[j].end, predIntrons[i].end ) )
			{
				refExons[j].matched = true ;
				predExons[i].matched = true ;
				++matchedInternalExons ;
				flag = true ;
				break ;
			}
			else if ( ( refExons[j].start == predExons[i].start || ( predExons[i].soft & 2 ) ) 
				&& ( refExons[j].end == predExons[i].end || (predExons[i].soft & 1 ) ) )
			{
				refExons[j].matched = true ;
				predExons[i].matched = true ;
				flag = true ;
				//break ;
			}
		}
		//printf( "%d\n", j ) ;
	}

	printf( "\n" ) ;
	printf( "Exon evaluation\n" ) ;
	//printf( "\tmatched\t%d\n", matchedExons ) ;
	matchedExons = 0 ;
	for ( i = 0 ; i < reCnt ; ++i )
		if ( refExons[i].matched )
			++matchedExons ;
	printf( "\trecall\t%lf\t%d/%d\n", (double)matchedExons / reCnt, matchedExons, reCnt ) ;
	matchedExons = 0 ;
	for ( i = 0 ; i < peCnt ; ++i )
		if ( predExons[i].matched )
			++matchedExons ;
	printf( "\tprecision\t%lf\t%d/%d\n", (double)matchedExons / peCnt, matchedExons, peCnt ) ;

	printf( "\n" ) ;
	printf( "Internal exon evaluation\n" ) ;
	//printf( "\tmatched\t%d\n", matchedInternalExons ) ;
	printf( "\trecall\t%lf\t%d/%d\n", (double)matchedInternalExons / refInternalExons, matchedInternalExons, refInternalExons ) ;
	printf( "\tprecision\t%lf\t%d/%d\n", (double)matchedInternalExons / predInternalExons, matchedInternalExons, predInternalExons ) ;
	return 0 ;
}
