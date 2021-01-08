#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <string>
#include <map>
#include <vector>
#include <algorithm>

char usage[] = "Usage: ./add-genename annotation.gtf gtflist [OPTIONS]\n"
		"\t-o: the directory for output gtf files (default: ./)\n" ;

struct _interval
{
	char chrom[97] ;
	char geneName[97] ;

	int start, end ;
} ;

char line[10000] ;
char buffer[10000], buffer2[10000], buffer3[10000] ;

int CompInterval( const struct _interval &a, const struct _interval &b )
{
	int chromCmp = strcmp( a.chrom, b.chrom ) ;
	if ( chromCmp )
		return chromCmp ;
	else if ( a.start != b.start )
		return a.start - b.start ;
	else
		return a.end - b.end ;
}

bool Comp( const struct _interval &a, const struct _interval &b )
{
	int cmp = CompInterval( a, b ) ;
	if ( cmp < 0 )
		return true ;
	else
		return false ;

	return false ;
}

void SortAndCleanIntervals( std::vector<struct _interval> &it )
{
	std::sort( it.begin(), it.end(), Comp ) ; 

	int i, k ;
	int size = it.size() ;
	if ( size == 0 )
		return ;

	k = 1 ;
	for ( i = 1 ; i < size ; ++i )
	{
		if ( CompInterval( it[k - 1], it[i] ) )
		{
			it[k] = it[i] ;
			++k ;
		}
	}
	it.resize( k ) ;
} 

int GetGTFField( char *ret, char *line, const char *fid )
{
	char *p = strstr( line, fid ) ;
	if ( p == NULL )
		return 0 ;
	//printf( "%s %s\n", line, fid ) ;
	for ( ; *p != ' ' ; ++p )
		;
	p += 2 ;
	
	sscanf( p, "%s", ret ) ;
	p = ret + strlen( ret ) ;
	while ( p != ret && *p != '\"' )
		--p ;
	*p = '\0' ;
	return 1 ;
}

void UpdateIdToAnnoId( std::vector<struct _interval> &transcript, char *tid, int exonTag, int intronTag, std::vector<struct _interval> &exons, std::vector<struct _interval> &introns,
	std::map<std::string, std::string> &txptIdToAnnoId, std::map<std::string, std::string> &geneIdToAnnoId )
{
	int i, k ;
	bool flag = false ;
	// First, try whether intron works.
	int ecnt = transcript.size() ;
	int intronSize = introns.size() ;
	int exonSize = exons.size() ;

	struct _interval itron ;
	strcpy( itron.chrom, transcript[0].chrom ) ;
	k = intronTag ;
	for ( i = 0 ; i < ecnt - 1 && !flag ; ++i )
	{
		itron.start = transcript[i].end ;
		itron.end = transcript[i + 1].start ;
		for ( ; k < intronSize ; ++k )
		{
			//printf( "hi %d: %s %d %d; %d %s %d %d\n", 0, itron.chrom, itron.start, itron.end, k, introns[k].chrom, introns[k].start, introns[k].end ) ;
			if ( strcmp( introns[k].chrom, itron.chrom ) )
				break ;

			int cmp = CompInterval( introns[k], itron ) ;
			if ( cmp == 0 )
			{
				txptIdToAnnoId[ std::string( tid ) ] = introns[k].geneName ;
				geneIdToAnnoId[ std::string( transcript[0].geneName ) ] = introns[k].geneName ;
				flag = true ;
				break ;
			}
			else if ( cmp > 0 )
				break ;
		}
	}

	// Next, try whether exon works
	k = exonTag ;
	for ( i = 0 ; i < ecnt && !flag ; ++i )
	{
		for ( ; k < exonSize ; ++k )
		{
			if ( strcmp( exons[k].chrom, transcript[i].chrom ) )
				break ;

			if ( exons[k].end < transcript[i].start )
				continue ;
			else if ( exons[k].start > transcript[i].end )
				break ;
			else
			{
				txptIdToAnnoId[ std::string( tid ) ] = exons[k].geneName ;
				geneIdToAnnoId[ std::string( transcript[0].geneName ) ] = exons[k].geneName ;
				flag = true ;
				break ;
			}
		}
	}
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	FILE *fp ;
	FILE *fpList ;
	char outputPath[1024] = "./" ;
	char prevTid[97] = "" ;

	std::vector<struct _interval> introns, exons ;
	std::map<std::string, std::string> txptIdToAnnoId, geneIdToAnnoId ;
	std::map<std::string, int> chromIdMapIntrons, chromIdMapExons ; // the index for the starting of a chrom.
	
	if ( argc < 3 )
	{	
		fprintf( stderr, "%s", usage ) ;
		exit( 1 ) ;
	}

	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-o" ) )
		{
			strcpy( outputPath, argv[i + 1 ] ) ;	
			++i ;
		}
		else
		{
			fprintf( stderr, "Unknown argument: %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	// Get the exons, introns from the annotation.
	fp = fopen( argv[1], "r" ) ;
	while ( fgets( line, sizeof( line ), fp ) != NULL )
	{
		if ( line[0] == '#' )
			continue ;
		char tid[97] ;
		char gname[97] ;
		char chrom[97] ;
		char type[50] ;
		int start, end ;

		sscanf( line, "%s %s %s %d %d", chrom, buffer, type, &start, &end ) ;
		if ( strcmp( type, "exon" ) )
			continue ;
		if ( GetGTFField( gname, line, "gene_name" ) )
		{
			struct _interval ne ;
			strcpy( ne.chrom, chrom ) ;
			strcpy( ne.geneName, gname ) ;
			ne.start = start ; ne.end = end ;

			exons.push_back( ne ) ;
		}
		else
		{
			fprintf( stderr, "%s has no field of gene_name.\n", line ) ;
			exit( 1 ) ;
		}

		if ( GetGTFField( tid, line, "transcript_id" ) )
		{
			if ( !strcmp( tid, prevTid ) )
			{	
				struct _interval ni ;
				
				strcpy( ni.chrom, chrom ) ;
				strcpy( ni.geneName, gname ) ;	

				int size = exons.size() ;
				if ( size < 2 )
					continue ;

				struct _interval &e1 = exons[ size - 2 ] ;
				struct _interval &e2 = exons[ size - 1 ] ;
				if ( e1.start < e2.start )
				{
					ni.start = e1.end ;
					ni.end = e2.start ;
				}
				else
				{
					ni.start = e2.end ;
					ni.end = e1.start ;
				}

				introns.push_back( ni ) ;			
			}
			else
				strcpy( prevTid, tid ) ;
		}
		else
		{
			fprintf( stderr, "%s has no field of transcript_id.\n", line ) ;
			exit( 1 ) ;
		}
	}
	fclose( fp ) ;
	SortAndCleanIntervals( exons ) ;
	SortAndCleanIntervals( introns ) ;
	int exonSize = exons.size() ;
	int intronSize = introns.size() ;
	// Determine the offset for each chrom on the list of features.	
	if ( exonSize )
	{
		chromIdMapExons[ std::string( exons[0].chrom ) ] = 0 ;
		for ( i = 1 ; i < exonSize ; ++i )
		{
			if ( strcmp( exons[i].chrom, exons[i - 1].chrom ) )
				chromIdMapExons[ std::string( exons[i].chrom ) ] = i ;
		}
	}

	if ( intronSize )
	{
		chromIdMapIntrons[ std::string( introns[0].chrom ) ] = 0 ;
		for ( i = 1 ; i < intronSize ; ++i )
		{
			if ( strcmp( introns[i].chrom, introns[i - 1].chrom ) )
			{
				//printf( "%s %d\n", introns[i].chrom, i ) ;
				chromIdMapIntrons[ std::string( introns[i].chrom ) ] = i ;
			}
		}

		//for ( i = 0 ; i < intronSize ; ++i )
		//	printf( "%s\t%d\t%d\n", introns[i].chrom, introns[i].start, introns[i].end ) ;
	}
	//printf( "%d %d\n", exonSize, intronSize ) ;

	// Go through all the GTF files to find the map between PsiCLASS's gene_id to annotation's gene name
	fpList = fopen( argv[2], "r ") ;
	std::vector<struct _interval> transcript ;
	int exonTag = 0 ;
	int intronTag = 0 ;

	while ( fscanf( fpList, "%s", line ) != EOF )
	{
		// Set the output file
		//printf( "hi: %s\n", line ) ;
		char *p ;
		p = line + strlen( line ) ;
		while ( p != line && *p != '/' )
		{
			if ( *p == '\n' )
				*p = '\0' ;
			--p ;
		}
		sprintf( buffer, "%s/%s", outputPath, p ) ;

		// Test whether this will overwrite the input gtf file.
		if ( realpath( buffer, buffer2 ) != NULL )
		{
			if ( realpath( line, buffer3 ) && !strcmp( buffer2, buffer3 ) )
			{
				fprintf( stderr, "Output will overwrite the input files. Please use -o to specify a different output directory.\n" ) ;
				exit( 1 ) ;
			}
		}

		fp = fopen( line, "r" ) ;
		transcript.resize( 0 ) ; // hold the exons in the transcript.
		prevTid[0] = '\0' ;
		exonTag = 0 ;
		intronTag = 0 ;
		int farthest = 0 ;

		while ( fgets( line, sizeof( line ), fp ) )
		{
			char tid[97] ;
			//char gname[97] ;
			char chrom[97] ;
			char type[50] ;
			int start, end ;

			if ( line[0] == '#' )
				continue ;

			sscanf( line, "%s %s %s %d %d", chrom, buffer, type, &start, &end ) ;
			if ( strcmp( type, "exon" ) )
				continue ;

			if ( GetGTFField( tid, line, "transcript_id" ) )
			{
				struct _interval ne ;
				strcpy( ne.chrom, chrom ) ;
				ne.start = start ;
				ne.end = end ;
				GetGTFField( ne.geneName, line, "gene_id" ) ;
				
				if ( !strcmp( tid, prevTid ) )
				{
					transcript.push_back( ne ) ;	
					if ( end > farthest )
						farthest = end ;
				}
				else if ( transcript.size() > 0 )
				{
					// the non-existed chrom will be put to offset 0, and that's fine
					if ( strcmp( transcript[0].chrom, introns[intronTag].chrom ) )
						intronTag = chromIdMapIntrons[ std::string( transcript[0].chrom ) ] ; 
					if ( strcmp( transcript[0].chrom, introns[intronTag].chrom ) )
						exonTag = chromIdMapIntrons[ std::string( transcript[0].chrom ) ] ;

					UpdateIdToAnnoId( transcript, prevTid, exonTag, intronTag, exons, introns, txptIdToAnnoId, geneIdToAnnoId ) ;

					// Adjust the offset if we are done with a gene cluster
					// We don't need to worry about the case of changing chrom here.
					if ( !strcmp( ne.geneName, transcript[0].geneName ) &&
						start > farthest ) // Use farthest to avoid inteleaved gene.
					{
						while ( intronTag < intronSize && !strcmp( introns[intronTag ].chrom, ne.chrom )
							&& introns[intronTag].end < ne.start )		
							++intronTag ;

						while ( exonTag < exonSize && !strcmp( exons[ exonTag ].chrom, ne.chrom ) 
							&& exons[ exonTag ].end < ne.start )
							++exonTag ;

						farthest = end ;
					}

					transcript.resize( 0 ) ;
					transcript.push_back( ne ) ;
					
					// Find the overlaps.
					strcpy( prevTid, tid ) ;
				}
				else
					strcpy( prevTid, tid ) ;
			}
			else
			{
				fprintf( stderr, "Could not find transcript_id field in GTF file: %s", line ) ;
				exit( 1 ) ;
			}
		}

		if ( transcript.size() > 0 )
		{
			if ( strcmp( transcript[0].chrom, introns[intronTag].chrom ) )
				intronTag = chromIdMapIntrons[ std::string( transcript[0].chrom ) ] ; 
			if ( strcmp( transcript[0].chrom, introns[intronTag].chrom ) )
				exonTag = chromIdMapIntrons[ std::string( transcript[0].chrom ) ] ;

			UpdateIdToAnnoId( transcript, prevTid, exonTag, intronTag, exons, introns, txptIdToAnnoId, geneIdToAnnoId ) ;
		}
		fclose( fp ) ;
	}
	fclose( fpList ) ;

	// Add the gene_name field.
	fpList = fopen( argv[2], "r" ) ;
	int novelCnt = 0 ;
	while ( fscanf( fpList, "%s", line ) != EOF )
	{
		FILE *fpOut ;
		// Set the output file
		char *p ;
		p = line + strlen( line ) ;
		while ( p != line && *p != '/' )
		{
			if ( *p == '\n' )
				*p = '\0' ;
			--p ;
		}
		sprintf( buffer, "%s/%s", outputPath, p ) ;
		fpOut = fopen( buffer, "w" ) ;

		// Process the input file
		fp = fopen( line, "r" ) ;

		transcript.resize( 0 ) ;
		prevTid[0] = '\0' ;
		while ( fgets( line, sizeof( line ), fp ) )
		{
			char gname[97] ;
			if ( line[0] == '#' )
			{
				fprintf( fpOut, "%s", line ) ;
				continue ;
			}
			
			if ( GetGTFField( buffer, line, "transcript_id" ) )
			{
				if ( txptIdToAnnoId.find( std::string( buffer ) ) != txptIdToAnnoId.end() )
				{
					int len = strlen( line ) ;
					if ( line[len - 1] == '\n' )
					{
						line[len - 1] = '\0' ;
						--len ;
					}
					fprintf( fpOut, "%s gene_name \"%s\";\n", line, txptIdToAnnoId[std::string( buffer )].c_str() ) ;
					continue ;
				}
			}

			if ( GetGTFField( gname, line, "gene_id" ) )
			{
				int len = strlen( line ) ;
				if ( line[len - 1] == '\n' )
				{
					line[len - 1] = '\0' ;
					--len ;
				}
				if ( geneIdToAnnoId.find( std::string( gname ) ) != geneIdToAnnoId.end() )
				{
					fprintf( fpOut, "%s gene_name \"%s\";\n", line, geneIdToAnnoId[std::string(gname)].c_str() ) ;
				}
				else
				{
					sprintf( buffer, "novel_%d", novelCnt ) ;
					geneIdToAnnoId[ std::string( gname ) ] = std::string( buffer ) ;
					fprintf( fpOut, "%s gene_name \"%s\";\n", line, buffer ) ;
					++novelCnt ;
				}

			}
			else
				fprintf( fpOut, "%s", line ) ;

		}
			
		fclose( fp ) ;
		fclose( fpOut ) ;
	}

	return 0 ;
}
