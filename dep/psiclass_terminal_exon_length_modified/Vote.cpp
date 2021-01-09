// The program that vote to pick the trusted transcripts.
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string.h>
#include <string>

#include "defs.h"
#include "TranscriptDecider.hpp"

char usage[] = "./transcript-vote [OPTIONS] > output.gtf:\n"
	"Required:\n"
	"\t--lg: path to the list of GTF files.\n"
	"Optional:\n" 
	"\t-d FLOAT: threshold of average coverage depth across all the samples. (default: 1)\n"
	//"\t-n INT: the number of samples a transcript showed up. (default: 3)\n"
	;


/*struct _transcript
{
	int chrId ;
	int geneId ; 
	char strand ;
	struct _pair32 *exons ;
	int ecnt ;
	int sampleId ;
} ;*/

void GetGTFField( char *s, const char *field, char *ret )
{
	char *p = strstr( s, field ) ;
	if ( p == NULL )
		return ;
	for ( ; *p != ' ' ; ++p )
		;
	p += 2 ; // add extra 1 to skip \"
	sscanf( p, "%s", ret ) ;
	//printf( "+%s %d\n", tid, strlen( tid )  ) ;
	p = ret + strlen( ret ) ;
	while ( *p != '\"' )
		--p ;
	*p = '\0' ;
}

int GetTailNumber( char *s )
{
	int len = strlen( s ) ;
	int ret = 0 ;
	int i ;
	int factor = 1 ;
	for ( i = len - 1 ; i >= 0 && s[i] >= '0' && s[i] <= '9' ; --i, factor *= 10 )
	{
		ret += factor * ( s[i] - '0' ) ;
	}
	return ret ;
}

int CompDouble( const void *p1, const void *p2 )
{
	double a = *(double *)p1 ;
	double b = *(double *)p2 ;
	if ( a > b )
		return 1 ;
	else if ( a < b )
		return -1 ;
	else
		return 0 ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	double minAvgDepth = 1.0 ;
	double fraction = 1.0 ;
	int minSampleCnt = 3 ;
	std::map<std::string, int> chrNameToId ;
	std::map<int, std::string> chrIdToName ;

	FILE *fpGTFlist = NULL ;
	FILE *fp = NULL ;
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "--lg" ) )
		{
			fpGTFlist = fopen( argv[i + 1], "r" ) ;
			++i ;
		}
		else if ( !strcmp( argv[i], "-d" ) )
		{
			minAvgDepth = atof( argv[i + 1] ) ;
			++i ;
		}
		/*else if ( !strcmp( argv[i], "-n" ) )
		{
			minSampleCnt = atoi( argv[i + 1] ) ;
			++i ;
		}*/
		else 
		{
			printf( "%s", usage ) ;
			exit( 1 ) ;
		}
	}
	if ( fpGTFlist == NULL )
	{
		printf( "Must use --lg option to speicfy the list of GTF files.\n%s", usage ) ;
		exit( 1 ) ;
	}
	
	std::vector<struct _outputTranscript> transcripts ;
	std::vector<struct _outputTranscript> outputTranscripts ;

	char buffer[4096] ;
	char line[10000] ;
	char chrom[50], tool[20], type[40], strand[3] ;
	//char tid[50] ;
	int start, end ;
	std::vector<struct _pair32> tmpExons ;
	int sampleCnt = 0 ;
	int gid = -1 ;
	int tid = -1 ;
	int chrId = -1 ;
	int chrIdUsed = 0 ;
	char cStrand = '.' ;
	char prefix[50] ;
	double FPKM = 0, TPM = 0, cov = 0 ;

	while ( fgets( buffer, sizeof( buffer ), fpGTFlist ) != NULL )
	{
		int len = strlen( buffer ) ;	
		if ( buffer[len - 1] == '\n' )
		{
			buffer[len - 1] = '\0' ; 
			--len ;
		}
		fp = fopen( buffer, "r" ) ;


		while ( fgets( line, sizeof( line ), fp ) != NULL )
		{
			if ( line[0] == '#' )
				continue ;

			sscanf( line, "%s %s %s %d %d %s %s", chrom, tool, type, &start, &end, buffer, strand ) ;

			if ( strcmp( type, "exon" ) )
			{
				if ( tmpExons.size() > 0 )
				{
					struct _outputTranscript nt ;
					nt.sampleId = sampleCnt ; 
					nt.chrId = chrId ;
					nt.geneId = gid ;
					nt.transcriptId = tid ;
					nt.strand = cStrand ;
					nt.ecnt = tmpExons.size() ;
					nt.exons = new struct _pair32[nt.ecnt] ;
					nt.FPKM = FPKM ;
					nt.TPM = TPM ;
					nt.cov = cov ;
					for ( i = 0 ; i < nt.ecnt ; ++i )
						nt.exons[i] = tmpExons[i] ;
					tmpExons.clear() ;
					transcripts.push_back( nt ) ;
				}
				continue ;
			}

			if ( chrNameToId.find( std::string( chrom ) ) == chrNameToId.end() )
			{
				chrId = chrIdUsed ;

				std::string s( chrom ) ;
				chrNameToId[s] = chrIdUsed ;
				chrIdToName[ chrIdUsed ] = s ;
				++chrIdUsed ;
			}
			else
				chrId = chrNameToId[ std::string( chrom ) ] ;

			GetGTFField( line, "gene_id", buffer ) ;
			gid = GetTailNumber( buffer ) ;
			
			GetGTFField( line, "transcript_id", buffer ) ;
			tid = GetTailNumber( buffer ) ;

			GetGTFField( line, "FPKM", buffer ) ;
			FPKM = atof( buffer ) ;
			
			GetGTFField( line, "TPM", buffer ) ;
			TPM = atof( buffer ) ;

			GetGTFField( line, "cov", buffer ) ;
			cov = atof( buffer ) ;

			cStrand = strand[0] ;

			// Look for the user-defined prefix.
			int len = strlen( buffer ) ;
			j = 0 ;
			for ( i = len - 1 ; i >= 0 ; --i )
			{
				if ( buffer[i] == '.' )
				{
					++j ;
					if ( j >= 2 )
						break ;
				}
			}

			for ( j = 0 ; j < i ; ++i )
			{
				prefix[j] = buffer[j] ;
			}
			if ( i > 0 )
			{
				prefix[j] = '.' ;
				prefix[j + 1] = '\0' ;
			}
			else
				prefix[0] = '\0' ;

			struct _pair32 ne ;
			ne.a = start ;
			ne.b = end ;
			tmpExons.push_back( ne ) ;
		}
		if ( tmpExons.size() > 0 )
		{
			struct _outputTranscript nt ;
			nt.sampleId = sampleCnt ; 
			nt.chrId = chrId ;
			nt.geneId = gid ;
			nt.transcriptId = tid ;
			nt.strand = cStrand ;
			nt.ecnt = tmpExons.size() ;
			nt.exons = new struct _pair32[nt.ecnt] ;
			nt.FPKM = FPKM ;
			nt.TPM = TPM ;
			nt.cov = cov ;
			for ( i = 0 ; i < nt.ecnt ; ++i )
				nt.exons[i] = tmpExons[i] ;
			tmpExons.clear() ;
			transcripts.push_back( nt ) ;
		}

		++sampleCnt ;
		fclose( fp ) ;
	}
	fclose( fpGTFlist ) ;

	// Coalesce same transcripts
	std::sort( transcripts.begin(), transcripts.end(), MultiThreadOutputTranscript::CompSortTranscripts ) ;
	std::vector<int> sampleSupport ;
	int size = transcripts.size() ;
	if ( minSampleCnt > sampleCnt )
		minSampleCnt = sampleCnt ;

	for ( i = 0 ; i < size ; )
	{
		int l ;
		for ( j = i + 1 ; j < size ; ++j )
		{
			if ( MultiThreadOutputTranscript::CompTranscripts( transcripts[i], transcripts[j] ) )
				break ;
		}
		// [i,j) are the same transcripts.
		/*if ( j - i >= fraction * sampleCnt && j - i >= minSampleCnt )
		{
			outputTranscripts.push_back( transcripts[i] ) ;
		}*/

		double sumFPKM = 0 ;
		double sumTPM = 0 ;
		double sumCov = 0 ;
		for ( l = i ; l < j ; ++l )
		{
			sumFPKM += transcripts[l].FPKM ;
			sumTPM += transcripts[l].TPM ;
			sumCov += transcripts[l].cov ;
		}

		transcripts[k] = transcripts[i] ;
		for ( l = i + 1 ; l < j ; ++l )
			delete[] transcripts[l].exons ;

		transcripts[k].FPKM = sumFPKM / ( j - i ) ;
		transcripts[k].TPM = sumTPM / ( j - i ) ;
		transcripts[k].cov = sumCov / ( j - i ) ;
		
		if ( ( j - i < int( fraction * sampleCnt ) || j - i <  minSampleCnt ) && sumCov < minAvgDepth * sampleCnt )
		//if ( transcripts[k].score * ( j - i ) < 1 && ( j - i ) < sampleCnt * 0.5 )
		//if ( sumTPM < sampleCnt || sumFPKM < sampleCnt )
		//if ( sumCov < minAvgDepthn * sampleCnt )
		{
			transcripts[k].FPKM = -1 ;
		}
		else
			sampleSupport.push_back( j - i ) ;
		++k ;
		i = j ;
	}
	transcripts.resize( k ) ;	
	
	// The motivation to separate the coalscing and voting is to allow
	// a gene with no txpt passing the vote to pick a representative.
	size = k ;
	for ( i = 0 ; i < size ; )
	{
		int cnt = 0 ;
		if ( transcripts[i].FPKM != -1 )
			++cnt ;

		for ( j = i + 1 ; j < size ; ++j )
		{
			if ( transcripts[i].geneId != transcripts[j].geneId )
				break ;
			if ( transcripts[j].FPKM != -1 )
				++cnt ;
		}
		int l ;
		/*double *freq = new double[cnt] ;
		for ( l = 0 ; l < cnt ; ++l )
			freq[l] = transcripts[l].FPKM / sampleCnt ;
		qsort( freq, cnt, sizeof( double ), CompDouble ) ;

		for ( l = 0 ; l < cnt ; ++l )
			printf( "%lf\n", freq[l] ) ;
		printf( "===========\n" ) ;
		delete[] freq ;*/
		/*if ( cnt == 0 )
		{
			int maxEcnt = -1 ;
			int maxCnt = 0 ;
			int maxtag ;
			for ( l = i ; l < j ; ++l )
				if ( transcripts[l].ecnt > maxEcnt )
					maxEcnt = transcripts[l].ecnt ;
			
			if ( maxEcnt >= 2 )
			{
				int maxSampleSupport = -1 ;
				maxtag = i ;
				for ( l = i ; l < j ; ++l )
				{
					if ( transcripts[l].ecnt > 1 && transcripts[l].ecnt >= 2  )
					{
						if ( sampleSupport[l] > maxSampleSupport )
						{
							maxSampleSupport = sampleSupport[l] ;
							maxtag = l ;				
							maxCnt = 1 ;
						}
						else if ( sampleSupport[l] == maxSampleSupport )
							++maxCnt ;
					}
				}
				
				if ( maxSampleSupport >= 3 && maxSampleSupport >= ( fraction / 10 ) * sampleCnt && transcripts[maxtag].TPM > 0)
				{
					if ( maxCnt > 1 )
					{
						int maxTPM = -1 ;
						for ( l = i ; l < j ; ++l )
						{
							if ( sampleSupport[l] != maxSampleSupport )
								continue ;
							if ( transcripts[l].TPM > maxTPM )
							{
								maxTPM = transcripts[l].TPM ;
								maxtag = l ;
							}
						}

					}
					transcripts[maxtag].FPKM = 0 ;
					//fprintf( stderr, "recovered: %d %d\n", sampleSupport[ maxtag ], transcripts[ maxtag ].ecnt ) ;
					outputTranscripts.push_back( transcripts[maxtag] ) ;
				}
			}
		}
		else*/
		{
			for ( l = i ; l < j ; ++l )
			{
				//if ( transcripts[l].FPKM >= fraction * sampleCnt && transcripts[l].FPKM >= minSampleCnt )
				if ( transcripts[l].FPKM != -1 )
				{
					outputTranscripts.push_back( transcripts[l] ) ;
				}
			}
		}
		
		i = j ;
	}
	
	// Output
	size = outputTranscripts.size() ;
	//printf( "%d\n", size ) ;
	int transcriptId = 0 ;
	int prevGid = -1 ;
	for ( i = 0 ; i < size ; ++i )
	{
		struct _outputTranscript &t = outputTranscripts[i] ;
		const char *chrom = chrIdToName[t.chrId].c_str() ;
		/*if ( t.geneId != prevGid )
			transcriptId = 0 ;
		else
			++transcriptId ;*/
		transcriptId = outputTranscripts[i].transcriptId ;
		fprintf( stdout, "%s\tPsiCLASS\ttranscript\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; transcript_id \"%s%s.%d.%d\"; FPKM \"%.6lf\"; TPM \"%.6lf\"; cov \"%.6lf\"; sample_cnt \"%d\";\n",
				chrom, t.exons[0].a, t.exons[t.ecnt - 1].b, t.strand,
				prefix, chrom, t.geneId,
				prefix, chrom, t.geneId, transcriptId, t.FPKM, t.TPM, t.cov, sampleSupport[i] ) ;
		for ( j = 0 ; j < t.ecnt ; ++j )
		{
			fprintf( stdout, "%s\tPsiCLASS\texon\t%d\t%d\t1000\t%c\t.\tgene_id \"%s%s.%d\"; "
					"transcript_id \"%s%s.%d.%d\"; exon_number \"%d\"; FPKM \"%.6lf\"; TPM \"%.6lf\"; cov \"%.6lf\"; sample_cnt \"%d\";\n",
					chrom, t.exons[j].a, t.exons[j].b, t.strand,
					prefix, chrom, t.geneId,
					prefix, chrom, t.geneId, transcriptId,
					j + 1, t.FPKM, t.TPM, t.cov,
					sampleSupport[i] ) ;
		}

		prevGid = t.geneId ;
	}
	return 0 ;
}
