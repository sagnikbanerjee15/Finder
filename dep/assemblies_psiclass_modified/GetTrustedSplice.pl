#!/bin/perl

use strict ;
use warnings ;
use List::Util qw[min max];

die "usage: a.pl path_to_list_of_splice_file > trusted.splice\n" if ( @ARGV == 0 ) ;

my %spliceSupport ;
my %spliceSampleSupport ;
my %spliceUniqSupport ;
my %spliceSecSupport ;
my %uniqSpliceSites ;
my %spliceSiteSupport ;

my $sampleCnt = 0 ;
open FP1, $ARGV[0] ;
while ( <FP1> )
{
	++$sampleCnt ;
}
close FP1 ;

open FP1, $ARGV[0] ;
while ( <FP1> )
{
	chomp ;
	open FP2, $_ ;
	while ( <FP2> )
	{
		chomp ;
		my $line = $_ ;
		my @cols = split /\s+/, $line ;
		my $key = $cols[0]." ".$cols[1]." ".$cols[2]." ".$cols[4] ;
		if ( $cols[3] <= 0 )
		{
			$cols[3] = 0.1 ;
		}
		elsif ( $cols[3] == 1 && $sampleCnt > 5 )
		{
			$cols[3] = 0.75 ;
		}

		if ( ! defined $spliceSupport{$key} )
		{
			$spliceSupport{ $key } = $cols[3] ;
			$spliceSampleSupport{ $key } = 1 ;
			$spliceUniqSupport{ $key } = $cols[5] ;
			$spliceSecSupport{ $key } = $cols[6] ;
		}
		else
		{
			$spliceSupport{ $key } += $cols[3] ;
			$spliceSampleSupport{ $key } += 1 ;
			$spliceUniqSupport{ $key } += $cols[5] ;
			$spliceSecSupport{ $key } += $cols[6] ;
		}
		
		for ( my $i = 1 ; $i <=2 ; ++$i )
		{
			$key = $cols[0]." ".$cols[$i] ;
			if ( defined $spliceSiteSupport{ $key } )
			{
				$spliceSiteSupport{ $key } += $cols[3] ;
			}
			else
			{
				$spliceSiteSupport{ $key } = $cols[3] ;
			}

			if ( defined $uniqSpliceSites{ $key } && $uniqSpliceSites{ $key } != $cols[2 - $i + 1] )
			{
				$uniqSpliceSites{ $key } = -1 ;
			}
			else
			{
				$uniqSpliceSites{ $key } = $cols[2 - $i + 1] ;
			}
		}
	}
	close FP2 ;
}
close FP1 ;

foreach my $key (keys %spliceSupport)
{
	next if ( $spliceSupport{ $key } / $sampleCnt < 0.5 ) ;
	#next if ( $spliceUniqSupport{$key} / ( $spliceSecSupport{$key} + $spliceUniqSupport{$key} ) < 0.01 ) ;
	next if ( $spliceUniqSupport{$key} <= 2 && ( $spliceSupport{ $key } / $sampleCnt < 1 || $spliceSampleSupport{$key} < min( 2, $sampleCnt ) ) ) ;
	
	my @cols = split /\s+/, $key ;
	my $flag = 0 ;
	#if ( $cols[2] - $cols[1] + 1 >= 10000 )
	#{
	#	$flag = 1 if ( $spliceSupport{ $key } / $sampleCnt < 1 ) ;
	#}
	my $siteSupport = max( $spliceSiteSupport{ $cols[0]." ".$cols[1] }, $spliceSiteSupport{ $cols[0]." ".$cols[2] } ) ;
	
	if ( $spliceSupport{ $key } < $siteSupport / 10.0 )
	{
		#print  $spliceSupport{ $key } / $siteSupport,  " ", -log( $spliceSupport{ $key } / $siteSupport ) / log( 10.0 ), "\n" ;
		#if ( $cols[1] == 73518141 && $cols[2] == 73518206 )
		#{
		#	print "test: ", $spliceSupport{$key}, " $siteSupport ",  -log( $spliceSupport{ $key } / $siteSupport ), "\n";
		#}
		my $needSample = min( -log( $spliceSupport{ $key } / $siteSupport ) / log( 10.0 ) + 1, $sampleCnt ) ;
		next if ( $spliceSampleSupport{ $key } < $needSample ) ;
	}
	
	if ( $cols[2] - $cols[1] + 1 >= 100000 )
	{
		my $needSample = int( ( $cols[2] - $cols[1] + 1 ) / 100000 ) + 1 ;
		$needSample = $sampleCnt if ( $needSample > $sampleCnt ) ;
		$flag = 1 if ( $spliceUniqSupport{$key} / ( $spliceSecSupport{$key} + $spliceUniqSupport{$key} ) < 0.1 
				|| ( $spliceUniqSupport{ $key } / $sampleCnt < 1 ) 
				|| $spliceSampleSupport{ $key } < $needSample ) ;
		next if ( $flag == 1 && $cols[2] - $cols[1] + 1 >= 300000 ) ;
	}
	if ( $flag == 1 && ( ( $uniqSpliceSites{ $cols[0]." ".$cols[1] } == -1 || $uniqSpliceSites{ $cols[0]." ".$cols[2] } == -1 ) 
		|| $spliceSampleSupport{ $key } <= 1 ) )
	{
		next ;
	}
	print $cols[0], " ", $cols[1], " ", $cols[2], " 10 ", $cols[3], " 10 0 0 0\n" ;
}
