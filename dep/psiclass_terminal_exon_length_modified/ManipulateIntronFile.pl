#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl intronA intronB [op]\n" if (@ARGV == 0 ) ;

my %intronInfo ;
my %chromRank ;

sub sortIntron
{
	my @cols1 = split /\s+/, $a ;
	my @cols2 = split /\s+/, $b ;

	if ( $cols1[0] ne $cols2[0] )
	{
		$chromRank{ $cols1[0] } cmp $chromRank{ $cols2[0] } ;
	}
	elsif ( $cols1[1] != $cols2[1] )
	{
		$cols1[1] <=> $cols2[1] ;
	}
	else
	{
		$cols1[2] <=> $cols2[2] ;
	}
}

open FP1, $ARGV[0] ;
my $cnt = 0 ;
while ( <FP1> )
{
	chomp ;
	my $line = $_ ;
	my @cols = split /\s+/ ;
	push @cols, 1 ;
	@{ $intronInfo{ $cols[0]." ".$cols[1]." ".$cols[2] } } = @cols ;
	if ( !defined $chromRank{ $cols[0]} )
	{
		$chromRank{ $cols[0] } = $cnt ;
		++$cnt ;
	}
}
close FP1 ;

open FP1, $ARGV[1] ;
while ( <FP1> )
{
	chomp ;
	my $line = $_ ;
	my @cols = split /\s+/ ;
	my $key = $cols[0]." ".$cols[1]." ".$cols[2] ;
	next if ( !defined $intronInfo{ $key } ) ;
	my @infoCols = @{ $intronInfo{ $key } } ;
	$infoCols[4] = $cols[4] if ( $infoCols[4] ne "+" || $infoCols[4] ne "-" ) ;
	$infoCols[9] |= 2  ;

	@{ $intronInfo{ $key } } = @infoCols ;
}
close FP1 ;

foreach my $key (sort sortIntron keys %intronInfo )
{
	my @infoCols = @{ $intronInfo{ $key } } ;
	#print join( " ", @infoCols ), "\n" ;

	next if ( $infoCols[9] != 3 ) ;
	pop @infoCols ;
	print join( " ", @infoCols ), "\n" ;
}
