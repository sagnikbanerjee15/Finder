#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;


my $prog = basename ($0);
my $printUniqOnly = 0;
my $anchor = 0;
my $shuffle = 0;

my $verbose = 0;
my $useRNAStrand = 0; # use the strand of the RNA instead of the read

GetOptions (
	"u|uniq"=>\$printUniqOnly,
	"a|anchor:i"=>\$anchor,
	"shuffle"=>\$shuffle,
	"r|use-RNA-strand"=>\$useRNAStrand,
#	"s|separate-bed"=>\$separateBed, 
	"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print STDERR "extract junction reads from SAM and shuffle the order of exonic sequences\n";
	print STDERR "Usage: $prog [options] <in.sam> <out1.bed> [out2.bed]\n";
	print STDERR " <in.sam> 	: gzip compressed input file with .gz extension is allowed\n";
	print STDERR " <out.fastq> 	: output fastq file\n";
	print STDERR " You can also use - to specify STDIN for input or STDOUT for output\n\n";
	print STDERR "options:\n";
	print STDERR "-u,--uniq          :  print uniquely mapped reads only\n";
	print STDERR "-a,--anchor        :  anchor size ($anchor)\n";
#	print STDERR "-r,--use-RNA-strand:  force to use the strand of the RNA based on the XS tag \n";
	print STDERR "-v,--verbose       :  verbose\n";
	exit (1);
}

my ($inSAMFile, $outFastqFile) = @ARGV;


my ($fin, $fout);

if ( $inSAMFile eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($inSAMFile =~/\.gz$/)
	{
		open ($fin, "gunzip -c $inSAMFile | ") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
	else
	{
    	open ($fin, "<$inSAMFile") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
}
if ( $outFastqFile eq "-")
{
     $fout = *STDOUT;
}
else
{
    open ($fout, ">$outFastqFile") || Carp::croak "cannot open file $outFastqFile to write\n";
}


my $i = 0;
my $found = 0;

while (my $line = <$fin>)
{
	chomp $line;

	next if $line=~/^\s*$/;
	next if $line=~/^\@/;

	print STDERR "$i ...\n" if $verbose && $i % 50000 == 0;
	$i++;

	my $sam = lineToSam ($line);
	
	my $CIGAR = $sam->{'CIGAR'};
	next unless $CIGAR =~/\d+M\d+N/;

	my $uniq = 0;
	$uniq = 1 if $sam->{"TAGS"}=~/XT:A:U/;

	if ($printUniqOnly)
	{
		next unless $uniq;
	}
	
	my $QNAME = $sam->{'QNAME'};
	my $QUAL = $sam->{'QUAL'};
	my $SEQ = $sam->{'SEQ'};

	
	#consider only those without indels
	next if $CIGAR=~/[^\d+|M|N]/g;

	my @blockSizes;

	my $currLen = 0;

	my (@seqBlocks, @qualBlocks);
	while ($CIGAR =~/(\d+)M/g)
	{
		push @blockSizes, $1;
		$currLen += $1;
	}
	next if $currLen != length($SEQ);

	#Carp::croak Dumper (\@blockSizes), "\n";
	$currLen = 0;
	
	my $anchorCheckPass = 1;
	foreach my $len (@blockSizes)
	{
		push @seqBlocks, substr ($SEQ, $currLen, $len);
		push @qualBlocks, substr ($QUAL, $currLen, $len);
		$currLen += $len;
		if ($currLen == $len || $currLen == length ($SEQ))
		{
			#first block or last block
			$anchorCheckPass = 0 unless $len >= $anchor;
		}
	}
	next unless $anchorCheckPass == 1;
	
	if ($shuffle)
	{
		@seqBlocks = reverse @seqBlocks;
		@qualBlocks = reverse @qualBlocks;
	}

	$SEQ = join ("", @seqBlocks);
	$QUAL = join ("", @qualBlocks);
	
	print $fout "@", $QNAME, "\n";
	print $fout $SEQ, "\n";
	print $fout "+", "\n";
	print $fout $QUAL, "\n";
}

print STDERR "Done! Totally $i lines processed! \n" if $verbose;

close ($fin) if $inSAMFile ne '-';
close ($fout) if $outFastqFile ne '-';




#subroutines in Align.pm
sub lineToSam
{
	my $line = $_[0];
	my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, $MRNM, $MPOS, $ISIZE, $SEQ, $QUAL, $TAGS) = split (/\s+/, $line, 12);
	return {
	QNAME => $QNAME,
	FLAG=> $FLAG,
	RNAME=>$RNAME,
	POS=>$POS,
	MAPQ=>$MAPQ,
	CIGAR=>$CIGAR,
	MRNM=>$MRNM,
	MPOS=>$MPOS,
	ISIZE=>$ISIZE,
	SEQ=>$SEQ,
	QUAL=>$QUAL,
	TAGS=>$TAGS
	};
}



