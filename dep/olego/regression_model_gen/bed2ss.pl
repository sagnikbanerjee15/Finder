# convert bed file into splice sites ( +-15 nt)
use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;

my $prog = basename ($0);

my $verbose = 0;

my $upflanking = 15;
my $dnflanking = 15;

GetOptions (
    "v|verbose"=>\$verbose,
);

if (@ARGV != 2)
{
    print "bed format to junctions\n";
    print "usage: $prog <in.bed> <out.bed> \n";
    exit(1);
}

my ($inputfilename, $outputfilename) = @ARGV;

my ($fin, $fout);

my %junchash;

open($fin, $inputfilename) or Carp::croak "cannot open file $inputfilename to read!\n";
while(my $line = <$fin>)
{
    chomp($line);
    my @a =split("\t", $line);
    my @blockSizes = split(",", $a[10]);
    my @blockStarts = split(",", $a[11]);

    for (my $i=0; $i<$a[9]-1; $i++)
    {
	    my $start = $a[1] + $blockStarts[$i] + $blockSizes[$i];
	    my $end = ($a[1] + $blockStarts[$i+1]);
    	$junchash{join(",", $a[0], $start, $end, $a[5])}++;
	    
    }

}
close($fin);

open($fout, ">$outputfilename") or Carp::croak "cannot open file $outputfilename to write!\n";
foreach my $key (keys %junchash)
{
    my @a = split(",", $key);
    print $fout $a[0],"\t";
    print $fout $a[1]-$upflanking, "\t";
    print $fout $a[2]+$dnflanking, "\t";
    print $fout "\.\t";
    print $fout $junchash {$key},"\t";
    print $fout $a[3],"\t";
    print $fout $a[1]-$upflanking, "\t";
    print $fout $a[2]+$dnflanking, "\t";
    print $fout "255,255,0\t", "2\t", "30,30,\t", "0,", $a[2]-$a[1], "\n"; 
}
close($fout);
