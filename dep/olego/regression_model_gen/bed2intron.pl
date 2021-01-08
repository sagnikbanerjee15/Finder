use strict;
use warnings;

use Carp;

my $verbose = 0;

if (@ARGV != 2)
{
    print "extract introns\n";
    print "usage:  <in.bed> <out.bed> \n";
    exit(1);
}

my ($inputfilename, $outputfilename) = @ARGV;

my ($fin, $fout);


open($fin, $inputfilename) or Carp::croak "cannot open file $inputfilename to read!\n";
open($fout, ">$outputfilename") or Carp::croak "cannot open file $outputfilename to write!\n";

while(my $line = <$fin>)
{
    chomp($line);
    my @a =split("\t", $line);
    my @blockSizes = split(",", $a[10]);
    my @blockStarts = split(",", $a[11]);
    my $chr = $a[0];
    my $start = $a[1];
    my $end = $a[2];
    my $strand = $a[5];
    for (my $i=0; $i<$a[9]-1; $i++)
    {
	    my $blockstart = $start + $blockStarts[$i] + $blockSizes[$i];
	    my $blockend = $start + $blockStarts[$i+1];
	    print $fout join("\t", $chr, $blockstart, $blockend, ".", "0", $strand),"\n";
    }

}
close($fout);
close($fin);

