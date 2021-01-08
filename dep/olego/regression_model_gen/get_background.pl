# input: fa file
# print ACGT percent in the sequences
use strict;

my $inputfilename = $ARGV[0];

my %background;

$background{'A'} = 0;
$background{'C'} = 0;
$background{'T'} = 0;
$background{'G'} = 0;

my @keys = keys %background;

my $totalnum = 0;

open(input, $inputfilename);
while(my $line = <input>)
{
    next if ( substr($line, 0, 1) eq ">" );
# chomp($line);
#    for( my $i=0; $i<length($line); $i++)
#	{
#	    $background{uc (substr($line, $i, 1)) }++;
#	    $totalnum ++;
#	}
    chomp($line);
    $totalnum = $totalnum + length($line) ;
   foreach my $key (@keys)
   {
       my $count = ($line =~ s/($key)//g);
      # print $key,"\t", $count,"\n";
       $background{$key} = $background{$key}+$count;
   }
   
	
}
close(input);
foreach my $nt (keys %background)
{
    print $nt,"\t", $background{$nt}/$totalnum,"\n";
}
