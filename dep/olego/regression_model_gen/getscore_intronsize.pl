use strict;
#input 1 is a fa file

#>WWW|chr1|4774176|4775663|+|2|2
#GACGTCTTTCCTACAAGAAAAAGCGTTTACCCGTTTTCTG
#GTCCACATGCGTGAGTATTTTTCCACCTAGCCCTGTGATG

#input 2 is matrix
#input 3 is the background
my %matrix;
open(input2, $ARGV[1]);
my $lineid=-1;
while(my $line = <input2>)
{

    chomp($line);
    my @a = split(/\s+/, $line );
    next if $a[0] ne "BS";
    #ACGT
    $lineid++;
    $matrix{$lineid}{'A'} = $a[1];
    $matrix{$lineid}{'C'} = $a[2];
    $matrix{$lineid}{'G'} = $a[3];
    $matrix{$lineid}{'T'} = $a[4];
}

my %background;

#$background{'A'} =0.272;
#$background{'C'} =0.211;
#$background{'G'} =0.217;
#$background{'T'} =0.300;

close(input2);

open(input3, $ARGV[2]);

while(my $line = <input3>)
{
    chomp($line);
    my @a = split(/\s+/, $line );
    $background{$a[0]} = $a[1];
}

close(input3);

open(input1, $ARGV[0]);
while(my $line = <input1>)
{
    chomp($line);
    my @a = split('\|', $line);
    #my $marker = substr($a[0], 1,3);
    my $intronsize = $a[3]-$a[2];

    $line = <input1>;
    chomp($line);
    my $score=0;
    for(my $i=0; $i<length($line); $i++)
    {
	next if (substr($line, $i,1) eq 'N');
	$score = $score+ log( $matrix{$i}{substr($line, $i,1)} /$background{substr($line, $i,1)} );
    }
    print $intronsize,"\t", $score,"\n";
	
}
close(input1);
