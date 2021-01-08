# this scripts pick random entries from 5ss pos and 3ss pos to make false splice sites, so there are three inputs
#   5ss frag 
#   3ss frag
#   true sites (excluding these sites)
#   strand p/n

use strict;

my $strand = "+";
$strand = "-" if ($ARGV[3] eq "n");


my $max_intron = 500000;

my $s5in_filename = $ARGV[0];
my $s3in_filename = $ARGV[1];
my $truesites_filename =$ARGV[2];

open(s5in, $s5in_filename);
my %s5;

while(my $line =<s5in> )
{
    chomp($line);
    my @a = split("\t", $line);
    die "check strand\n" if($a[2] ne $strand);
    $s5{$a[0]}{$a[1]} =1;
}
close(s5in);

open(truesites, $truesites_filename);
my %truesites_hash;
while(my $line =<truesites> )
{
    chomp($line);
    my ($chr, $start, $end, $name, $score, $strand_here) = split("\t", $line);
    next if($strand_here ne $strand );
    $truesites_hash{$chr}{$start + 15}{$end -15} = 1;	#TODO
}
close(truesites);

open(s3in, $s3in_filename);

while(my $line =<s3in>)
{
    #chomp($line);
    my ($chr, $poss3) = split("\t", $line);
    
    foreach my $poss5 (keys %{$s5{$chr}})
    {
	my $intron_size = $poss3-$poss5;
	if($intron_size < $max_intron and $intron_size> 50 and not exists $truesites_hash{$chr}{$poss5}{$poss3})
	{
	    print join("\t",$chr, $poss5, $poss3, ".",0, $strand ),"\n";
	}
    }
}

