use strict;
# this script scan fa format files and give the number of ACGT in the background


# also gives GT/AG(+)  CT/AC (-) sites in a prob of 1/1000 (to avoid too many samples)
my $prob =1000; 
my %background;
my $chrlistfilename = $ARGV[0];
my $outprefix = $ARGV[1];

open(chrlistfile, $chrlistfilename);

my @chromList = <chrlistfile>;

close(chrlistfile);

chomp(@chromList);

my $len_chromList = @chromList;
my %chr_hash;
foreach my $one (@chromList){
     if($one =~/\/([^\/]*?)\.*fa$/i){
            $chr_hash{$1} = $one;
        }
}

open(f_gt, ">".$outprefix.".GT");
open(f_ct, ">".$outprefix.".CT");
open(fag, ">".$outprefix.".AG");
open(fac, ">".$outprefix.".AC");

foreach my $chrom ( keys(%chr_hash))
{
    my $chromFastaFile = $chr_hash{$chrom};
    open(fin, "$chromFastaFile") or die "can't open the chrom file : $!";
    local ($/) = undef;
    my $contigSeqStr = <fin>;
    close (fin);
    $contigSeqStr =~s/^\>.*?\n//g;
    $contigSeqStr =~s/\s|\n//g;
    my $len_contigSeqStr = length $contigSeqStr;

    for(my $i = 0; $i< $len_contigSeqStr; $i++)
    {
	$background{uc (substr($contigSeqStr, $i, 1)) }++;
	my $j = $i+2;
	my $dstr = uc (substr($contigSeqStr, $i, 2));
	print f_gt join("\t",$chrom, $i, "+"),"\n" if($dstr eq "GT" and int(rand($prob))==500);
	print f_ct join("\t",$chrom, $i, "-"),"\n" if($dstr eq "CT" and int(rand($prob))==500);
	print fag join("\t",$chrom, $j, "+"),"\n" if($dstr eq "AG" and int(rand($prob))==500);
	print fac join("\t",$chrom, $j, "-"),"\n" if($dstr eq "AC" and int(rand($prob))==500);
    }
}

close(f_gt);
close(f_ct);
close(fag);
close(fac);

open(fback, ">".$outprefix.".bg");
foreach my $nt (keys %background)
{
    print fback $nt,"\t", $background{$nt},"\n";
}
close(fback);
