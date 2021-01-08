# Author: wuj@cshl.edu

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
my $PROG = $0;
my $PROG_ABS_PATH = dirname(Cwd::abs_path($PROG));

######################

my @programs = ('R','sort', 'ls');
foreach my $program (@programs)
{
       die ("CHECK: $program not found\n") if(system("hash $program >/dev/null"));
}


my $annofilename = "";
my $genomedir = "";
my $outputprefix = "userdefined";

GetOptions (
        "g:s"=>\$genomedir,
        "a:s"=>\$annofilename,
        "o:s"=>\$outputprefix
);

my $InputParaDes="      Usage of the script:
	    $PROG   [options]
        -g      genome fasta file location
        -a      annotation file (bed format)
        -o      output prefix
";

if($genomedir eq "" or $annofilename eq "")
{
        print $InputParaDes;
        exit;
}

$genomedir = Cwd::abs_path($genomedir);
$annofilename = Cwd::abs_path($annofilename);

my $cachefolder = $outputprefix.".cache";

if (! -e $cachefolder)
{
        mkdir $cachefolder or die "could not create cache folder $cachefolder\n";
}


print "Extracting junctions from $annofilename\n";

system("perl $PROG_ABS_PATH/bed2ss.pl $annofilename $cachefolder/true.ss.bed");

print "Sorting $cachefolder/true.ss.bed\n";

system("sort -k1,1 $cachefolder/true.ss.bed >$cachefolder/true.ss.sort.bed");

system("ls $genomedir/*.fa >$cachefolder/chr.list");

print "Extracting sequences from $cachefolder/true.ss.sort.bed\n";

system("perl $PROG_ABS_PATH/get_bed_fa_j.pl $cachefolder/true.ss.sort.bed $cachefolder/chr.list $cachefolder/out.bed $cachefolder/true.ss.fasta");

print "Extracting matrix from $cachefolder/true.ss.fasta\n";

system("perl $PROG_ABS_PATH/get_motif_matrix.pl $cachefolder/true.ss.fasta >$cachefolder/true.ss.mat");

print "Extracting introns from $annofilename\n";

system("perl $PROG_ABS_PATH/bed2intron.pl $annofilename $cachefolder/true.intron.bed");

print "Sorting $cachefolder/true.intron.bed\n";

system("sort -k1,1 $cachefolder/true.intron.bed |uniq >$cachefolder/true.intron.sort.bed");

print "Extracting sequences from $cachefolder/true.intron.sort.bed\n";

system("perl $PROG_ABS_PATH/get_bed_fa_j.pl $cachefolder/true.intron.sort.bed $cachefolder/chr.list $cachefolder/out.bed $cachefolder/true.intron.fasta");

print "Extrancting backgroud from $cachefolder/true.intron.fasta\n";

system("perl $PROG_ABS_PATH/get_background.pl $cachefolder/true.intron.fasta > $cachefolder/true.intron.bg");

print "Extracting intron sizes and splice site scores from true dataset\n";

system("perl $PROG_ABS_PATH/getscore_intronsize.pl $cachefolder/true.ss.fasta $cachefolder/true.ss.mat $cachefolder/true.intron.bg >$cachefolder/true.intronsize_score");


if (! -e $cachefolder."/frag_tmp")
{
      mkdir $cachefolder."/frag_tmp" or die "could not create cache folder $cachefolder/frag_tmp\n";
}
print "Extracting GTAG sites from genome... \n";

system("perl $PROG_ABS_PATH/scanfa.pl $cachefolder/chr.list $cachefolder/frag_tmp/frag");

print "Combining random GTAG sites as false set...\n";

system("perl  $PROG_ABS_PATH/combinerandfrags.pl $cachefolder/frag_tmp/frag.GT $cachefolder/frag_tmp/frag.AG $cachefolder/true.ss.sort.bed p >$cachefolder/false.junc");

system("perl $PROG_ABS_PATH/combinerandfrags.pl $cachefolder/frag_tmp/frag.CT $cachefolder/frag_tmp/frag.AC $cachefolder/true.ss.sort.bed n >>$cachefolder/false.junc");

system("perl $PROG_ABS_PATH/randomselectfalse.pl $cachefolder/true.ss.sort.bed $cachefolder/false.junc >$cachefolder/false.samenum.junc");

system("perl $PROG_ABS_PATH/junc2splicesites.pl $cachefolder/false.samenum.junc >$cachefolder/false.samenum.ss.bed");
system("sort -k1,1 -k2,2n $cachefolder/false.samenum.ss.bed >$cachefolder/false.samenum.ss.sort.bed");

print "Extracting sequences from flase splice sites...\n";

system("perl $PROG_ABS_PATH/get_bed_fa_j.pl $cachefolder/false.samenum.ss.sort.bed $cachefolder/chr.list $cachefolder/out.bed $cachefolder/false.samenum.ss.sort.fa");
system("perl $PROG_ABS_PATH/getscore_intronsize.pl $cachefolder/false.samenum.ss.sort.fa $cachefolder/true.ss.mat $cachefolder/true.intron.bg >$cachefolder/false.intronsize_score");

print "Generating file for regression..\n";

open(logitfile, ">$cachefolder/logit.csv");
print logitfile "class,size,score\n";

open(input, "$cachefolder/true.intronsize_score");

while(my $line = <input>)
{
    my @a = split("\t", $line);
    print logitfile join(",", 1, $a[0], $a[1]);
    
}

close(input);

open(input, "$cachefolder/false.intronsize_score");

while(my $line = <input>)
{
    my @a = split("\t", $line);
    print logitfile join(",", 0, $a[0], $a[1]);

}

close(input);


close(logitfile);

print "Using R to regress the Logistic model....\n";

system("R --slave --args $cachefolder/logit.csv $cachefolder/logit.rout <$PROG_ABS_PATH/logit.R");

print "Finishing and generating the model file...\n";

open(output,">$cachefolder/$outputprefix.cfg");
print output "# MATRIX\n";
open(input, "$cachefolder/true.ss.mat");
for (my $i = 0; $i<4; $i++)
{
    my $line = <input>;
    print output $line;
}

close(input);

open(input, "$cachefolder/logit.rout");
while(my $line=<input>)
{
    if($line =~ /^Coefficients/)
    {
	$line = <input>;
	$line = <input>;

	my @a = split(/\s+/, $line );
	print output "# LOGIT_A\n";
	print output $a[1],"\n";
	
	$line = <input>;
	
	my @a = split(/\s+/, $line );
	print output "# LOGIT_B\n";
	print output $a[1],"\n";
	
	$line = <input>;
	my @a = split(/\s+/, $line );
	
	print output "# LOGIT_C\n";
	print output $a[1],"\n";
	
	last;	
   }
}

close(input);
print output "# FOREGROUND_NUM\n1\n";

my %background;
open(input, "$cachefolder/true.intron.bg");
while(my $line = <input>)
{
    chomp($line);
    my @a = split(/\s+/, $line );
    $background{$a[0]} = $a[1];
    
}
close(input);

print output "# BACKGROUND\n";

foreach my $nt ('A', 'C', 'G', 'T')
{
    print output $background{$nt},",";
}
print output "\n";

close(output);

print "Done! Model file generated : $cachefolder/$outputprefix.cfg\n";
