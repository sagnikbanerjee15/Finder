open(input, $ARGV[0]);
$totalnum = 0;
while($line = <input>)
{
    next if (substr($line, 0, 1) eq ">" or $line =~/N/);
    chomp($line); 
    for( $i=0; $i<60; $i++)
    {
	$matrix{substr($line, $i, 1)}[$i]++;
    }
    $totalnum++;
}
close(input);

#print "The following is for the c++ matrix\n";
foreach $nt ("A", "C", "G", "T")
{
#    print "{";
    for( $i=0; $i<60; $i++)
    {
	print $matrix{$nt}[$i]/$totalnum,",";
    }
    print "\n";
#    print "}\n";
}
#print "Total: $totalnum\n";

for ($i=0; $i<60; $i++)
{
    print "BS\t";
    foreach $nt ("A", "C", "G", "T")
    {
	print $matrix{$nt}[$i]/$totalnum,"\t";
    }
    print "\n";
}
