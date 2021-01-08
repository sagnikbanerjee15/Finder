open(input, $ARGV[0]);
while($line =<input>)
{
    chomp($line);
    my ($chr, $start, $end, $name, $score, $strand) =split("\t", $line);
    print join ("\t", $chr, $start-15, $end+15, $name, $score, $strand, $start-15, $end+15, "255,255,0", 2, "30,30,", "0,".($end-$start) ), "\n";
}
close(input);
