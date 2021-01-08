open(tpfile,$ARGV[0]);
@tpentries = <tpfile>;
close(tpfile);
$tpnum = @tpentries;



open(fpfile, $ARGV[1]);
@fpentries = <fpfile>;
$fpnum = @fpentries;
my %picked;
for($i=0; $i<$tpnum; $i++)
{
    $randnum = int(rand($fpnum));
    if (not exists $picked{$randnum})
    {
	print $fpentries[$randnum];
	$picked{$randnum}=1;
    }
    else
    {
	redo;
    }
}

close(fpfile);
