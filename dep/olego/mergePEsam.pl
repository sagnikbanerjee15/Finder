#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Carp;


my $verbose = 0;

my $maxdist = 5000000; #cross 11 exons

my $samestrand = 0;
my $oppostrand = 1;
my $nostrand   = 0;
my $nocheckinput = 0;

GetOptions (
    "d:i"=>\$maxdist,
    "ss|same-strand"=>\$samestrand,
    "ns|no-strand"=>\$nostrand,
    "nci|no-check-input"=>\$nocheckinput,
    "v|verbose"=>\$verbose
);

die "--ss and --ns cannot be used together.\n" if($nostrand == 1 and $samestrand == 1);
    
$oppostrand = 0 if ($nostrand == 1 or $samestrand == 1);

my $prog = $0;

if ( @ARGV != 3)
{
    print STDERR "Merge sam format output from PE reads\n";
    print STDERR "Usage: $prog [options] <end1.sam> <end2.sam> <out.sam>\n\n";
	print STDERR "  gz files accepted for input; Also, you can use - to direct the output into STDOUT\n";
    print STDERR "	Use --ss or --ns to specify the strategy of using strand information to filter read pairs.\n\n";
    print STDERR "	-d                      max distance between the two ends on genome. [$maxdist]\n";
    print STDERR "	--ss, --same-strand     the two ends come from the same strand, instead of requring opposite strands by default \n";
    print STDERR "	--ns, --no-strand       do not use strand information. \n";
    print STDERR "	--nci, --no-check-input do not check the input file . \n"; 
    print STDERR "	-v			verbose\n";
    exit 1;
}


my ($inSAMFile1,$inSAMFile2, $outFile) = @ARGV;
my ($fin1,$fin2, $fout); 

if ($inSAMFile1 =~/\.gz$/)
{
	open ($fin1, "gunzip -c $inSAMFile1 | ")||Carp::croak "cannot open file $inSAMFile1 to read\n";
}
else
{
	open ($fin1, "<$inSAMFile1") || Carp::croak "cannot open file $inSAMFile1 to read\n";
}

if ($inSAMFile2 =~/\.gz$/)
{
	open ($fin2, "gunzip -c $inSAMFile2 | ")||Carp::croak "cannot open file $inSAMFile2 to read\n";
}
else
{
	open ($fin2, "<$inSAMFile2") || Carp::croak "cannot open file $inSAMFile2 to read\n";
}

if($outFile ne "-")
{
    open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
}
else
{
    $fout = *STDOUT;
}

my $linenum = 0;

while (my $line1 = <$fin1>)
{
    chomp $line1;
    my $line2 = <$fin2>;
    $linenum ++;
    print STDERR "$linenum reads done ...\n" if $verbose && $linenum % 10000 == 0;
    chomp $line2;
    if ($line1=~/^\@/)
    {
	print $fout $line1,"\n";
	next;
    }
    #get information from file 1
    my ($QNAME1, $FLAG1, $RNAME1, $POS1, $MAPQ1, $CIGAR1, $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, $TAG1) = split (/\s+/, $line1, 12);

    my ($QNAME2, $FLAG2, $RNAME2, $POS2, $MAPQ2, $CIGAR2, $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, $TAG2) = split (/\s+/, $line2, 12);
    # put the infor into the arrays 
    die "In input $linenum, read names not matched" if (( !$nocheckinput) and ( substr($QNAME1, 0, length($QNAME1)-2) ne substr($QNAME2, 0, length($QNAME2)-2) ) );
    my (@chr1, @chr2);
    my (@pos1, @pos2);
    my (@size1, @size2);
    my (@tag1, @tag2);
    my (@cigar1, @cigar2);
    my (@nm1, @nm2);
    my (@strand1, @strand2);
    my (@sensestrand1, @sensestrand2); # the strand of the RNA

    push (@chr1, $RNAME1);
    push (@chr2, $RNAME2);
    
    push (@pos1, $POS1);
    push (@pos2, $POS2);

    push (@cigar1, $CIGAR1);
    push (@cigar2, $CIGAR2);

#   modify FLAGs

    $FLAG1 = $FLAG1 | 0x0001 ; #the read is paired
    $FLAG2 = $FLAG2 | 0x0001 ;
    $FLAG1 = $FLAG1 | 0x0040 ; # the read is the first read in the pair
    $FLAG2 = $FLAG2 | 0x0080 ; # the read is the second read in the pair
    
    my $flaginfo1 = decodeSAMFlag ($FLAG1);
    my $flaginfo2 = decodeSAMFlag ($FLAG2);
 
    # do not look at strand now
    # unless ($TAG1 and $TAG2)
    # only need to look at those pairs both of which have tags
    if ( ($FLAG1&0x0004)==0x0004 or ($FLAG2&0x0004)==0x0004 )
    # if at least one end unmapped, simply output the original lines
    {
	# update MRNM and MPOS
	$MRNM2 = $RNAME1;
	$MRNM1 = $RNAME2;
	
	$MPOS1 = $POS2;
	$MPOS2 = $POS1;
	
	# TLEN is already 0, no need to update 
	
	# update the flags
	# 0x2 0, not properly paired
	# 0x8 the other end not mapped
	 
	$FLAG1 = $FLAG1 | 0x0008 if ( ($FLAG2&0x0004) == 0x0004);
	$FLAG2 = $FLAG2 | 0x0008 if ( ($FLAG1&0x0004) == 0x0004);
	
	# 0x20, strand infor
	$FLAG1 = $FLAG1 | 0x0020 if ( ($FLAG2&0x0010) == 0x0010);
	$FLAG2 = $FLAG2 | 0x0020 if ( ($FLAG1&0x0010) == 0x0010);
	
	# update some other columns
	
	if ( $TAG1)
	{		
		print $fout join("\t",$QNAME1, $FLAG1, $RNAME1, $POS1, $MAPQ1, $CIGAR1, $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, $TAG1),"\n" ;
	}
	else
	{
		print $fout join("\t",$QNAME1, $FLAG1, $RNAME1, $POS1, $MAPQ1, $CIGAR1, $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1),"\n";
	}

	if ( $TAG2)
	{
		print $fout join("\t",$QNAME2, $FLAG2, $RNAME2, $POS2, $MAPQ2, $CIGAR2, $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, $TAG2), "\n";
	}
	else
	{
                print $fout join("\t",$QNAME2, $FLAG2, $RNAME2, $POS2, $MAPQ2, $CIGAR2, $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2), "\n";	
	}
	next;
    }
       push (@size1, getChromSize($CIGAR1));
    push (@size2, getChromSize($CIGAR2));
 
    push (@strand1, $flaginfo1->{'query_strand'});
     push (@strand2, $flaginfo2->{'query_strand'});
     
    push (@tag1, $TAG1); 

    if( $TAG1=~/NM\:\S*\:(\d+)/)
    {
	push(@nm1, $1);
    }
    else
    {
	push(@nm1, -1);
    }
    
    push (@tag2, $TAG2);

    if( $TAG2=~/NM\:\S*\:(\d+)/)
    {
	push(@nm2, $1);
    }
    else
    {
	push(@nm2, -1);	
    }
    
    if ( $TAG1=~/XS\:\S*\:([\+\-\.])/)
    {
	push(@sensestrand1, $1);
    }
    else
    {
	die "Error: input 1 $linenum: No XS tag!\n"; 
    }

   if ( $TAG2=~/XS\:\S*\:([\+\-\.])/)
    {
	push(@sensestrand2, $1);
    }
    else
    {
	die "Error: input 2 $linenum: No XS tag!\n"; 
    }

      
    # scan XA tags
    # XA should be the last tag

    if ($TAG1=~/XA\:\S\:(\S*)$/)
    {
	my @XAstrs = split(";",$1);
	for(my $i=0; $i<@XAstrs; $i++)
	{
#	    XA:Z:chr4,+149621574,100M,0;chr2,+80678177,100M,0;
	    if ($XAstrs[$i]=~/^(\w*?),([\+\-])(\d*?),(\w*?),(\d*?),([\+\-\.])$/)
	    {
		push(@chr1, $1);
		push(@strand1, $2);
		push(@pos1, $3);
		push(@cigar1, $4);
		push(@size1, getChromSize($4));
		push(@nm1, $5);
		push(@sensestrand1, $6);
	    }
	    else
	    {
		die "Error: input 1 $linenum: Wrong XA format.\n";
	    }
	}
    }
    if ($TAG2=~/XA\:\S\:(\S*)$/)
      {
	  my @XAstrs = split(";",$1);
	  for(my $i=0; $i<@XAstrs; $i++)
	  {
#	    XA:Z:chr4,+149621574,100M,0;chr2,+80678177,100M,0;
	      if ($XAstrs[$i]=~/^(\w*?),([\+\-])(\d*?),(\w*?),(\d*?),([\+\-\.])$/)
	      {
		  push(@chr2, $1);
		  push(@strand2, $2);
		  push(@pos2, $3);
		  push(@cigar2, $4);
		  push(@size2, getChromSize($4));
		  push(@nm2, $5);
		  push(@sensestrand2, $6);
	      }
	      else
	      {
		  die "Error: input 2 $linenum: Wrong XA format.\n";

		}
	  }
      }
    # scan for the best match
#    my $foundmatches=0, $bestdistance = $distance, $bestnm = 0;
    my %distanceij;
    my %rankij;
    for (my $i = 0; $i<@chr1; $i++)
     {
	 for(my $j = 0; $j<@chr2; $j++)
	{
	    if($chr1[$i] eq $chr2[$j] and abs($pos2[$j]-$pos1[$i])<$maxdist)
	    {
		next if($oppostrand and ($strand1[$i] eq $strand2[$j]) );
		next if($samestrand and ($strand1[$i] ne $strand2[$j]) );
		
		#$bestdistance =  abs($pos2[$j]-$pos2[$i]);
		#$bestnm = $nm1[$i] + $nm2[$j];
		$distanceij{$i.",".$j} = abs($pos2[$j]-$pos1[$i]);
		$rankij{$i.",".$j} = $i+$j;
		#print join("\t", $i,$j),"\n";	
	    }
	}
     }
     
     if(scalar (keys %distanceij) >0 )
    {
	# only output reliable hits if matches found
	# but some entries could be output several times. 
	my $ontop = 1;
	my $outputline1;
	my $outputline2;
	
	foreach my $ij (sort {$rankij{$a} <=> $rankij{$b}} keys %rankij)
	{
	    my ($i,$j) = split(",", $ij);
	    if($ontop == 1)
	    {
		# update MRNM and MPOS
		$MRNM1 = "=";
		$MRNM2 = "=";
		$MPOS1 = $pos2[$j];
		$MPOS2 = $pos1[$i];
		# update TLEN
		if ($pos1[$i] > $pos2[$j])
		{
		    $ISIZE2 = $pos1[$i] + $size1[$i] - $pos2[$j];
		    $ISIZE1 = 0-$ISIZE2;
		}
		else
		{
		    $ISIZE1 = $pos2[$j] + $size2[$j] - $pos1[$i];
		    $ISIZE2 = 0-$ISIZE1;
		}
		# update FLAGS
		# properly paired
		$FLAG1 = $FLAG1 | 0x0002;
		$FLAG2 = $FLAG2 | 0x0002;
		# strand infor
		if($strand1[$i] eq "+")
		{
		    # set as 0
		    $FLAG1 = $FLAG1 & (~0x0010);
		    $FLAG2 = $FLAG2 & (~0x0020);
		}
		else
		{
		    # set as 1
		    $FLAG1 = $FLAG1 | 0x0010;
		    $FLAG2 = $FLAG2 | 0x0020;
		    
		}
		if($strand2[$j] eq "+")
		{
		    $FLAG2 = $FLAG2 & (~0x0010);
		    $FLAG1 = $FLAG1 & (~0x0020);
		}
		else
		{
		    $FLAG2 = $FLAG2 | 0x0010;
		    $FLAG1 = $FLAG1 | 0x0020;
		}
		if ( $strand1[$i] ne $strand1[0])
		{
			$SEQ1 = reverse $SEQ1;
			$SEQ1 =~ tr/ACGTacgt/TGCAtgca/;
			$QUAL1 = reverse $QUAL1;
		}
		$outputline1 = join("\t", $QNAME1, $FLAG1, $chr1[$i], $pos1[$i], $MAPQ1, $cigar1[$i], $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, "NM:i:".$nm1[$i]);
		$outputline1 = $outputline1."\tXS:A:".$sensestrand1[$i];
		if(scalar (keys %distanceij) >1)
		{
		    $outputline1 = $outputline1."\tXT:A:R";
		    $outputline1 = $outputline1."\tXA:Z:";
		}
		else
		{
		    $outputline1 = $outputline1."\tXT:A:U";
		}
		if ( $strand2[$j] ne $strand2[0])
		{
			$SEQ2 = reverse $SEQ2;
                        $SEQ2 =~ tr/ACGTacgt/TGCAtgca/;
                        $QUAL2 = reverse $QUAL2;
		}
		$outputline2 = join("\t", $QNAME2, $FLAG2, $chr2[$j], $pos2[$j], $MAPQ2, $cigar2[$j], $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, "NM:i:".$nm2[$j]);
		$outputline2 = $outputline2."\tXS:A:".$sensestrand2[$j];
		if(scalar (keys %distanceij) >1)
		{
		    $outputline2 = $outputline2."\tXT:A:R";
		    $outputline2= $outputline2."\tXA:Z:" ;
		}
		else
		{
		    $outputline2 = $outputline2."\tXT:A:U";
		}
		$ontop = 0;
	    }
	    else
	    {
		$outputline1 = $outputline1 .join(",",$chr1[$i], $strand1[$i].$pos1[$i],$cigar1[$i], $nm1[$i], $sensestrand1[$i] ).";";
		$outputline2 = $outputline2 .join(",",$chr2[$j], $strand2[$j].$pos2[$j],$cigar2[$j], $nm2[$j], $sensestrand2[$j] ).";";
	    }
	    
        }
	print $fout $outputline1, "\n";
	print $fout $outputline2, "\n";

    }
    else
    {
	# no match, output the original hits
	# update FLAGS
	# 0x2: not properly paired
	# strand
	$FLAG1 = $FLAG1 | 0x0020 if ( ($FLAG2&0x0010) == 0x0010);
	$FLAG2 = $FLAG2 | 0x0020 if ( ($FLAG1&0x0010) == 0x0010);

	# update the names
	if ($RNAME1 eq $RNAME2)
	{
	    $MRNM2 = "=";
	    $MRNM1 = "=";
	}
	else
	{
	    $MRNM2 = $RNAME1;
	    $MRNM1 = $RNAME2;
	}
		
	$MPOS1 = $POS2;
	$MPOS2 = $POS1;
	if($RNAME1 eq $RNAME2 )
	{
	    if ($POS1 > $POS2)
	    {
		$ISIZE2 = $POS1 + $size1[0] - $POS2;
		$ISIZE1 = 0-$ISIZE2;
	    }
	    else
	    {
		$ISIZE1 = $POS2 + $size2[0] - $POS1;
		$ISIZE2 = 0-$ISIZE1;
	    }
	}

					
	print $fout join("\t",$QNAME1, $FLAG1, $RNAME1, $POS1, $MAPQ1, $CIGAR1, $MRNM1, $MPOS1, $ISIZE1, $SEQ1, $QUAL1, $TAG1),"\n" ;
	print $fout join("\t",$QNAME2, $FLAG2, $RNAME2, $POS2, $MAPQ2, $CIGAR2, $MRNM2, $MPOS2, $ISIZE2, $SEQ2, $QUAL2, $TAG2), "\n";
	
    }
}

print STDERR "Done! Totally $linenum lines processed.\n" if $verbose;

close($fin1);
close($fin2);


sub decodeSAMFlag
{
        my $flag = $_[0];

        #print "flag = $flag\n";
        $flag = sprintf ("%012b", $flag);

        #print "flag binary = $flag\n";
        my @flags = split (//, $flag);

        my $flagInfo = {
                PE=>$flags[11],
                PE_map=>$flags[10],
                query_map=>$flags[9],
                mate_map=>$flags[8],
                query_strand=>$flags[7] == 0 ? '+' : '-',
                mate_strand=>$flags[6] == 0 ? '+' : '-',
                read_1_or_2=> $flags[5] == 1 ? 1 : 2 };
        return $flagInfo;
}

sub getChromSize
{
    my $CIGAR = $_[0];
    if ($CIGAR=~/[^\d+|M|N|I|D]/g)
    {
	Carp::croak "unexpected CIGAR string: $CIGAR in \n";
	
    }
    my $currentSize = 0;
    while ($CIGAR=~/(\d+)([M|N|I|D])/g)
    {
	my ($size, $type) = ($1, $2);
	$currentSize += $size if ($type eq "D" or $type eq "M" or $type eq "N");
    }
    return $currentSize;
    
}
