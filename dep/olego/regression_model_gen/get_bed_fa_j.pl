# Adapted from Chenghai Xue's script

$starttime=time();

$input_file_1 = $ARGV[0];	# exon junction file
$input_file_2 = $ARGV[1];	# genome file list
$output_file_1 = $ARGV[2];	# exon junction bed (might be less than input_file_1
$output_file_2 = $ARGV[3];	# exon junction fa
#$leftLen = $ARGV[4];
#$rightLen = $ARGV[5];

open(IN_1, "$input_file_1") or die "can't open the input file : $!";
open(IN_2, "$input_file_2") or die "can't open the input file : $!";
open OUT_1, ">$output_file_1" or die "Can not open output_file : $!";
open OUT_2, ">$output_file_2" or die "Can not open output_file : $!";

@chromList = (<IN_2>);
chomp(@chromList);
$len_chromList = @chromList;
print "BED2FA: in $input_file_2, found $len_chromList chromosomes\n";
foreach $one (@chromList){
	if($one =~ /\/([^\/]*?)\.*fa$/i){
		$chr_hash{$1} = $one;
		#print $1,"\n";	
	}
}
@key_chr_hash = keys(%chr_hash);
$len_key_chr_hash = @key_chr_hash;
#@sort_key_chr_hash = sort_chromNo(@key_chr_hash);
#$len_sort_key_chr_hash = @sort_key_chr_hash;
#for($i=0; $i<$len_sort_key_chr_hash; $i++){
#	print "$sort_key_chr_hash[$i]	$chr_hash{$sort_key_chr_hash[$i]}\n";
#}

$num_1=0;
$num_2=0;
$num_count_chrom=0;
my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts);
$current_chrom = "";
while(<IN_1>){
	$num_1++;
	$line = $_;
	chomp $line;
	#print $line,"\n";
	@cols = split ("\t", $line);
        if(scalar(@cols)==12)
	{
	($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts) = @cols;
	}
	if(scalar(@cols)!=12)
	{
		($chrom, $chromStart, $chromEnd, $name, $score, $strand)=@cols;
		$thickStart=$chromStart;
	#	print $thickStart,"\n";
		$thickEnd = $chromEnd;
		$blockCount=1;
		$blockSizes=$chromEnd-$chromStart;
		$blockStarts = 0;	
	}
	$strand="+" if !$strand;
	@a_blockSizes = split (/\,/, $blockSizes);
	@a_blockStarts = split (/\,/, $blockStarts);
	if($chrom ne $current_chrom){
		if($num_1 != 1){
			print "$num_chr_1	$num_chr_2	$len_contigSeqStr\n";			
		}
		print "BED2FA: $chrom:	";
		
		$num_chr_1=0;
		$num_chr_2=0;		

		if(exists $chr_hash{$chrom}){
			$num_count_chrom++;
			$current_chrom = $chrom;
			#print $current_chrom,"\n";
#=pod			
			$chromFastaFile = $chr_hash{$chrom};
			#print $chromFastaFile,"\n";
			open($fin, "<$chromFastaFile") or die "can't open the chrom file : $!";
			local ($/) = undef;
			$contigSeqStr = <$fin>;
			close ($fin);
			#print $contigSeqStr,"mark\t";
			$contigSeqStr =~s/^\>.*?\n//g;
                        #print $contigSeqStr,"mark2\t";

			$contigSeqStr =~s/\s|\n//g;
                        #print $contigSeqStr,"mark3\n";

			$len_contigSeqStr = length $contigSeqStr;
#=cut
		}
		else{
			$num_chr_1++;
			next;
		}
	}
	$num_chr_1++;
	
# modify from here................................
	my @Starts;
	my @Ends;
	my @JuncSeq;
	my $ssStrTag=1;
	for($i_wuj=0;$i_wuj<$blockCount;$i_wuj++)
	{
		$Starts[$i_wuj] = $chromStart + $a_blockStarts[$i_wuj];
		$Ends[$i_wuj] = $Starts[$i_wuj] + $a_blockSizes[$i_wuj];
		$JuncSeq[$i_wuj] = uc substr ($contigSeqStr,$Starts[$i_wuj], $a_blockSizes[$i_wuj]);
		if($strand eq "-"){
		  #$JuncSeq[$i_wuj] = uc string_reverse_complement(lc $JuncSeq[$i_wuj]);
		  $JuncSeq[$i_wuj] = uc (reverse $JuncSeq[$i_wuj]);
		    $JuncSeq[$i_wuj] =~ tr/ACGT/TGCA/;
		    
		}
	}	 
       # for($i_wuj=0;$i_wuj<$blockCount-1;$i_wuj++)
#	{
#	        $ssStr = uc substr ($contigSeqStr, $Ends[$i_wuj], 2) . substr ($contigSeqStr, $Starts[$i_wuj+1]  - 2, 2);
#	        if($strand eq "-"){
 #               $ssStr = uc string_reverse_complement(lc $ssStr);
                #$ssStr = $rc_ssStr;
#	        }
#		$ssStrTag = 0 if ($ssStr ne "GTAG");
		
 #       }
#	if($ssStrTag ==1){
        if(1){
		$num_2++;
		$num_chr_2++;
		print OUT_1 "$line\n";
		#print OUT_2 ">$name\|$chrom\|$chromStart\|$chromEnd\|$strand\|$ssStr\|$num_2\n$junctionSeqStrLeft$junctionSeqStrRight\n";
#		print OUT_2 ">$name\|$chrom\|$chromStart\|$chromEnd\|$strand\|GTAG\|$num_2\|$blockCount\n";
		print OUT_2 ">$name\|$chrom\|$chromStart\|$chromEnd\|$strand\|$num_2\|$blockCount\n";
		if($strand eq "+")
		{
			for($i_wuj=0;$i_wuj<$blockCount;$i_wuj++)
			{
				print OUT_2 $JuncSeq[$i_wuj]; 
			}
		}
		else
		{
			for($i_wuj=$blockCount-1;$i_wuj>-1;$i_wuj--)
                        {
                                print OUT_2 $JuncSeq[$i_wuj];
                        }

		}
		print OUT_2 "\n";
	}
	
}
print "$num_chr_1	$num_chr_2	$len_contigSeqStr\n";
print "BED2FA: in file1, $num_count_chrom chroms, $num_1 beds, $num_2 saved.\n";

close IN_1 or die "can't close the input file : $!";
close IN_2 or die "can't close the input file : $!";
close OUT_1 or die "can't close the output file : $!"; 
close OUT_2 or die "can't close the output file : $!"; 

#######################################
$complete_time = time()-$starttime;
print "BED2FA: Run $complete_time seconds...Done!\n";

#######################################
# sub fuctions
sub string_reverse_complement{
	local($string) = @_;
	local($len_str, $ret, $i, $char);
	
	$len_str = length $string;
	$ret = "";
	for($i=0; $i<$len_str; $i++){
		$char = substr($string, $i, 1);
		if($char eq 'a'){
			$char = 't';	
		}	
		elsif($char eq 't'){
			$char = 'a';	
		}
		elsif($char eq 'c'){
			$char = 'g';	
		}
		elsif($char eq 'g'){
			$char = 'c';	
		}
		else{
			$char = 'n';
		}
		$ret = $char.$ret;		
	}
	
	return $ret;
}

sub sort_chromNo{
	local(@chrom) = @_;
	local($len_key_chr_hash, $i, @sort_chr_hash);
	local(@digit_random, @words_random, @digit_other_1, @digit_other_2, @words_other_1, @words_other_2, @digit, @words);
	local(@sort_digit, @sort_words, @sort_digit_random, @sort_words_random, @sort_digit_other, @sort_words_other);
	local($len_digit, $len_words, $len_digit_random, $len_words_random, $len_digit_other, $len_words_other, $term);
	
	$len_key_chr_hash = @chrom;
	# sort via chr number for printing result
	for($i=0; $i<$len_key_chr_hash; $i++){
		if($key_chr_hash[$i] =~ /chr(\d+)\_random/){
			push(@digit_random, $1);
		}
		elsif($key_chr_hash[$i] =~ /chr(\w+)\_random/){
			push(@words_random, $1);
		}
		elsif($key_chr_hash[$i] =~ /chr(\d+)\_([\w\d\_]+)/){
			push(@digit_other_1, $1);
			push(@digit_other_2, $2);
		}
		elsif($key_chr_hash[$i] =~ /chr(\w+)\_([\w\d\_]+)/){
			push(@words_other_1, $1);
			push(@words_other_2, $2);
		}
		elsif($key_chr_hash[$i] =~ /chr(\d+)/){
			push(@digit, $1);
		}
		elsif($key_chr_hash[$i] =~ /chr(\w+)/){
			push(@words, $1);
		}
		else{
			print "BED2FA: There is unknown type of chromosomes: $key_chr_hash[$i]\n";
		}
	}
	@sort_digit = sort by_mostly_numeric @digit;
	@sort_words = sort by_mostly_string @words;
	@sort_digit_random = sort by_mostly_numeric @digit_random;
	@sort_words_random = sort by_mostly_string @words_random;
	@sort_digit_other = sort_2_array_number_string(\@digit_other_1, \@digit_other_2);
	@sort_words_other = sort_2_array_string_string(\@words_other_1, \@words_other_2);
	
	$len_digit = @sort_digit;
	for($i=0; $i<$len_digit; $i++){
		$term = "chr".$sort_digit[$i];
		push(@sort_chr_hash, $term);
	}
	$len_words = @sort_words;
	for($i=0; $i<$len_words; $i++){
		$term = "chr".$sort_words[$i];
		push(@sort_chr_hash, $term);
	}
	$len_digit_random = @sort_digit_random;
	for($i=0; $i<$len_digit_random; $i++){
		$term = "chr".$sort_digit_random[$i]."_random";
		push(@sort_chr_hash, $term);
	}
	$len_words_random = @sort_words_random;
	for($i=0; $i<$len_words_random; $i++){
		$term = "chr".$sort_words_random[$i]."_random";
		push(@sort_chr_hash, $term);
	}	
	$len_digit_other = @sort_digit_other;
	for($i=0; $i<$len_digit_other; $i=$i+2){
		$term = "chr".$sort_digit_other[$i]."_".$sort_digit_other[$i+1];
		push(@sort_chr_hash, $term);
	}
	$len_words_other = @sort_words_other;
	for($i=0; $i<$len_words_other; $i=$i+2){
		$term = "chr".$sort_words_other[$i]."_".$sort_words_other[$i+1];
		push(@sort_chr_hash, $term);
	}	
	
	return @sort_chr_hash;
}

sub sort_2_array_number_string{
	local($a, $b) = @_;
	local($len_a, $len_b, $i, %family, $one, $two);
	local(@ret);
	
	$len_a = @$a;
	$len_b = @$b;
	if($len_a == $len_b){
		for($i=0; $i<$len_a; $i++){
			$family{$$a[$i]}{$$b[$i]} = 0;									
		}
		for $one (sort by_mostly_numeric keys %family) {
			for $two (sort by_mostly_string keys %{ $family{$one} }) {
					push(@ret, $one);
					push(@ret, $two);
			}	
		}
	}
	else{		
		print "ERROR: Sort array is not same size\n";
		print "a $len_a, b $len_b\n";
	}	
	
	return @ret;
}

sub sort_2_array_string_string{
	local($a, $b) = @_;
	local($len_a, $len_b, $i, %family, $one, $two);
	local(@ret);
	
	$len_a = @$a;
	$len_b = @$b;
	if($len_a == $len_b){
		for($i=0; $i<$len_a; $i++){
			$family{$$a[$i]}{$$b[$i]} = 0;									
		}
		for $one (sort by_mostly_string keys %family) {
			for $two (sort by_mostly_string keys %{ $family{$one} }) {
					push(@ret, $one);
					push(@ret, $two);
			}	
		}
	}
	else{		
		print "ERROR: Sort array is not same size\n";
		print "a $len_a, b $len_b\n";
	}	
	
	return @ret;
}

sub by_mostly_numeric{
#	( $a <=> $b ) || ( $a cmp $b );
	( $a <=> $b );
}

sub by_mostly_string{
#	( $a <=> $b ) || ( $a cmp $b );
	( $a cmp $b );
}

