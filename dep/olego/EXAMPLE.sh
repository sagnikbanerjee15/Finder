# More example script can be found at
# http://zhanglab.c2b2.columbia.edu/index.php/OLego_Documentation#Examples

Forward_fastq=r1.fq
Reverse_fastq=r2.fq
Model=model/mm.cfg
Junction_db=mm10.intron.hmr.bed
Index=mm10.fa
olego -v -r $Model -j $Junction_db -o r1.sam $Index $Forward_fastq
olego -v -r $Model -j $Junction_db -o r2.sam $Index $Reverse_fastq
mergePEsam.pl r1.sam r2.sam merge.sam
sam2bed.pl merge.sam merge.bed
bed2junc.pl merge.bed merge.junc.bed
olego -v -r $Model -j merge.junc.bed --non-denovo -o r1.remap.sam $Index $Forward_fastq
olego -v -r $Model -j merge.junc.bed --non-denovo -o r2.remap.sam $Index $Reverse_fastq
mergePEsam.pl r1.remap.sam r2.remap.sam merge.remap.sam
