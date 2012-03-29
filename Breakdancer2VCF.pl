### tool to convert break dancer output to vcf output 
    ##  contact : Saurabh Baheti
    ##	email	: baheti.saurabh@mayo.edu
    ##	date	: Feb 2012
	## used this as a refernece http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41

use strict;
use warnings;
use Getopt::Std;

our($opt_i, $opt_f, $opt_o, $opt_s,$opt_t);
print "RAW parameters: @ARGV\n";
getopt('ifost');
if( (!defined $opt_i) && (!defined $opt_f) && (!defined $opt_o) && (!defined $opt_s) && (!defined $opt_t)){
    die ("Usage: $0 \n\t -i [ input BD file ] \n\t -f [ FASTA reference file ] \n\t -o [ output vcf file ] \n\t -s [ sample name to be placed in vcf ] \n\t -t [ path to samtools] \n NOTE: This script works only when reference genome and indexed genome is in the same folder\n");
}
else	{
	my $ref=$opt_f;
	my $infile=$opt_i;
	my $outfile=$opt_o;
	my $sample=$opt_s;
    my $samtools=$opt_t;
	open IN, "<$infile" or die "can not read BD input file $infile :$! \n";
	open OUT, ">$outfile" or die "can not write VCF file $outfile :$! \n";
	open FAIL, ">${outfile}.fail" or die "can not write failed calls file $outfile :$! \n";
	my $old_fh = select(OUT); $| = 1; select($old_fh);
	&print_header($sample,$ref);
	while(my $l = <IN>)	{
		next if($l =~ /^#/);
		chomp $l;
		## header
		## Chr1   Pos1    Orientation1    Chr2    Pos2    Orientation2    Type    Size    Score   num_Reads       
		## chr1    6015088 	2+2-    chr1    6015288 	2+2-    INS     -102    43      2
		## assigning values to the required columns
		## Columns 1-3 and 4-6 are used to specify the coordinates of the two SV breakpoints. The orientation is a string that records the number of reads mapped to the plus (+) or the minus (-) strand in the anchoring regions.
		## Column 7 is the type of SV detected: DEL (deletions), INS (insertion), INV (inversion), ITX (intra-chromosomal translocation), CTX (inter-chromosomal translocation), and Unknown.
		## Column 8 is the size of the SV in bp.  It is meaningless for inter-chromosomal translocations.
		## Column 9 is the confidence score associated with the prediction.
		## Column 10 Total number of supporting read pairs	
		my ($BP1_chr,$BP1_pos,$left,$BP2_chr,$BP2_pos,$right,$type,$size,$score,$num_reads) = split (/\t/, $l);
        my $base1=getFASTABaseSamtools($BP1_chr,$BP1_pos,$ref,$samtools);
        my $base2=getFASTABaseSamtools($BP2_chr,$BP2_pos,$ref,$samtools);
		
        
        ## using chr pos to get the base from the FASTA file using subroutine";
		#my $row=GetChrPos($ref,$BP1_chr);
		#my $base1=GetBaseFasta($ref,$row,$BP1_pos-1,1);
		#$row=GetChrPos($ref,$BP2_chr);
		#my $base2=GetBaseFasta($ref,$row,$BP2_pos-1,1);
		if ( (length($base1) == 1) && (length($base2) == 1) )	{
			## 1st is for + strand reads
			## 2nd is for - strand reads
			my @left_reads=split(/[+-]/,$left);
			my @right_reads=split(/[+-]/,$right);
			## conditions for different types of structural variants
			my $orientation_left;
			my $orientation_right;
			if ( $type =~ /CTX/ ){
				## to figure the direction of translocation
				if($left_reads[0] > $left_reads[1]){
					$orientation_left="+";
				}
				else {
					$orientation_left="-";
				}
				if($right_reads[0] > $right_reads[1]){
					$orientation_right="+";
				}
				else {
					$orientation_right="-";
				}
				### to print to vcf according to the direction found above
				###
				
				
				my ($Alt_left,$INFO_Left,$Alt_right,$INFO_Right);
				if (($orientation_left eq "+")&&($orientation_right eq "+")){
					$Alt_left=join('',$base1,"[",$BP2_chr,":",$BP2_pos,"["); #BND_W
					$Alt_right=join('',$base2,"]",$BP1_chr,":",$BP1_pos,"]"); #BND_Y
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,$score,"PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,$score,"PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}
				elsif (($orientation_left eq "+")&&($orientation_right eq "-")){
					$Alt_left=join('',$base1,"[",$BP2_chr,":",$BP2_pos,"[");	### BND_U
					$Alt_right=join('',"]",$BP1_chr,":",$BP1_pos,"]",$base2);	### BND_V
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,$score,"PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,$score,"PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
					
				}
				elsif (($orientation_left eq "-")&&($orientation_right eq "+")){
					$Alt_left=join('',"]",$BP2_chr,":",$BP2_pos,"]",$base1); ##BND_U
					$Alt_right=join('',$base2,"[",$BP1_chr,":",$BP1_pos,"["); ##BND_V
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,$score,"PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,$score,"PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}
				elsif (($orientation_left=="-")&&($orientation_right=="-")){
					$Alt_left=join('',"[",$BP2_chr,":",$BP2_pos,"[",$base1); #BND_X
					$Alt_right=join('',"[",$BP1_chr,":",$BP1_pos,"[",$base2); #BND_Z
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_",,"CTX_",$.,"_2","DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1","DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,$score,"PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,$score,"PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}	
			}
			if ( $type =~ /INS/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<INS>",$score,"PASS","IMPRECISE;SVTYPE=INS;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			if ( $type =~ /DEL/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<DEL>",$score,"PASS","IMPRECISE;SVTYPE=DEL;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			if ( $type =~ /INV/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<INV>",$score,"PASS","IMPRECISE;SVTYPE=INV;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			if ( $type =~ /ITX/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<ITX>",$score,"PASS","IMPRECISE;SVTYPE=DEL;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}	
		}
		else	{
			print FAIL "$l\n";
		}
	}
	close IN;
	close OUT;
}



## to get chr pos from fai file
sub getFASTABaseSamtools
{
    my $chrID = shift;
    my $basepos = shift;
    my $fastapath = shift;
    my $samtools = shift;
    my @result = `$samtools/samtools faidx $fastapath $chrID:$basepos-$basepos`;
    chomp($result[1]);
    return uc($result[1]);
} 

sub GetChrPos	{
	my $fai=$_[0].".fai";
	open FAI, "$fai" or die "";
	my $r;
	while(<FAI>){
		my ($chr)=split(/\t/,$_);
		if ($chr eq "$_[1]"){
			$r="$.";
		}
	}
	close FAI;
	return ($r);
	
}
## to get base from the fasta file only one base
sub GetBaseFasta	{
	local $/ = "\n>";
	open FASTA, "$_[0]" or die "";
	my $line=$_[1];
	my $pos=$_[2];
	my $num_bases=$_[3];
	my $base;
	while(<FASTA>){
		if ($. == $line)	{
			my $l=uc($_);
			$l =~ s/^>*.+\n//;  ## remove fasta header
			$l =~ s/\n//g;  # remove endlines
			$base=substr($l,$pos,$num_bases);
			last;
		}		
	}
	close FASTA;
	undef $\;
	return ($base);
}
## to print the date
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}
## to print the header for the vcf file
sub print_header{
my $date=&spGetCurDateTime();
my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=Breakdancer2VCF.pl
##reference=$_[1]
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 132">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=SC,Number=1,Type=Integer,Description="Somatic Score">
##FORMAT=<ID=CQ,Number=1,Type=Integer,Description="Consensus Quality">
##FORMAT=<ID=SNVQ,Number=1,Type=Integer,Description="SNV Quality">
##FORMAT=<ID=RMS,Number=1,Type=Integer,Description="RMS mapping Quality">
##FORMAT=<ID=BRQ,Number=1,Type=Integer,Description="Mean base Quality of reads supporting reference">
##FORMAT=<ID=MRQ,Number=1,Type=Integer,Description="Mean mapping Quality of reads supporting reference">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Depth of reads supporting reference">
##FORMAT=<ID=BVQ,Number=1,Type=Integer,Description="Mean base Quality of reads supporting variant">
##FORMAT=<ID=MVQ,Number=1,Type=Integer,Description="Mean mapping Quality of reads supporting variant">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of reads supporting variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$_[0]\n};
print OUT $header;
}


