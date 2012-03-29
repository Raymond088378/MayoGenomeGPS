### tool to convert break dancer output to vcf output 
    ##  contact : Saurabh Baheti
    ##	email	: baheti.saurabh@mayo.edu
    ##	date	: Feb 2012
	## used this as a refernece http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
use strict;
use warnings;
use Getopt::Std;

our($opt_i, $opt_f, $opt_o, $opt_s, $opt_t);
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
		### header information from crest
		##left_chr, left_pos,
		## left_strand, # of left soft-clipped reads, right_chr, right_pos, right_strand,
		# right soft-clipped reads, SV type, coverage at left_pos, coverage at
		##right_pos, assembled length at left_pos, assembled length at right_pos,
		##average percent identity at left_pos, percent of non-unique mapping reads at
		##left_pos, average percent identity at right_pos, percent of non-unique mapping
		##reads at right_pos, start position of consensus mapping to genome,
		##starting chromosome of consensus mapping, position of the genomic mapping of
		##consensus starting position, end position of consensus mapping to genome,
		##ending chromsome of consnesus mapping, position of genomic mapping of
		##consensus ending posiiton, and consensus sequences.  For inversion(INV), the
		##last 7 fields will be repeated to reflect the fact two different breakpoints
		##are needed to identify an INV event.
	
		#chr16   57602935        +       0       chr1    1116624 +       5       CTX     326     161     0       46      0.985768200349447       0.0764705882352941      0.989214869050972       0.0776119402985075      1       chr16   57602890	100     chr1    1116677 ACCTCCTCCACCTTCATCATCACCATCACCTCCTCCACCTCCATCACCACCTGTCAATCTGCCTTCACCCCAAACATCCATCCAATCCATCCACAATCCC
		my ($BP1_chr,$BP1_pos,$left,$left_reads,$BP2_chr,$BP2_pos,$right,$right_reads,$type) = split (/\t/, $l);
		## using chr pos to get the base from the FASTA file using subroutine";
		my $base1=getFASTABaseSamtools($BP1_chr,$BP1_pos,$ref,$samtools);
        my $base2=getFASTABaseSamtools($BP2_chr,$BP2_pos,$ref,$samtools);
        #my $row=GetChrPos($ref,$BP1_chr);
		#my $base1=GetBaseFasta($ref,$row,$BP1_pos-1,1);
		#$row=GetChrPos($ref,$BP2_chr);
		#my $base2=GetBaseFasta($ref,$row,$BP2_pos-1,1);
		my $num_reads=$left_reads+$right_reads;
		my $size=$BP1_pos - $BP2_pos;$size=abs($size);
		my $orientation_left;
		my $orientation_right;
		if ( (length($base1) == 1) && (length($base2) == 1) )	{
			
			if ( $type =~ /CTX/ ){
				## to figure the direction of translocation
				$orientation_left=$left;
				$orientation_right=$right;	
				my ($Alt_left,$INFO_Left,$Alt_right,$INFO_Right);
				if (($orientation_left eq "+")&&($orientation_right eq "+"))	{	##correct 
					$Alt_left=join('',$base1,"[",$BP2_chr,":",$BP2_pos,"["); #BND_W
					$Alt_right=join('',$base2,"]",$BP1_chr,":",$BP1_pos,"]"); #BND_Y
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,".","PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,".","PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}
				elsif (($orientation_left eq "+")&&($orientation_right eq "-")){
					$Alt_left=join('',$base1,"[",$BP2_chr,":",$BP2_pos,"[");	### BND_U
					$Alt_right=join('',"]",$BP1_chr,":",$BP1_pos,"]",$base2);	### BND_V
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,".","PASS",$INFO_Left,"\t");
                    print OUT join ("\t","GT","1/1","\n");
					print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,".","PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
						
				}
				elsif (($orientation_left eq "-")&&($orientation_right eq "+")){
					$Alt_left=join('',"]",$BP2_chr,":",$BP2_pos,"]",$base1); ##(BND_U)
					$Alt_right=join('',$base2,"[",$BP1_chr,":",$BP1_pos,"["); ##BND_V
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,".","PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,".","PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}
				elsif (($orientation_left=="-")&&($orientation_right=="-")){
					$Alt_left=join('',"[",$BP2_chr,":",$BP2_pos,"[",$base1); #BND_X
					$Alt_right=join('',"[",$BP1_chr,":",$BP1_pos,"[",$base2); #BND_Z
					$INFO_Left=join('',"SVTYPE=BND;MATE_ID=bnd_",,"CTX_",$.,"_2",";DP=$num_reads");
					$INFO_Right=join('',"SVTYPE=BND;MATE_ID=bnd_","CTX_",$.,"_1",";DP=$num_reads");
					print OUT join("\t",$BP1_chr,$BP1_pos,"bnd_CTX_$._l",$base1,$Alt_left,".","PASS",$INFO_Left,"\t");
					print OUT join ("\t","GT","1/1","\n");
                    print OUT join ("\t",$BP2_chr,$BP2_pos,"bnd_CTX_$._2",$base2,$Alt_right,".","PASS",$INFO_Right,"\t");
					print OUT join ("\t","GT","1/1","\n");
				}	
			}		
			elsif ( $type =~ /INS/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<INS>",".","PASS","IMPRECISE;SVTYPE=INS;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			elsif ( $type =~ /DEL/ ){
				
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<DEL>",".","PASS","IMPRECISE;SVTYPE=DEL;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			elsif ( $type =~ /INV/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<INV>",".","PASS","IMPRECISE;SVTYPE=INV;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}
			elsif ( $type =~ /ITX/ ){
				print OUT join("\t",$BP1_chr,$BP1_pos,".",$base1,"<ITX>",".","PASS","IMPRECISE;SVTYPE=ITX;END=$BP2_pos;SVLEN=$size;DP=$num_reads","GT","1/1","\n");
			}	
		}
		else	{
			print FAIL "$l\n";
		}
	}	
	close IN;
	close OUT;
}

############### SUB fucntions#########
####################################
## to get chr pos from fai file
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
##source=CREST2VCF.pl
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