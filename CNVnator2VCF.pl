#!/usr/local/biotools/perl/5.10.0/bin/perl

### tool to convert cnvnator output to vcf output 
    ##  contact : Saurabh Baheti
    ##	email	: baheti.saurabh@mayo.edu
    ##	date	: work in progress

##  deletion
##  chr10   792601  795400  sample  0.624158        +
##  chr10   802801  803600  sample  0.301309        +

##  duplication
##  chr20   1389801 1554600 sample  2.97556 +
##  chr20   1556601 2541000 sample  2.97814 +


use strict;
use warnings;
use Getopt::Std;

our($opt_i, $opt_f, $opt_o, $opt_s,$opt_t);
print "RAW parameters: @ARGV\n";
getopt('ifost');
if( (!defined $opt_i) && (!defined $opt_f) && (!defined $opt_o) && (!defined $opt_s) ){
    die ("Usage: $0 \n\t -i [ input CNV file ] \n\t -f [ FASTA reference file ] \n\t -o [ output vcf file ] \n\t -s [ sample name to be placed in vcf ] \n\t -t [ path to samtools]\n NOTE: This script works only when reference genome and indexed genome is in the same folder\nOuputs default GT in vcf = 1/1\n");
}
else	{
    my $ref=$opt_f;
    my $infile=$opt_i;
    my $outfile=$opt_o;
    my $sample=$opt_s;
    my $samtools=$opt_t;
    open IN, "<$infile" or die "can not read CNVnator input file $infile :$! \n";
    open OUT, ">$outfile" or die "can not write VCF file $outfile :$! \n";
    open FAIL, ">${outfile}.fail" or die " can not write fail calls :$!\n";
    &print_header($sample,$ref);
    while(my $l = <IN>)	{
	next if($l =~ /^#/);
	chomp $l;
	## chr start stop sample normalized read depth strand
	my ($chr,$start,$stop,$sample,$normalized_depth,$strand)=split(/\t/,$l);
    my $base2=getFASTABaseSamtools($chr,$stop,$ref,$samtools);
        #    my $row=GetChrPos($ref,$chr);
        #my $base2=GetBaseFasta($ref,$row,$stop-1,1);
	if (length($base2) == 1)
	{
	    my $length=$stop-$start;
	    $normalized_depth=sprintf("%.2f",$normalized_depth);
	    if ($normalized_depth < 1){		## deletion
		if ($strand =~ /\+/){
		    print OUT join ("\t",$chr,$start,".",$base2,"<DEL>",".","PASS","IMPRECISE;SVTYPE=DEL;END=$stop;SVLEN=$length","GT:CN","1/1:$normalized_depth\n");
		}
		else{
		    print OUT join ("\t",$chr,$start,".",$base2,"<DEL>",".","PASS","IMPRECISE;SVTYPE=DEL;END=$stop;SVLEN=-$length","GT:CN","1/1:$normalized_depth\n")
		}
	    }
	    else	{		## duplication
		if ($strand =~ /\+/){
		    print OUT join ("\t",$chr,$start,".",$base2,"<DUP>",".","PASS","IMPRECISE;SVTYPE=DUP;END=$stop;SVLEN=$length","GT:CN","1/1:$normalized_depth\n");
		}
		else{
		    print OUT join ("\t",$chr,$start,".",$base2,"<DUP>",".","PASS","IMPRECISE;SVTYPE=DUP;END=$stop;SVLEN=-$length","GT:CN","1/1:$normalized_depth\n");
		}
	    }  
	}
	else{
	    print FAIL "$l\n";
	}
    }	
    close IN;
    close OUT;
    close FAIL;
}

### SUB ROUTINES

sub getFASTABaseSamtools
{
    my $chrID = shift;
    my $basepos = shift;
    my $fastapath = shift;
    my $samtools = shift;
    my @result = `$samtools/samtools faidx $fastapath $chrID:$basepos-$basepos`;
    chomp($result[1]) if defined $result[1];
    return uc($result[1]);
} 
## to get chr pos from fai file
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
##source=CNVnator2VCF.pl
##reference=$_[1]
##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$_[0]\n};
print OUT $header;
}

        
