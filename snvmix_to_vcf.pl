#!/usr/bin/perl

# Converts raw svnmix output to VCF 4.1 format
# Jared Evans  evans.jared@mayo.edu
# 1/31/2012

use strict;
use warnings;
use Getopt::Long;


my($input, $output, $gzipped, $sample, $filter_depth, $filter_prob, $help);
GetOptions("in|i=s"	=> \$input,
	"out|o:s"	=> \$output,
	"gzipped|z"	=> \$gzipped,
	"sample|s:s"	=> \$sample,
	"depth|d:i"	=> \$filter_depth,
	"prob|p:f"	=> \$filter_prob,
	"help|h|?|"	=> \&help);

if(not $input){
	print "Missing input file!\n";
	help();
	exit 1;
}
if($gzipped){
	open IN, "gunzip -c $input |" or die "opening gzipped $input\n";
}else{
	open IN, "<$input" or die "opening $input\n";
}
open OUT, ">$output" or die "opening $output\n" if defined $output;

# ensure there is a sample name
if(not defined $sample){
	$sample = "Sample1";
	if(defined $input){
		# try to extract one from the input file name
		$input =~ s{.*/}{};	# remove path  
		$input =~ s{\..*}{}; # remove extensions
		$sample = $input if $input ne "";
	}
}

# print VCF header
if(defined $output){
	print OUT header($sample);
}else{
	print header($sample);
}

while(<IN>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my @locus = split(":",$line[0]);
	my @stats = split(",",$line[3]);
	my @ref = split(":",$stats[0]);
	my @alt = split(":",$stats[1]);
	my $total_depth = $ref[1] + $alt[1];
	my @bases=qw (A C G T);

	my $genotype = "0/0";
	my $prob = $stats[2]; #hom ref
	if($stats[5] == 2){
		$genotype = "0/1";
		$prob = $stats[3]; #het
	}elsif($stats[5] == 3){
		$genotype = "1/1";
		$prob = $stats[4]; #hom alt
	}

	# check filters
	if(defined $filter_depth){
		next if $total_depth < $filter_depth;
	}
	if(defined $filter_prob){
		next if $prob < $filter_prob;
	}
	
	next if (! grep (/$line[1]/,@bases));
	if ($line[2] eq 'N')	{
		$line[2] = '.';
	}	
	
	if(defined $output){
		print OUT join("\t",$locus[0],$locus[1],".",$line[1],$line[2],".",".");
		print OUT "\tAN=2;DP=".$total_depth.";NS=1\tGT:AD:DP:GQ\t".$genotype.":".$ref[1].",".$alt[1].":".$total_depth.":".phred(1-$prob)."\n";
	}else{
		print join("\t",$locus[0],$locus[1],".",$line[1],$line[2],".",".");
		print "\tAN=2;DP=".$total_depth.";NS=1\tGT:AD:DP:GQ\t".$genotype.":".$ref[1].",".$alt[1].":".$total_depth.":".phred(1-$prob)."\n";
	}
}

close IN;
close OUT if defined $output;

sub phred{
	my $val = $_[0];
	return 101 if $val == 0; # can't take the log of 0
	my $phredQual = -10*(log($val)/log(10));
	return sprintf("%.2f",$phredQual);
}

sub header{
	my $sm = $_[0];
	return '##fileformat=VCFv4.1
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth for This Sample">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'.$sm."\n";
}

sub help{

	print 'DESCRIPTION:
	snvmix_to_vcf.pl converts raw snvmix output to VCF format. The snvmix variant call
	probability is log transformed to Phred scale genotype quality (GQ).

USAGE:
	snvmix_to_vcf.pl -i sample.snvs.raw.snvmix.gz -o sample.vcf -z

OPTIONS:
	--in,-i		Path to snvmix file. Required parameter. Input can be gzipped as long 
			as the -z flag is also used. Input files should not have a header line.

	--out,-o 	Optional path to VCF output file. If no file is defined then output will
			be sent to STDOUT. 

	--gzipped,-z 	A flag used if the input file is gzipped.

	--sample,-s 	An optional sample name can be defined. If no sample name is defined 
			then the first part of the input filename will be used as the sample
			name.

	--depth,-d	Optional parameter to filter variants by read depth. A variant will be
			skipped if the total read depth at the position of the variant is less 
			than this value.

	--prob,-p	Optional parameter to filter variants by snvmix probability. A variant will
			be skipped if the snvmix probability is less than this value. Acceptable
			range: 0-1.0

	--help,-h,-?	Display this documentation.

';


}

