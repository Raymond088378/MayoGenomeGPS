#!/usr/bin/perl
use warnings;
use strict;
#use Time::Local;
## modified to generate one more file with tri allelelic snps 10/05/2011

my $usage = "usage: perl ss2vcf.pl [somatic sniper (ss) output file] [output vcf file] [ dbsnp132 single class] [normal sample name] [paired tumor sample name] [outputfile triallelic snps]\n";
my @par = @ARGV;
(@par==6) or die $usage;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my $dbsnp=$ARGV[2];
my $normal_name = $ARGV[3];
my $tumor_name = $ARGV[4];
my $triAllele = $ARGV[5];
my %snphash;

open TRI, ">$triAllele" or die "can't open $triAllele\n";
open DB, "<$dbsnp" or die "can't open $dbsnp\n";
while (my $line = <DB>){
	chomp $line;
	my ($bin,$chr,$start,$end,$rsid,$other) = split("\t",$line);
	$chr =~ s/chr//;
	my $key=$chr.'|'.$end;
	$snphash{$key} = $rsid;	
}
close DB;

open IN, "<$infile" or die "can't open $infile\n";
open OUT, ">$outfile" or die "can't open $outfile\n";
&print_header($normal_name,$tumor_name);

while (my $line2 = <IN>){
	next if ($line2 =~ /^$/);
	chomp $line2;
	my @col = split("\t",$line2); 
	my $id = &get_rsid( $col[0],$col[1] );
#	my $db = '';
#	$db = ';DB' if ($id ne '.');
	
	my ($n_alt,$n_gt,$n_pass) = split(/:/,&IUB2base($col[3],$col[2]));
	my ($t_alt,$t_gt,$t_pass) = split(/:/,&IUB2base($col[4],$col[2]));
	
	if ($n_pass eq 'TRI' || $t_pass eq 'TRI'){
		print TRI "$line2\n";
	}
	else	{
		my $alt;
		if($n_alt ne $col[2]){
			$alt = $n_alt;
		}
		elsif($t_alt ne $col[2]){
			$alt = $t_alt;
		}
		my $total_dp = $col[12]+$col[13];
		my $s_num = 2;
		my $pass = 'PASS';
		if ($n_pass eq 'FAIL' || $t_pass eq 'FAIL'){
			$pass = 'ploidy';
		}
		print OUT "$col[0]\t$col[1]\t$id\t$col[2]\t$alt\t.\t$pass\t";
		print OUT "NS=$s_num;DP=$total_dp";
		print OUT "\t";
		print OUT "GT:SC:CQ:SNVQ:RMS:DP:BRQ:MRQ:DR:BVQ:MVQ:DV\t";
		print OUT "$n_gt:0:$col[9]:$col[10]:$col[11]:$col[13]:$col[20]:$col[21]:$col[22]:$col[23]:$col[24]:$col[25]\t";
		print OUT "$t_gt:$col[5]:$col[6]:$col[7]:$col[8]:$col[12]:$col[14]:$col[15]:$col[16]:$col[17]:$col[18]:$col[19]\n";		
	}
}
close OUT;
close IN;

sub print_header{
my $date=&get_date();
my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=ss2vcf.pl
##reference=file:///RandD/GATK/resources/allchr_hg19.fa
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
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Depth of reads supporting variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$_[0]\t$_[1]\n};
	print OUT $header;
}
sub get_rsid {
	my $chr = $_[0];
	my $pos = $_[1];
	my $k=$chr.'|'.$pos;
	if ( exists $snphash{$k}){
		return $snphash{$k};
	}
	else{
		return ".";
	}
}

#IUB parsing
#M (A or C)
#R (A or G)
#W (A or T)
#S (C or G)
#Y (C or T)
#K (G or T)
#B (C, G or T/U)
#D (A, G or T/U)
#H (A, C or T/U)
#V (A, C or G)
#N (A, C, G or T/U)

sub IUB2base {
	my $c = $_[0];
	my $ref = $_[1];
	if ($c =~ /[AGTC]/){
		if ($c eq $ref){	# homo-ref
			$c.':0/0:PASS';
		}
		else{			# hetero
			$c.':1/1:PASS';
		}
	}
	elsif ($c =~ /[BDHVN]/){
		'.:./.:ploidy';
	}
	elsif ($c eq 'M'){
		my $alt1 = 'A';
		my $alt2 = 'C';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{							#Tri Allelic
			print "warning: M, but neither A nor C\n";
			'.:./.:TRI';
		}
	}
	elsif ($c eq 'R'){
		my $alt1 = 'A';
		my $alt2 = 'G';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{
			print "warning: R, but neither A nor G\n";
			'.:./.:TRI';
		}
	}
	elsif ($c eq 'W'){
		my $alt1 = 'A';
		my $alt2 = 'T';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{
			print "warning: W, but neither A nor T\n";
			'.:./.:TRI';
		}
	}
	elsif ($c eq 'S'){
		my $alt1 = 'C';
		my $alt2 = 'G';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{
			print "warning: S, but neither C nor G\n";
			'.:./.:TRI';
		}
	}
	elsif ($c eq 'Y'){
		my $alt1 = 'C';
		my $alt2 = 'T';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{
			print "warning: Y, but neither C nor T\n";
			'.:./.:TRI';
		}
	}
	elsif ($c eq 'K'){
		my $alt1 = 'G';
		my $alt2 = 'T';
		if ( $alt1 eq $ref){	# hetero
			$alt2.':0/1:PASS';
		}
		elsif ( $alt2 eq $ref){			# hetero
			$alt1.':0/1:PASS';
		}
		else{
			print "warning: K, but neither G nor T\n";
			'.:./.:TRI';
		}
	}	
}

sub get_date{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$year += 1900;
	$mon++;
	$mon = sprintf("%2d", $mon);
	$mon =~ tr/ /0/;
	return "$year$mon$mday";
}
