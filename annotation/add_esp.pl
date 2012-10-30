#!/usr/local/biotools/perl/5.10.0/bin/perl

use warnings;
#use strict;
use Getopt::Std;

our ($opt_i, $opt_d, $opt_o);
#print "Raw parameters: @ARGV\n";
getopt('ido');

unless ( (defined $opt_i) && (defined $opt_d) && (defined $opt_o) ){
	die ("Usage: $0 [-i input file containing chr and pos in tab delimited format -full path] [-d db file] [-o output file path]\n");
}

my $infile = $opt_i;
my $outfile = $opt_o;
my $db_file=$opt_d;
my %db_hash;
my %db_ref;


####################################################
# parse db source file and store info in a hash
####################################################

#`dos2unix $infile`;
open SNP, "<$db_file" or die "can not open file $db_file: $!\n";
my $line = <SNP>;
$line .= <SNP>;
while ( $line = <SNP> ) {
	chomp $line;
	my ($chr_pos,$rsid,$db_v,$allele,$eur_al,$afr_al,$o3,$maf,$other ) = split (/[ ]/, $line);
	my($chr,$pos) = split(/:/,$chr_pos);
	my($eur1,$eur2) = split ("/",$eur_al);
	my($afr1,$afr2) = split ("/",$afr_al);
	my($eur_al_1,$eur_af_1) = split ("=",$eur1);
	my($eur_al_2,$eur_af_2) = split ("=",$eur2);
	my($afr_al_1,$afr_af_1) = split ("=",$afr1);
	my($afr_al_2,$afr_af_2) = split ("=",$afr2);
	my ($eur_minor,$eur_major, $afr_minor,$afr_major);
	if($eur_af_1 > $eur_af_2){
		$eur_minor = $eur_al_2;
		$eur_major = $eur_al_1;
	}
	else{
                $eur_minor = $eur_al_1;
                $eur_major = $eur_al_2;
	}
	if($afr_af_1 > $afr_af_2){
                $afr_minor = $afr_al_2;
                $afr_major = $afr_al_1;
        }
	else{
                $afr_minor = $afr_al_1;
                $afr_major = $afr_al_2;
	}	
	my($eur_maf,$afr_maf,$all) = split ("/",$maf);
	my $eur_af = 100 - $eur_maf;
	my $afr_af = 100 - $afr_maf;
	$eur_maf=$eur_maf/100;
	$afr_maf=$afr_maf/100;
	$eur_af=$eur_af/100;
	$afr_af=$afr_af/100;
	$chr = uc($chr) if $chr eq "x" || $chr eq "y";
	$db_hash{'chr'.$chr}{$pos}{$allele}{'eur'} = $eur_minor.'/'.$eur_major.','.$eur_maf.'/'.$eur_af;
	$db_hash{'chr'.$chr}{$pos}{$allele}{'afr'} = $afr_minor.'/'.$afr_major.','.$afr_maf.'/'.$afr_af;
	
}
close SNP;

# add db allele and maf column to input file
open IN, "<$infile" or die "can't open $infile\n";
open OUT, ">$outfile" or die "can't open $outfile\n";

my $header = <IN>;chomp $header;
print OUT "$header\tESP5400_EUR_maf\tESP5400_AFR_maf\n";
while (my $line2=<IN>){
	chomp $line2;
	my ($chr,$pos) = split(/\t/,$line2);
	if (exists $db_hash{$chr}{$pos} ) {
		my $af;
		for my $al ( keys %{$db_hash{$chr}{$pos}}) {
			$string = join("\t",$line2,$db_hash{$chr}{$pos}{$al}{'eur'},$db_hash{$chr}{$pos}{$al}{'afr'});
		}
	}
	else {
		$string = join("\t",$line2,"-","-");
	}
	print OUT "$string\n";
}
close IN;
close OUT;

exit;
