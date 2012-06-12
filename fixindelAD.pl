#!/usr/local/biotools/perl/5.10.0/bin/perl

## Purpose: 
## Fix the AD field from the merged somatic variants from multiple callers 
## 
## 1. Find which format field is AD and which one is SC
## 2.   The first number in AD, the reference reads, should be:
##          the sum of the first two numbers in SC
## 3.   The second number in AD, the alt reads, should be:
##          the sum of the last two numbers in SC


use strict;

## check argument count, print usage statement
my $nargs = @ARGV;
if($nargs !=2){
  printf("usage: fixindelAD.pl  <indel.vcf>  <indelfixed.vcf>\n");
  exit(1);
}

my $indel = @ARGV[0];
my $fixindel = @ARGV[1];


open(INDEL, $indel) || die "indel vcf file could not be opened\n";
open(FIX, ">$fixindel") || die "could not open output vcf\n";

while(<INDEL>) {
    my $l=$_;
    ## print header lines, stop at CHROM header line
    if( $l =~ m/^\x23+/) {
	print FIX $l;	
    } else {
        # split the row, if chrom and position match, add the flag
	my @row = split(' ', $l);
        my @format = split(":", @row[8]);
        my $idxAD=0;
        my $idxSC=0;
        for(my $k=0; $k< ($#format+1); $k++) {
	    if(@format[$k] eq "AD") { 
		$idxAD=$k;
	    }
	    if(@format[$k] eq "SC") { 
		$idxSC=$k;
	    }
	}
	
	## printf("indices: %s, %s\n", $idxAD, $idxSC);

     ## Split into Blood genotype and Blood AD/SC
      
	my @bloodGT = split(":", @row[9]);
	my @bloodSC = split(",", @bloodGT[$idxSC]);
	my @bloodAD = split(",", @bloodGT[$idxAD]);
	## make AD the sum of components of SC
	@bloodAD[0] = $bloodSC[0]+$bloodSC[1];
	@bloodAD[1] = $bloodSC[2]+$bloodSC[3];

      ## Split into Blood genotype and Blood AD/SC
	my @tissueGT = split(":", @row[10]);
	my @tissueSC = split(",", @tissueGT[$idxSC]);
	my @tissueAD = split(",", @tissueGT[$idxAD]);
	@tissueAD[0] = $tissueSC[0]+$tissueSC[1];
	@tissueAD[1] = $tissueSC[2]+$tissueSC[3];

      ## PUT BACK TOGETHER  Blood GT, then Tissue GT
	my $joinbloodAD=join(",", @bloodAD);
	@bloodGT[$idxAD] = $joinbloodAD;
	my $joinbloodGT = join(":", @bloodGT);
	@row[9] = $joinbloodGT;

  ## check steps 
  #	print "bloodAD $joinbloodAD \n";
  #	print "bloodGT $joinbloodGT \n";

	my $jointissueAD=join(",", @tissueAD);
	@tissueGT[$idxAD] = $jointissueAD;
	my $jointissueGT = join(":", @tissueGT);
	@row[10] = $jointissueGT;

	my $joinrow = join("\t", @row); 
	print FIX "$joinrow\n";
    } # end if/else
}  ## end while

