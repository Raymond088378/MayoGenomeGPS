#!/usr/local/biotools/perl/5.14.2/bin/perl -w

use strict;

## Create an Excel Spreadsheet with three worksheets
## 3) "complete_annotation": a full set of all columns, filtered to not 
##     have only minimal "junk" variants"
## 2) "key_annotation" columns for the same variants in 3) 
## 1) "high-quality" same columns as "key", but only the highest-quality variants, which may be empty, especially for indels.
## Add comments from a data-dictionary file for the headers, such that
## the comment appears when the mouse hovers over the header cell.

my $nargs = @ARGV;
if($nargs !=3){
  printf("usage: addCancerDrugs.pl  <cancerdrug.tab> <treatfile.tab> <newtreat.tab>\n");
  exit(1);
}

my $drug = $ARGV[0];
my $treat = $ARGV[1];
my $newtreat = $ARGV[2];


open(DRUG, $drug) || die "Cancer Drug file could not be opened\n";
$_ = <DRUG>; ## skip header line
my $gene="nogene";
my %genedrug = ();  # hash table to store gene : "drug(inter);drug;drug(inter)"
my $firstgene = 1;
my $drugstr;
my @drow=();
while(<DRUG>){
    my $l=$_;
    my @drow = split("\t",$l);
    ## check if new gene
    if($gene ne $drow[0]) {
	if($firstgene>0) {
	    # on first gene, so don't write to genedrug
	    # add to drugstr
	    $firstgene=0;
	   # print "first gene: $drow[0], drug $drow[2]\n";
	    	    
	} else {
	    ## write to hash the drugstr
	    $genedrug{$gene} = $drugstr;
	    #printf "gene %s; drug %s \n", $gene, $genedrug{$gene};
	    $drugstr="";
	}
	# set current gene, start new drugstr
	$gene=$drow[0];	
	if(length($drow[1]) > 1) {
	    $drugstr = $drow[2]."(".$drow[1].")";
	} else {
	    $drugstr = $drow[2];
	}

    } else {
	## on same drug
	## add to drugstr
	if(length($drow[1]) > 1) {
	    $drugstr .= ";".$drow[2]."(".$drow[1].")";
	} else {
	    $drugstr .= ";".$drow[2];
	}
    }
}
# write last gene entry  
$genedrug{$gene} = $drugstr;
 
close(DRUG);

open(TREAT, "<", $treat)  || die "Treat file could not be opened\n";
open(NEWTREAT, ">$newtreat")  || die "New Treat file could not be opened\n";
my @row;
my $nf=0;
my $geneidx = 0;
my $candrug=" ";
while(<TREAT>) {
    my $l = $_;  # chomp \n from end of $_ 
    chomp $l;
    if($l =~ /SNPEFF/){
	## first header, repeat in out file with empty field at end
	print NEWTREAT "$l\t \n";
	next;
    }

    @row=split("\t", $l);
    if($l =~ /geneList/) {
	for(my $j=0; $j<=$#row;$j++) {
	    if($row[$j] eq "geneList") {
		$geneidx = $j;	    
	    }
	}
#	print "geneList index is $geneidx\n";
	print NEWTREAT "$l\tNCCN_Cancer_Drug\n";
	next;
    }

    ## all other lines
    @row=split("\t", $l);
    $candrug=" ";
	#print "$geneidx => $row[$geneidx] --> $genedrug{$row[$geneidx]}\n";
    if( exists( $genedrug{$row[$geneidx]} ) ) {
	$candrug = $genedrug{$row[$geneidx]};
    }
    printf NEWTREAT  "%s\t%s\n", $l, $candrug;
    $nf += 1;
    
}


print "treat file has $nf lines \n";
close(TREAT);
close(NEWTREAT);

