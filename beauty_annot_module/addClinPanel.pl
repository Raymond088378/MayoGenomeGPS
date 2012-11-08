#!/usr/local/biotools/perl/5.14.2/bin/perl -w

use strict;

## Add flags for variants in genes on known clinical panels by gene symbol
## 1) Read gene names into hash table
## 2) Write field to variants in these panel genes
my $nargs = @ARGV;
if($nargs !=3){
  printf("usage: addClinPanel.pl  <panel.txt> <treatfile.tab> <newtreat.tab>\n");
  exit(1);
}

my $panel = $ARGV[0];
my $treat = $ARGV[1];
my $newtreat = $ARGV[2];


open(PANEL, $panel) || die "Clinical panel gene file could not be opened\n";
my $gene="nogene";
my %genepanel = ();  # hash table to store gene : "drug(inter);drug;drug(inter)"
my $firstgene = 1;
while(<PANEL>){
    my $l=$_;
    chomp $l;
    $gene=$l;	
   
    ## write to hash
    $genepanel{$gene} = "YES";
 #   print "gene $gene; panel: $genepanel{$gene}";
    
}
 
close(PANEL);


open(TREAT, $treat)  || die "Treat file could not be opened\n";
open(NEWTREAT, ">$newtreat")  || die "New Treat file could not be opened\n";
my @row;
my $nf=0;
my $geneidx = 0;
my $annotpanel=" ";
while(<TREAT>) {
    my $l = $_;  # chomp \n from end of $_ 
    chomp $l;
    if($l =~ /SNPEFF/){
	## first header, repeat in out file with empty field at end
	print NEWTREAT "$l\t\n";
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
	print NEWTREAT "$l\tClinical_Panel_Gene\n";
	next;
    }

    ## all other lines
    @row=split("\t", $l);
    $annotpanel=" ";
    if(exists($genepanel{$row[$geneidx]})) {
	$annotpanel = $genepanel{$row[$geneidx]};
    }
    printf NEWTREAT  "%s\t%s\n", $l, $annotpanel;
    $nf += 1;
    
}

print "treat file has $nf lines \n";
close(TREAT);
close(NEWTREAT);

