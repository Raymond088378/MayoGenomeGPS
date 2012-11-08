#!/usr/local/biotools/perl/5.14.2/bin/perl -w

use strict;

## Add kinome targets by gene symbol
## 1) Read gene names into hash table
## 2) Write field to variants in a kinase gene
my $nargs = @ARGV;
if($nargs !=3){
  printf("usage: addKinome.pl  <kinome.tab> <treatfile.tab> <newtreat.tab>\n");
  exit(1);
}

my $kinome = $ARGV[0];
my $treat = $ARGV[1];
my $newtreat = $ARGV[2];


open(KINOME, $kinome) || die "Gene Kinase file could not be opened\n";
$_ = <KINOME>; ## skip header line
my $gene="nogene";
my %genekinome = ();  # hash table to store gene : "drug(inter);drug;drug(inter)"
my $firstgene = 1;
my $kinomestr;
my @drow=();
while(<KINOME>){
    my $l=$_;
    my @drow = split("\t",$l);

    $gene=$drow[5];	
    $kinomestr = $drow[7];
    $kinomestr =~ s/\"//g; 
    ## write to hash the kinomestr
    $genekinome{$gene} = $kinomestr;
 #   printf "gene %s; kinome %s \n", $gene, $genekinome{$gene};
    
}
# write last gene entry  
$genekinome{$gene} = $kinomestr;
 

close(KINOME);


open(TREAT, $treat)  || die "Treat file could not be opened\n";
open(NEWTREAT, ">$newtreat")  || die "New Treat file could not be opened\n";
my @row;
my $nf=0;
my $geneidx = 0;
my $annotkinome=" ";
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
	print NEWTREAT "$l\tKinase_2007\n";
	next;
    }

    ## all other lines
    @row=split("\t", $l);
    $annotkinome=" ";
    if(exists($genekinome{$row[$geneidx]})) {
	$annotkinome = $genekinome{$row[$geneidx]};
    }
    printf NEWTREAT  "%s\t%s\n", $l, $annotkinome;
    $nf += 1;
    
}

print "treat file has $nf lines \n";
close(TREAT);
close(NEWTREAT);

