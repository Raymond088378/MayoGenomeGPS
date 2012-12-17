#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;


my $in=shift @ARGV;
my $poly=shift @ARGV;
my $output=shift @ARGV;
## chr4:367199         Q15928                  benign

open FH, "$in" or die " cannot open $in :$!\n";
open ANNOT, "$poly" or die " can not open $poly : $!\n";
open OUT, ">$output" or die " ca not open $output : $!\n";
my $in_l=<FH>;chomp $in_l;
print OUT "$in_l\t\t\n";

$in_l=<FH>;chomp $in_l;
print OUT "$in_l\tUniprotID\tpolyphen2\n";

$in_l=<FH>;
my $an_l=<ANNOT>;
my $old_fh = select(OUT); $| = 1; select($old_fh);

while(defined $in_l || defined $an_l)	{
	chomp $in_l if defined $in_l;
	chomp $an_l if defined $an_l;
	if(!defined $an_l){
		print OUT "$in_l\t-\tnone\n";
		$in_l=<FH>;	
	}
	elsif (!defined $in_l)	{
		last;	
	}
	else    {
		my @input=split(/\t/,$in_l);
		my @annotation=split(/\t/,$an_l);
		my ($chr,$pos)=split(/:/,$annotation[0]);
		$annotation[1] =~ s/^\s+//;
		$annotation[2] =~ s/^\s+//;
		$annotation[1] =~ s/\s+$//;
		$annotation[2] =~ s/\s+$//;
		if ($input[1] == $pos)	{
			print OUT "$in_l\t$annotation[1]\t$annotation[2]\n";
			$in_l=<FH>;
			$an_l=<ANNOT>;
		}
		elsif ( $input[1] < $pos)	{
			print OUT "$in_l\t-\tnone\n";
			$in_l=<FH>;
		}
		else	{
			$an_l=<ANNOT>;
		}
	}
}
close FH;
close ANNOT;
close OUT;

			
