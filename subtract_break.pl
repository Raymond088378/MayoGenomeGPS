#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;

open TUMOR, "$ARGV[0]"  or die "";
open NORMAL, "$ARGV[1]" or die "";

my %var;

while(my $l = <NORMAL>)	{
	chomp $l;
	my ($chr,$pos)=split(/\t/,$l);
	my $id=$chr."_".$pos;
	$var{$id}=$l;
}
close NORMAL;

while(my $l = <TUMOR>)	{
	chomp $l;
	my ($chr,$pos)=split(/\t/,$l);
	my $id=$chr."_".$pos;
	next if defined $var{$id};
	print "$l\n";
}
close TUMOR;
	
