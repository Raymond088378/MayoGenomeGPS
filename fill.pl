#!/usr/local/biotools/perl/5.10.0/bin/perl


# Fills all lines in a tab-delimited file with '-''s in place of nulls
use strict;
use warnings;

open FH, "$ARGV[0]" or die "";
my $header=<FH>;
chomp $header;

# This might have been supposed to be \s+ for multiple spaces in the header
my @head=split('/s+',$header);
print "$header\n";

while (my $l = <FH>)	{
	chomp $l;
	my @a = split(/\t/,$l);
	for (my $i = 0 ; $i <=$#head;$i++)	{
		if (length($a[$i]) == 0)	{
			$a[$i]='-';
		}
	}
	print join("\t",@a);
	print "\n";
}
close FH;		

