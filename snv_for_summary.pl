use strict;
use warnings;

my $gene=shift @ARGV;
my @files = @ARGV;
	

my %hash_g=();
open FH, "$gene" or die "cannot open $gene : $! \n";

while(my $l = <FH>)	{
	chomp $l;
	$l =~ s/^\s+//; #remove leading spaces
	$l =~ s/\s+$//; #remove trailing spaces
	$hash_g{$l}=$l;
}
close FH;






		

