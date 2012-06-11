use strict;
use warnings;

my $vcf=shift @ARGV;
open FH, "$vcf" or die "can not open $vcf : $!\n";

while(my $l = <FH>)	{
	chomp $l;
	if ($l =~ m/^##/)	{
		print "$l\n";
	}
	elsif ($l =~ m/^#/)	{
		print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		print "$l\n";
	}
	else	{
		print "$l\n";
	}
}
close FH;	
