use strict;
use warnings;

my $vcf=shift @ARGV;
my $type=shift @ARGV;
open FH, "$vcf" or die "can not open $vcf : $!\n";

while(my $l = <FH>)	{
	chomp $l;
	if ($l =~ m/^##/)	{
		print "$l\n";
	}
	elsif ($l =~ m/^#/)	{
		print "##INFO=<ID=CAPTURE,Number=1,Type=Boolean,Description=\"variant in capture kit ot not\">\n";
		if ($type eq 'SNV')	{
			print "##INFO=<ID=CLOSE2INDEL,Number=1,Type=Boolean,Description=\"if a snp is close to indel for this sample\">\n";
		}
		print "$l\n";
	}
	else	{
		print "$l\n";
	}
}
close FH;	