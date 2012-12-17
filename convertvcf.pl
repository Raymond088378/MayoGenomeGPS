#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;

my $vcf = shift @ARGV;

open FH, "$vcf"  or die " failed to open $vcf : $!\n";
my @ar;
while (my $l = <FH>)	{
	if ($l =~ m/^##/)	{
		print "$l";
	}
	elsif ($l =~ m/^#/)	{
		chomp $l;
		@ar = split(/\t/,$l);
		print "$l\n";
	}
	else	{
		chomp $l;
		my @a = split(/\t/,$l);
		### 7th field is pass field
		if($a[6] ne "PASS")	{
			my @f=(0 .. 5);
			my @la=(7 .. $#ar);
			print join ("\t",@a[@f],"PASS",@a[@la]) . "\n";
		}
		else	{
			print "$l\n";
		}
	}
}
close FH; 
