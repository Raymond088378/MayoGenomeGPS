#!/usr/local/biotools/perl/5.10.0/bin/perl

use strict;
use warnings;

open FH, "$ARGV[0]" or die "";
my $sign;
my $bases;
while(my $l = <FH>){
	chomp $l;
	my @a = split(/\t/,$l);
	#chr1    761957  761957  1       0       A       AT      1       117     118
	#chr1    900717  900721  1       0       CTTAT   C       4       100     128
	
	if (length($a[5]) == 1)	{
		$sign="+";
		$bases=substr( $a[6],1,length($a[6]));
		print "$a[0]\t$a[1]\t$a[1]\t$sign$bases:$a[8]/$a[9]\n";
	}
	else	{
		$sign="-";
		$bases=substr( $a[5],1,length($a[5]));
		print "$a[0]\t$a[1]\t$a[2]\t$sign$bases:$a[8]/$a[9]\n";
	}
}
close FH;
	
	
		
		
