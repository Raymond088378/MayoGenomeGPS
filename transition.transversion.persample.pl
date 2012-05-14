# script to calculate transition to trnsversion ratio

#Transition [ A-G G-A C-T T-C ]
#Transversion [ C-G G-C G-T T-G A-T T-A C-A A-C ]

use strict;
use warnings;
my $file = shift @ARGV;
open (DAT,"$file");

my $ref= shift @ARGV;$ref=$ref-1;
my $alt = shift @ARGV; $alt=$alt-1;
my $transition=0;
my $transversion=0;
while ( my $l = <DAT>)		{
	chomp $l;
	my @a = split (/\t/,$l);
	$a[$ref] =~ tr/[a-z]/[A-Z]/;
	$a[$alt] =~ tr/[a-z]/[A-Z]/;
	if( (($a[$ref] eq 'A') && ($a[$alt] eq 'G')) || 
		(($a[$ref] eq 'G') && ($a[$alt] eq 'A')) ||
		(($a[$ref] eq 'C') && ($a[$alt] eq 'T')) ||
		(($a[$ref] eq 'T') && ($a[$alt] eq 'C')) )	{
			$transition++;
	}
	else	{
		$transversion++;
	}
}
my $ratio=$transition/$transversion;$ratio=sprintf ("%.2f",$ratio);
print "$ratio\n";
		
