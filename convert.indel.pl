use strict;
use warnings;

open FH, "$ARGV[0]" or die "";
my $sign;
my $bases;
while(my $l = <FH>){
	chomp $l;
	my @a = split(/\t/,$l);
	#chr1	5559708	5559708	A	AAGAT	4	5	7
	#chr1	6599719	6599720	AT	A	1	5	6	
	if (length($a[3]) == 1)	{
		$sign="+";
		$bases=substr( $a[4],1,length($a[4]));
		print "$a[0]\t$a[1]\t$a[1]\t$sign$bases:$a[6]/$a[7]\n";
	}
	else	{
		$sign="-";
		$bases=substr( $a[3],1,length($a[3]));
		print "$a[0]\t$a[1]\t$a[2]\t$sign$bases:$a[6]/$a[7]\n";
	}
}
close FH;
	
	
		
		
