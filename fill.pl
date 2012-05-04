use strict;
use warnings;

open FH, "$ARGV[0]" or die "";
my $header=<FH>;
chomp $header;
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

