#script to comibine all sift results to get combine over all samples
use strict;	
#Usage: perl merge.sift.results.pl <siftid file> <path to sift reports> > merged.sift.result

my $list=shift @ARGV;
open FH, "$list" or die "";
my $head;
my $sift_size;
my %sift=();
while(my $l = <FH>)	{
	chomp $l;
	open FILE, "$l" or die "";
	while(my $m = <FILE>)	{
		chomp $m;
		if($. == 1)	{
			$head=$m;
			my @b=split(/\t/,$m);
			$sift_size = scalar(@b);
		}
		else	{
			my @b=split(/\t/,$m);
			my @c;
			for(my $i = 0; $i < $sift_size; $i++)	{
				my $len = length($b[$i]);
				if($len == 0)	{
					$b[$i] = '-';
				}
			}	
			$sift{$b[0]} = join("\t",@b);
		}
	}	
	close FILE;	
}
close FH;
print "$head\n";
foreach my $key (keys %sift)	{
	print "$sift{$key}\n";
}
undef %sift;

	
