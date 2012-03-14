use strict;
use warnings;

open FH, "$ARGV[0]" or die "";
my $count=1;
while(my $l = <FH>){
	chomp $l;
	if ($l =~ /^#/)	{
	
	}
	else	{
        my @a = split(/\t/,$l);
        my $start=$a[1]-1;
        print "link$count\t$a[0]\t$start\t$a[1]\n";
            $start=$a[4]-1;
            print "link$count\t$a[3]\t$start\t$a[4]\n";
        $count++;
	}
}
close FH;
	
