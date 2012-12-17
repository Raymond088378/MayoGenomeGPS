#!/usr/local/biotools/perl/5.10.0/bin/perl                                                                                                                                           
use strict;

my $id = shift or die "";
my $threshold = shift or die "";

while (<>) {
    chomp;
   if ($. > 1)	{
    my @fields = split/\t/;
    next if ( @fields < 0);
    next if ( ($threshold>1) && ($fields[3]<$threshold));
    next if ( ($threshold<1) && ($fields[3]>$threshold));
    
    print "chr".$fields[0]."\t".$fields[1]."\t".$fields[2]."\t".$id."\t".$fields[3]."\t+\n";
	}
}
    
