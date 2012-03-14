#!/usr/bin/perl

use strict;
use warnings;

my $usage = "cnv2bed.pl label";

my $label = shift or die $usage;

my $i=1;
while (<STDIN>){
    chomp;
    my @fields = split /\t/;
    my @coord = split(/:|\-/,$fields[1]);

    print $coord[0]."\t".$coord[1]."\t".$coord[2]."\t".
	$label."\t".$fields[3]."\t+\n";
    $i++;
}

