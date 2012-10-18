#!/usr/local/biotools/perl/5.10.0/bin/perl        
use strict;
use warnings;

my $usage = "Usage:cat <File1.txt to compare> | report_original.pl <File2.txt with original format> \n";
my $file = shift or die $usage;
my %pos=();

#Form a unique key with chr_pos 
while(<STDIN>)	
{
    next if (/\#/);
    chomp;
    my @array = split('\t',$_);
    my $id = $array[0]."_".$array[1];
    $pos{$id} = 1;
}

#Form a unique key with chr_pos
open (FILE,"<$file") or die $!;
while(my $line = <FILE>)	
{
    next if ($line =~ m/\#/);
    chomp $line;
    my @array = split('\t',$line);
    my $id = $array[0]."_".$array[1];
    print $line."\n" if (exists $pos{$id});
}

close FILE;

