## Perl script to add gene annotation for the right and left chromosomal positions for SV. The key used here is the line number, chromosome and position. Line number used due to redundancy found in using just chr_position pair.

use strict;
use warnings;

#check for presence of all i/p and o/p files
if (scalar(@ARGV) != 5)	
{
	die ( "Usage:\n GeneAnnotation.SV.pl <Input file containing let and right positions (provide full path)> <Input file with left end gene annotation (provide full path)> <Input file with right end gene annotation  (provide full path)> <Output file for left annotation (provide full path)> <Output file for right annotation (provide full path)> \n" );
}

my %hashinput_leftend=();
my %hashinput_rightend=();
my %hashleftend=();
my %hashrightend=();

#Form a unique key with line number, chr, pos
open (INPUT,"$ARGV[0]");
while(my $line = <INPUT>)	
{
	chomp $line;
	my @array = split('\t',$line);
	my $uniq = $array[0]."_".$array[1]."_".$array[2];
	$hashinput_leftend{$uniq} = join ("^",$array[0],$array[1],$array[2]);
}
close INPUT;

#Form a unique key with line number, chr, pos
open (INPUT,"$ARGV[0]");
while(my $line = <INPUT>)	
{
	chomp $line;
	my @array = split('\t',$line);
	my $uniq = $array[3]."_".$array[4]."_".$array[5];
	$hashinput_rightend{$uniq} = join ("^",$array[3],$array[4],$array[5],$array[6],$array[7]);
}
close INPUT;

#Form a unique key with line number, chr, pos
open (LEFTEND,"$ARGV[1]");
while(my $line = <LEFTEND>)	
{
	chomp $line;
	my @array = split('\t',$line);
	my $uniq = $array[8]."_".$array[5]."_".$array[6];
	$hashleftend{$uniq} = join ("^",$array[3]);
}
close LEFTEND;

#Form a unique key with line number, chr, pos
open (RIGHTEND,"$ARGV[2]");
while(my $line = <RIGHTEND>)	
{
	chomp $line;
	my @array = split('\t',$line);
	my $uniq = $array[8]."_".$array[5]."_".$array[7];
	$hashrightend{$uniq} = join ("^",$array[3]);
}
close RIGHTEND;

# If the keys match between the hashes, write to output files
open (OUTLEFT,'>',"$ARGV[3]");
foreach my $find (keys %hashinput_leftend)
{
	if (defined $hashleftend{$find})
	{
		print OUTLEFT "$hashinput_leftend{$find}^$hashleftend{$find}\n";
	}
	else
	{
		print OUTLEFT "$hashinput_leftend{$find}^NOGENE\n"
	}
}
close OUTLEFT;

open (OUTRIGHT,'>',"$ARGV[4]");
foreach my $find (keys %hashinput_rightend)
{
	if (defined $hashrightend{$find})
	{
		print OUTRIGHT "$hashinput_rightend{$find}^$hashrightend{$find}\n";
	}
	else
	{
		print OUTRIGHT "$hashinput_rightend{$find}^NOGENE\n"
	}
}
close OUTRIGHT;




