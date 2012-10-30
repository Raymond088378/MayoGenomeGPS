#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;
use Getopt::Std;

our ($opt_i, $opt_o,$opt_f);
#print "RAW paramters: @ARGV\n";
getopt('iof');
if ( (!defined $opt_i) && (!defined $opt_o)  && (!defined $opt_f)) {
        die ("Usage: $0 \n\t-i [nput file] \n\t-o [utput file]\n\t-f [ormat (Start/Position)]\n");
}
else    {
	my $source=$opt_i;
	my $output=$opt_o;
	my $format= $opt_f;
	open FH, "$source" or die "can not open $source : $! \n";
	open OUT, ">$output" or die " can not open $output : $! \n";
	my %report=();
	my %chrvalue = ("chrX"=>23,"chrY"=>24,"chrM"=>25);
	for (my $i=1; $i<23; $i++) {
		$chrvalue{"chr".$i} = $i;
	}
	my $chr=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if (\$i == \"Chr\") {print i} } }' $source`;chomp $chr;$chr=$chr-1;
	
	#print "$chr\n";
	my $pos;
	if($format eq 'Start')    {
        $pos=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if (\$i == \"Start\") {print i} } }' $source`;chomp $pos;$pos=$pos-1;
	}	
	else    {
        $pos=`awk -F '\t' '{ for(i=1;i<=NF;i++){ if (\$i == \"Position\") {print i} } }' $source`;chomp $pos;$pos=$pos-1;
	}
	while(my $l = <FH>)	{
		chomp $l;
		if($. < 3)	{
			print OUT "$l\n";
		}	
		else	{
			my @a = split(/\t/,$l);
			push(@{$report{$a[$chr]}{$a[$pos]}},$l);
		}
	}
	close FH;
	#print "printing the sorted report\n";
	foreach my $c (sort {$chrvalue{$a}<=>$chrvalue{$b}} keys %report)	{
		foreach my $p (sort {$a<=> $b} keys %{$report{$c}})	{
			my $num_records=$#{$report{$c}{$p}};
			for(my $i=0; $i <=$num_records; $i++)	{
				print OUT "$report{$c}{$p}[$i]";
				print OUT "\n";
			}
		}
	}
}
