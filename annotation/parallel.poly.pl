#!/usr/local/biotools/perl/5.14.2/bin/perl

use threads;
use threads::shared;
use POSIX;

my $num_threads =shift @ARGV;
my $file = shift @ARGV;
my $genome = shift @ARGV;
my $poly = shift @ARGV;
my $out = shift @ARGV;
my @threads;
my $len=`cat $file | awk '\$0 !~ /^#/' | wc -l`;

for ( my $count = 1; $count <= $num_threads; $count++) {
	my $start=ceil(($count-1)*($len/$num_threads)+1);
	my $end=ceil(($len/$num_threads)*$count);
	if ($end > $len)	{
		$end = $len;
	}
	my $t = threads->create(\&pph, $count, $genome, $file, $poly, $start, $end );
	push(@threads,$t);
}
foreach (@threads) {
	my $num = $_->join;
}
open OUT , ">>$out" or die "";
`cp $file.1.out $out`;
`rm $file.1.out`;
for ( my $count = 2; $count <= $num_threads; $count++) {
	my $out1="$file.$count.out";
	open OUT1, "$out1" or die "";
	while(<OUT1>)	{
		next if ($_ =~ /^#/);
		print OUT $_;
	}
	close OUT1;
	`rm $out1`;
}

sub pph {
	my $num = shift;
	my $genome = shift;
	my $input = shift;
	my $poly = shift;
	my $output = $input . ".$num" . ".out";
	my $log = $input . ".$num" . ".log";
	my $start = shift;
	my $end= shift;
	my $total = $end -$start+1;
	`cat $input | head -n $end | tail -n $total > $input.$num`; 
	`/usr/local/biotools/perl/5.10.0/bin/perl $poly/bin/mapsnps.pl -v 0 -A -g $genome $input.$num > $output 2> $log`;
	`rm $log`;
	`rm $input.$num`;
}

