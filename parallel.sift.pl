#!/usr/local/biotools/perl/5.14.2/bin/perl

use threads;
use threads::shared;
use POSIX;

my $num_threads =shift @ARGV;
my $file = shift @ARGV;
my $sift_ref = shift @ARGV;
my $sift = shift @ARGV;
my $out = shift @ARGV;
my $output = shift @ARGV;
my @threads;
my $len=`cat $file | awk '\$0 !~ /^#/' | wc -l`;

for ( my $count = 1; $count <= $num_threads; $count++) {
	my $start=ceil(($count-1)*($len/$num_threads)+1);
	my $end=ceil(($len/$num_threads)*$count);
	if ($end > $len)	{
		$end = $len;
	}
	my $t = threads->create(\&soft, $count, $sift_ref, $file, $sift, $output, $start, $end );
	push(@threads,$t);
}
foreach (@threads) {
	my $num = $_->join;
}
open OUT , ">>$out" or die "";

`cp $output/1.predictions.tsv $out`;
`rm $output/1.predictions.tsv`;
for ( my $count = 2; $count <= $num_threads; $count++) {
	my $out1="$output/$count.predictions.tsv";
	open OUT1, "$out1" or die "";
	while(<OUT1>)	{
		next if ($_ =~ /^#/);
		print OUT $_;
	}
	close OUT1;
	`rm $out1`;
}

sub soft {
	my $num = shift;
	my $sift_ref = shift;
	my $input = shift;
	my $sift = shift;
	my $output = shift;
	my $log = $input . ".$num" . ".run";
	my $start = shift;
	my $end= shift;
	my $total = $end -$start+1;
	`cat $input | head -n $end | tail -n $total > $input.$num`;
	`cd $sift`;
	`perl $sift/SIFT_exome_nssnvs.pl -i $input.$num -d $sift_ref -o $output/ -A 1 -B 1 -J 1 -K 1 > $log`;
	open FH, "$log" or die "";
	my $head=<FH>;my $id;
	$id=$1 if ($head =~ /Your job id is (\d+)/);
	`mv $output/$id/${id}_predictions.tsv $output/$num.predictions.tsv`;
	`rm $log`;
	`rm $input.$num`;
	`rm -R $output/$id`;
	`cd $output`;
}

