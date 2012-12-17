#! usr/bin/perl

$input = $ARGV[0];
open(REF, "<", $ARGV[1])||die$!;
open(IN, "<", $input)||die$!;
open(OUT, ">", $input."_DFlAG")|| die $!;

%druggableSources;
$hg=<REF>;
$hg=~s/\n//g;
@header=split(/\t/,$hg);
while(<REF>){
	chomp;
	@rr=split(/\t/, $_);
	@arrAble=();
	for($n=1;$n<=$#rr;$n++){
		if($rr[$n] eq "YES"){
			push(@arrAble, $header[$n]);
		}
	}
	$druggableSources{$rr[0]}=join("|",@arrAble);
#	print "$rr[0]\t".join("|",@arrAble)."\n";
}

if($input =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@inHd=split(/\t/, $h);
	pop(@inHd);
	print OUT join("\t", @inHd)."\tPossibleDruggable(DruggableGenome)\n";
}



$i=0;
while(<IN>){
	chomp;
	@ln=split(/\t/,$_);
	$gene=$ln[5];
	if($gene ne ""){
		pop(@ln)
	}

	print OUT join("\t", @ln)."\t".$druggableSources{$gene}."\n";

#	if($i > 35){last;}
#	$i++;
}
close(IN);