#! /usr/bin/perl

### So far...this is dumb, it relys solely on RSID!
$input=$ARGV[0];
open(IN, "<", $input)|| die $!;
open(RSID, "<", $ARGV[1])|| die $!;
open(OUT, ">", $input."_RSID")|| die $!;

if($input =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@header=split(/\t/, $h);
	print OUT join("\t", @header)."\tPharmicogenetic_Flag (dbSNP)\n";
}

%pharmico;
while(<RSID>){
	next if($_ !~ /^rs\d+/);
	chomp;
	($rs,$gene)=split(/\t/, $_);
	$pharmico{$rs}=$gene;
}

#print "Number of Library Vars: ".(keys %pharmico)."\n";

#$fnd=0;
while(<IN>){
	chomp;
	@var=split(/\t/, $_);
	if($var[2] !~ /^rs\d+/){
		print OUT "$var[0]\t$var[1]\t$var[2]\t$var[3]\t$var[4]\tNO\n";
	}
	else{
		if( exists( $pharmico{$var[2]}) ){
			print OUT "$var[0]\t$var[1]\t$var[2]\t$var[3]\t$var[4]\tYES\n";
#			$fnd++;
		}
		else{
			print OUT "$var[0]\t$var[1]\t$var[2]\t$var[3]\t$var[4]\tNO\n";
		}
	}
}

#print "Found ".$fnd." Matches\n";
close(IN);close(RSID);close(OUT);