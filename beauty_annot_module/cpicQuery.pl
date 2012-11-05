#! /usr/bin/perl

use DBI;



$delin="->";

$filename=$ARGV[0];

$db='/data2/bsi/tertiary/Kocher_Jean-Pierre_m026645/s112224.Beauty_Annotation/Im_Clinic_Mod/db/cpicTemp.sqlite';

my $dbh = DBI->connect("dbi:SQLite:dbname=".$db, "", "", {RaiseError => 1, AutoCommit => 1});



open(IN, "<", $filename)||die$!;

$filename =~ s/\.[A-Za-z]{3}$//;

open(OUT, ">", $filename."_CPIC")||die$!;

if($filename =~ /\_1(\_GENE)*$/){

	$h=<IN>;

	$h=~s/\n//g;

	@header=split(/\t/, $h);

	print OUT join("\t", @header)."\tVIP_Gene(CPIC)\tVIP_Variant(CPIC)\n";

}



$baseGene=$dbh->prepare(" SELECT drug,status FROM baseCpic WHERE gene = ?");

$vipVar=$dbh->prepare("SELECT decriptions FROM fullGeneCpic WHERE type = 'rsid' AND class = 'VIP' AND typeId = ?");



##### LEGEND

%legend=(

	'DG'=>'Dosing Guideline information is available',

	'DL'=>'Drug Label information is available',

	'CA'=>'High-level Clinical Annotation is available',

	'VA'=>'Variant Annotation is available',

	'VIP'=>'VIP information is available',

	'PW'=>'Pathway is available'

);



$i=0;

while(<IN>){

	chomp;

	@line=split(/\t/, $_);

	

	$vipGENE="";

	if($line[5] ne ""){

		@tmpVip=();

		$baseGene->execute($line[5]);

		my $baseG = $baseGene->fetchall_arrayref();

		foreach $rowArr (@$baseG){

			push(@tmpVip, $$rowArr[0]."[".$$rowArr[1]."]");

		}

		$vipGENE=join("|", @tmpVip);

	}

	

	$vipVAR="";

	if($line[2] =~ /^rs\d+/){

		@tmpVipV=();

		$vipVar->execute($line[2]);

		my $baseV = $vipVar->fetchall_arrayref();

		foreach $rowArr (@$baseV){

			push(@tmpVipV, $$rowArr[0]);

		}

		$vipVAR=join("|", uniq(@tmpVipV));

	}

	

	

	print OUT join("\t", @line)."\t";

	print OUT "$vipGENE\t$vipVAR\n";

	

#	if($i>8){last;}

#	$i++;

}

close(IN);





sub uniq {

    return keys %{{ map { $_ => 1 } @_ }};

}