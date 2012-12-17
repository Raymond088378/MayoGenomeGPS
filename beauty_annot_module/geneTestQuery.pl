#! /usr/bin/perl
use DBI;

$filename=$ARGV[0];
$db='db/geneTests.sqlite';
my $dbh = DBI->connect("dbi:SQLite:dbname=".$db, "", "", {RaiseError => 1, AutoCommit => 1});

open(IN, "<", $filename)||die$!;
$filename =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $filename."TEST")||die$!;
if($filename =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@header=split(/\t/, $h); pop(@header);
	print OUT join("\t", (@header,"Disease (GeneTest)", "MIM Numbers (GeneTest)", 'Availability (GeneTest)'))."\n";
} 

$geneT=$dbh->prepare("SELECT disease_name,mim_nums,testAvailability FROM genePlus AS g JOIN geneTestCustomPlus AS t ON g.testRefId=t.testId WHERE g.gene = ?");

$i=0;
while(<IN>){
	chomp;
	@line=split(/\t/, $_);
	next if($line[5] eq "");

	$geneT->execute($line[5]);
	my $geneTestRet = $geneT->fetchall_arrayref();
	@disease=();@mims=();@avalibility=();
	foreach $rowArr (@$geneTestRet){
		push(@disease, $$rowArr[0]);
		push(@mims, split(/\|/, $$rowArr[1]) );
		if($$rowArr[2] !~ /(na)/){
			push(@avalibility, $$rowArr[2]);
		}
	}
	
	@disease2=uniq(@disease);
	$dis=join("|", @disease2);
	@mims2=uniq(@mims);
	$mim=join("|", @mims2);
	@avalibility2=uniq(@avalibility);
	$avail=join("|", @avalibility2);
	
	pop(@line);
	print OUT join("\t", @line)."\t";
	print OUT "$dis\t$mim\t$avail\n";
}


$dbh->disconnect();
close(OUT);

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}