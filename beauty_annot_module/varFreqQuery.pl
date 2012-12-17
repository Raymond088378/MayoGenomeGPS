#! /usr/bin/perl
use DBI;
use Getopt::Long;
use List::Util qw[max];
my %options;
GetOptions(\%options,"i:s","db:s","possible|?");

#########################
## REQUIRE DEFINITIONS ##
#########################
$myName="VarFreq";
@myColumns=('DiseaseByVariant','GeneFromVariant','Bases','VariantTag','Pubmed','AccessionNumber');

if ($options{'possible'}) { print "$myName|".join(',',@myColumns)."\n"; exit; }

$filename=$options{'i'};
$db='db/freqOtherMAF.sqlite';
if ($options{'db'}) { $db=$options{'db'}; }
my $dbh = DBI->connect("dbi:SQLite:dbname=".$db, "", "", {RaiseError => 1, AutoCommit => 1});

open(IN, "<", $filename)||die$!;
$filename =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $filename."_VARFREQ")||die$!;
if($filename =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@header=split(/\t/, $h);
	print OUT join("\t", @header)."\tMaxKGenome(kGenome)\tAnyGreaterThan10%(kGenome)\n";
}

$hRSID=$dbh->prepare("SELECT amrAf,asnAf,afrAf,eurAf FROM kgenomeFreq WHERE rsid = ?");
#$hPOS=$dbh->prepare("SELECT acc_num FROM hg19_coords WHERE chromosome = '?' AND (coordSTART <= ? AND coordEND >= ?)");

$i=0;
while(<IN>){
	chomp;
	@line=split(/\t/, $_);
	
	if( $line[2] =~ /^rs\d+/i){
		$hRSID->execute($line[2]) or warn "Couldn't execute statement: " . $sth->errstr;
		if( defined($hRSID) ){
			my $freqRet = $hRSID->fetchall_arrayref();
			$freq = shift @$freqRet;
			print ">>".join("?",@$freqRet)."\n";
			print join("++",@$freq)."\n";
			print "MAX=".max(@$freq)."\n";
		}
	}
	
	print "$line[2]\n";
	
#	if($i > 5){last;}
	$i++;
}


$dbh->disconnect();
close(OUT);

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
