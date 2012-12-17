#! /usr/bin/perl
use DBI;

$delin="->";
$filename=$ARGV[0];
$hgmd_db='db/hgmd.sqlite';
my $dbh = DBI->connect("dbi:SQLite:dbname=".$hgmd_db, "", "", {RaiseError => 1, AutoCommit => 1});

open(IN, "<", $filename)||die$!;
$filename =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $filename."_HGMD")||die$!;
if($filename =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@header=split(/\t/, $h);
	print OUT join("\t", @header)."\tDiseaseByVariantHGMD)\tGeneFromVariant(HGMD)\tBases(HGMD)\tVariantTag(HGMD)\tPubmed(HGMD)\tAccessionNumber(HGMD)\n";
}

$hgmdRSID=$dbh->prepare("SELECT hgmd_acc FROM dbsnp WHERE dbsnp_id = ?");
$hgmdPOS=$dbh->prepare("SELECT acc_num FROM hg19_coords WHERE chromosome = '?' AND (coordSTART <= ? AND coordEND >= ?)");

%tagDesc=(
	'DFP'=>'Disease-Associated Polymorphisms with Supporting Functional Evidence',
	'DM'=>'Disease-Causing Mutation',
	'DM?'=>'Disease-Causing Mutation (Degree of Doubt)',
	'DP'=>'Disease-Associated Polymorphism',
	'FP'=>'Functional Polymorphisms',
	'FTV'=>'Predicted to Alter Gene Product (No Disease Association)'
);

while(<IN>){
		
	chomp;
	@line=split(/\t/, $_);
	
		
	#### ROUND 1 HGMD RSIDS ####
	$acc_num="";
	if($line[2] ne "."){
		$hgmdRSID->execute($line[2]);
		my $dbsnp = $hgmdRSID->fetchall_arrayref();
		$row1 = shift @$dbsnp;
		if(@$row1[0] ne ""){
			$acc_num=@$row1[0];
		}	
	}
	
	#### ROUND 2 HGMD POSITIONS ####
	$acc_num2="";
	if($line[2] eq "."){
		$hgmdPOS->execute($line[0],$line[1],$line[1]);
		my $hg19 = $hgmdPOS->fetchall_arrayref();
		$row2 = shift @$hg19;
		if(@$row2[0] ne ""){
			$acc_num2=@$row2[0];
		}	
	}
	
	#### UNION ACC OF RSID & POS HGMD ####
	@condit=();@genes=();@bases=();@aminos=();@tags=();@pmids=();
	if($acc_num ne "" || $acc_num2 ne ""){
		if($acc_num eq ""){$acc_num=$acc_num2;}
		##### TABLES DO NOT HAVE THE SAME schema.
		my $mut = $dbh->selectall_arrayref("SELECT disease,gene,base,amino,tag,pmid FROM mutation WHERE acc_num = \'$acc_num\';");
			foreach $rowArr (@$mut){
				$ref1="";$alt1="";
				$$rowArr[2] =~ s/[a-z]//g;
				($refB,$altB)=split("\-", $$rowArr[2]);
				@refArr=split("", $refB);
				@altArr=split("", $altB);
				for($y=0;$y<3;$y++){
					if($refArr[$y] ne $altArr[$y]){ 
						$ref1=$refArr[$y];
						$alt1=$altArr[$y];
					}
				}
				if($line[3] ne $ref1 && $line[4] ne $alt1){next;} ##{print "SNP DO NOT MATCH!\n\n"}
				push(@condit, $$rowArr[0]);
				push(@genes, $$rowArr[1]);
				push(@bases, $ref1."".$delin."".$alt1);
				push(@aminos, $$rowArr[3]);
				push(@tags, $$rowArr[4]);
				push(@pmids, $$rowArr[5]);				
			}
		my $mut = $dbh->selectall_arrayref("SELECT disease,gene,deletion,tag,pmid FROM deletion WHERE acc_num = \'$acc_num\';");
			foreach $rowArr (@$mut){
				$change="";
				if($$rowArr[2] =~ m/([a-z]+)/){
					$change="$1".$delin."";
				}
				push(@condit, $$rowArr[0]);
				push(@genes, $$rowArr[1]);
				push(@bases, uc($change));
				push(@tags, $$rowArr[3]);
				push(@pmids, $$rowArr[4]);				
			}
		my $mut = $dbh->selectall_arrayref("SELECT disease,gene,insertion,tag,pmid FROM insertion WHERE acc_num = \'$acc_num\';");
			foreach $rowArr (@$mut){
				$change="";
				if($$rowArr[2] =~ m/([a-z]+)/){
					$change="".$delin."$1";
				}
				push(@condit, $$rowArr[0]);
				push(@genes, $$rowArr[1]);
				push(@bases, uc($change));
				push(@tags, $$rowArr[3]);
				push(@pmids, $$rowArr[4]);				
			}
		my $mut = $dbh->selectall_arrayref("SELECT disease,gene,wildtype,insertion,tag,pmid FROM indel WHERE acc_num = \'$acc_num\';");
			foreach $rowArr (@$mut){
				$change="";
				if($$rowArr[2] =~ m/([a-z]+)/){
					$change="$1".$delin."".$$rowArr[3];
				}
				push(@condit, $$rowArr[0]);
				push(@genes, $$rowArr[1]);
				push(@bases, uc($change));
				push(@tags, $$rowArr[4]);
				push(@pmids, $$rowArr[5]);				
			}
	}
	
	@condit2=uniq(@condit);
	$disease=join("|", @condit2);
	@genes2=uniq(@genes);
	$gn=join("|", @genes2);
	@bases2=uniq(@bases);
	$bs=join("|", @bases2);
	@aminos2=uniq(@aminos);
	$aa=join("|", @aminos2);
	@tags2=uniq(@tags);
	@cvTags=map{ $tagDesc{$_} } @tags2;
	$tg=join("|", @cvTags);
	@pmids2=uniq(@pmids);
	$pm=join("|", @pmids2);
	## mutation add Amino Acid change
	if($aa ne ""){ $bs = $bs."; ".$aa }
	
	
	print OUT join("\t", @line)."\t";
	print OUT "$disease\t$gn\t$bs\t$tg\t$pm\t$acc_num\n";
}


$dbh->disconnect();
close(OUT);

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
