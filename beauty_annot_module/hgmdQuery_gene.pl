#! /usr/bin/perl
use DBI;

$filename=$ARGV[0];
$hgmd_db='db/hgmd_dup.sqlite';
my $dbh = DBI->connect("dbi:SQLite:dbname=".$hgmd_db, "", "", {RaiseError => 1, AutoCommit => 1});


open(IN, "<", $filename)||die$!;
$filename =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $filename."HGMD")||die$!;
if($filename =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	@header=split(/\t/, $h);
	print OUT join("\t", @header)."\tDiseaseByGene (HGMD)\tOMIM Ref (HGMD)\tComments (HGMD)\tMutations Per Gene (HGMD)\n";
}

$hgmdGene=$dbh->prepare("SELECT disease,omimid,comments,mut_total FROM allgenes WHERE gene = ?");
#$hgmdGene=$dbh->prepare("SELECT disease,omimid,comments,mut_total,go_terms_name FROM allgenes WHERE gene = ?");

$i=0;
while(<IN>){
	chomp;
	@line=split(/\t/, $_);
	
	if($line[5] eq ""){
		pop(@line);
		print OUT join("\t", @line)."\n";
	}
	else{
		
		#### ROUND 1 HGMD RSIDS ####
		$acc_num="";
		$hgmdGene->execute($line[5]);
		my $gene = $hgmdGene->fetchall_arrayref();
		@condit=();@omim=();@comments=();@totMut=(); #@goANNOT=();
			push(@condit, $$rowArr[0]);
			push(@omim, $$rowArr[1]);
			push(@comments, $$rowArr[2]);
			push(@totMut, $$rowArr[3]);
			### GENE GO ANNOTATION?? WANT IT, uncomment 3 places & add to print out
			#push(@goANNOT, $$rowArr[2]);
			
		
		@condit2=uniq(@condit);
		$disease=join("|", @condit2);
		@omim2=uniq(@omim);
		$omim1=join("|", @omim2);
		@comments2=uniq(@comments);
		$comments1=join("|", @comments2);
		@totMut2=uniq(@totMut);
		$totMut1=join("|", @totMut2);
		
		pop(@line);
		print OUT join("\t", @line)."\t";
		print OUT "$disease\t$omim1\t$comments1\t$totMut1\n";
			
		#if($i>40){last;}
		$i++;
	}
}


$dbh->disconnect();
close(OUT);

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
