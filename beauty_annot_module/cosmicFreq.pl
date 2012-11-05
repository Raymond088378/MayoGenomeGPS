#! /usr/bin/perl
use DBI;

$input=$ARGV[0];
$db='/data2/bsi/tertiary/Kocher_Jean-Pierre_m026645/s112224.Beauty_Annotation/Im_Clinic_Mod/db/freqCosmic.sqlite';
open(IN, "<", $input)|| die $!;
$input =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $input."_COSMIC")||die $!;

if($input =~ /\_1(\_GENE)*$/){
	$h=<IN>;
	$h=~s/\n//g;
	print OUT $h."\tMutationCDS (Cosmic)\tMutationDescription (Cosmic)\tZygosity (Cosmic)\tConfirmation (Cosmic)\tMutationsPerGene (Cosmic)\tMutationFrequency (Cosmic)\n";
}

my $dbh = DBI->connect("dbi:SQLite:dbname=".$db, "", "", {RaiseError => 1, AutoCommit => 1});

$mutation=$dbh->prepare("SELECT bicSampleRefId,mutCDS,mutDesc,mutZygosity,somaticStatus FROM cosmicMutations WHERE chr = ? AND start = ?");
$getSample=$dbh->prepare("SELECT gene FROM cosmicSampleCnts WHERE bicSampleId = ?"); ### think about adding tissue & histology
$allSamples=$dbh->prepare("SELECT bicSampleId,sampleCnt FROM cosmicSampleCnts WHERE gene = ? ");
#$allMutations=$dbh->prepare("SELECT COUNT(*) FROM cosmicMutations WHERE bicSampleRefId IN ( ? )");

$i=0;
while(<IN>){
	chomp;
	@var=split(/\t/, $_);
	$chr=$var[0];
	$chr =~ s/^chr//i;
	$chr =~ s/x/23/i;
	$chr =~ s/y/24/i;
	$chr =~ s/m/25/i;

	@cds=();@mut=();@zyg=();@source=();@genes=();
	$mutation->execute($chr,$var[1]);
	my $return = $mutation->fetchall_arrayref();
	foreach $ret1 (@$return){
		if($$ret1[0] ne ""){
			($refC,$altC) = $$ret1[1] =~ /([A-Z])>([A-Z])$/;
			next if( $var[3] ne $refC || $var[4] ne $altC );

			push(@cds, $$ret1[1]);
			push(@mut, $$ret1[2]);
			push(@zyg, uc($$ret1[3]));
			push(@source, $$ret1[4]);

			$getSample->execute($$ret1[0]);
			my $gene = $getSample->fetchall_arrayref()->[0][0];
			push(@genes, $gene);
			$i++;
		}
	}

	@bicSampleIds=();$sampleCnt=0;@finGene=();@finMAF=();
	@useableGenes=uniq(@genes);
	foreach $gn (@useableGenes){
		$allSamples->execute($gn);
		my $retSampl = $allSamples->fetchall_arrayref();
		foreach $ret2 (@$retSampl){
			push(@bicSampleIds, $$ret2[0]);
			$sampleCnt += $$ret2[1];
		}
		
		$allids = join(",",@bicSampleIds);
		$mutCnt = $dbh->selectrow_array("SELECT COUNT(*) FROM cosmicMutations WHERE bicSampleRefId IN ( $allids )");
		#$allMutations->execute( $allids );
		#my $mutCnt = $allMutations->fetchall_arrayref()->[0][0];
		#print $mutCnt."\n";
		if($mutCnt == 0){ print "Error: Cannot divide by Zero\n$allids\n\n"; exit;}
		push(@finGene, $gn." ".$mutCnt."/".$sampleCnt);
		$calc=($mutCnt/$sampleCnt);
		push(@finMAF,$calc);
	}


	@cds2=uniq(@cds);
	$cds1=join("|", @cds2);
	@mut2=uniq(@mut);
	$mut1=join("|", @mut2);
	@zyg2=uniq(@zyg);
	$zyg1=join("|", @zyg2);
	@source2=uniq(@source);
	$source1=join("|", @source2);
	$genes1=join(";", @finGene);
	$mafs1=join(";", @finMAF);	

	print OUT join("\t", @var)."\t";
	print OUT "$cds1\t$mut1\t$zyg1\t$source1\t$genes1\t$mafs1\n";

	
	#if($i>20){last;}
}

$dbh->disconnect();

close(OUT);



sub uniq {

    return keys %{{ map { $_ => 1 } @_ }};

}

