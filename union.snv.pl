use strict;
use warnings;

die "Usage:$0 <file w/ list of files> <output file>\n" if ($#ARGV != 1); 
my $list=$ARGV[0];
my $out=$ARGV[1];
open FH, "$list" or die "can not open $list :$!\n";
open OUT, ">$out" or die "can not open $out :$!\n";

my $IGV=0;
my $CHR=1;
my $POS=2;
my $ALT=19;
my $dbSNP=3;
my $REF=18;
my $START_INFO=20;
my $STOP_INFO=25;
my $num_col_samples=$STOP_INFO -$START_INFO +1;
my $SIFT_START=26;
my $SIFT_STOP=54;
my $SNPEFF_START=55;
my %chrvalue = ("chrX"=>23,"chrY"=>24,"chrM"=>25);
for (my $i=1; $i<23; $i++) {
	$chrvalue{"chr".$i} = $i;
}
my @samples;  ## store the sample names and print in the same order
my $i =0;
my @annot_ref=($dbSNP .. $REF);
my @annot_sift=($SIFT_START .. $SIFT_STOP);
my @sample=($START_INFO .. $STOP_INFO);
my (%igv,%ref, %sift, %snpeff, %sample_info);
my ($head1,$head2,$head3,$head4);
my $count=1;
while(my $l = <FH>){
    chomp $l;
    open FILE, "$l" or die "can not open $l : $!\n";
	my $prev=0;
	print "Reading $l\n";
	while(my $k = <FILE>){
		chomp $k;
		my @a = split("\t",$k);
		my $last_col=$#a;
		my @annot_snpeff=($SNPEFF_START .. $last_col);
		if ( $. == 1){
			$samples[$i]=$a[$START_INFO];
			$i++;
		}
		elsif ( $. == 2){
			$head1=join("\t",@a[@annot_ref]);
			$head2=join("\t",@a[@annot_sift]);
			$head3=join("\t",@a[@annot_snpeff]);
			$head4=join("\t",@a[@sample]);
		}
		else	{
			my $id="$a[$CHR],$a[$POS],$a[$ALT]";
			##snpeff values	
			my $snpeff_value=join("\t",@a[@annot_snpeff]);
			push(@{$snpeff{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$snpeff_value);
			if ($id ne $prev)	{
				my $value=$#{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}};
				$value=$value+2;
				my $sample_value;
				my $sam;
				if ($value == $i){
					$sample_value=join("\t",@a[@sample]);
					push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sample_value);
				}
				else {
					for(my $ll=$value; $ll <$i;$ll++){
						$sam="n/a\tn/a\tn/a\tn/a\tn/a\tn/a";
						push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sam);	
					}
					$sample_value=join("\t",@a[@sample]);
					push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sample_value);	
				}	
				##refernce value
				my $ref_value=join("\t",@a[@annot_ref]);
				$ref{$a[$CHR]}{$a[$POS]}{$a[$ALT]}=$ref_value;
				## sift value
				my $sift_value=join("\t",@a[@annot_sift]);
				$sift{$a[$CHR]}{$a[$POS]}{$a[$ALT]}=$sift_value;
				$igv{$a[$CHR]}{$a[$POS]}{$a[$ALT]}=$a[$IGV];
			}
			$prev=$id;
		}
	}
	$count++;
	close FILE;
}
close FH;
print "Merging all the sample files\n";
print OUT "-". "\t" x 6 . "Allele Freuency";
my $a1=$ALT-$dbSNP-2;
print OUT "\t" x $a1;
print OUT join ("\t\t\t\t\t\t",@samples);
print OUT "\t\t\t\t\t\tSIFT Annotation";
my $a2 = $SIFT_STOP - $SIFT_START +1;
print OUT "\t" x $a2;
print OUT "SNPEFF Annotation";

print OUT "\n";
print OUT "IGV Link\tChr\tPosition\t$head1\tAlt\t";
for(my $i=0; $i <=$#samples; $i++)	{
	print OUT "$head4\t";
}
print OUT "$head2\t$head3\n";	
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
foreach my $c (sort {$chrvalue{$a}<=>$chrvalue{$b}} keys %ref)	{
	foreach my $p (sort {$a<=> $b} keys %{$ref{$c}})	{
		foreach my $a (sort keys %{$ref{$c}{$p}})	{	
			my $key = "$c,$p,$a";
			my @required_snpeff=uniq(@{$snpeff{$c}{$p}{$a}});
			my $num_snpeff=scalar(@required_snpeff);
			for(my $rows=0; $rows < $num_snpeff; $rows++)	{
				my @a=split(/,/,$key);
				my $num_samples=$#samples;
				print OUT "$igv{$c}{$p}{$a}\t$a[0]\t$a[1]\t";
				print OUT "$ref{$c}{$p}{$a}\t$a[2]\t";
				my $values=$#{$sample_info{$c}{$p}{$a}};
				print OUT join ("\t", @{$sample_info{$c}{$p}{$a}});
				if ($values != $num_samples){
					for(my $k=$values+1; $k<$num_samples;$k++){
						print OUT "\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a";
					}
					print OUT "\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a";
				}
				print OUT "\t$sift{$c}{$p}{$a}\t$snpeff{$c}{$p}{$a}[$rows]\t";	
				print OUT "\n";
			}
		}
	}
}
close OUT;
print "$out is generated\n";	
    
	    
	    
	    
