use strict;
use warnings;
die "Usage:$0 <file w/ list of files> <output file>\n" if ($#ARGV != 1);
my $list=$ARGV[0];
open FH, "$list" or die "can not open $list : $!\n";
my $out=$ARGV[1];
open OUT, ">$out" or die " can not open $out: $!\n";
my $IGV=0;
my $CHR=1;
my $START=2;
my $STOP=3;
my $ALT=10;
my $BASE=11;
my $dbSNP=4;
my $REF=9;
my $START_INFO=12;
my $STOP_INFO=13;
my $SNPEFF_START=14;
my @samples;  ## store the sample names and print in the same order
my $i =0;
my @annot_ref=($dbSNP .. $REF);
my @sample=($START_INFO .. $STOP_INFO);
my (%igv,%ref, %snpeff, %sample_info);
my ($head1,$head3,$head4);
my $count=1;
my %chrvalue = ("chrX"=>23,"chrY"=>24,"chrM"=>25);
for (my $i=1; $i<23; $i++) {
	$chrvalue{"chr".$i} = $i;
}
while(my $l = <FH>){
   	chomp $l;
	open FILE, "$l" or die "";
	my $prev=0;
   	print "Reading $l\n";
	while(my $k = <FILE>){
		chomp $k;
		my @a = split("\t",$k);
		my $last_col=$#a;
		my @annot_sseq=($SNPEFF_START .. $last_col);
		if ( $. == 1){
			$samples[$i]=$a[$START_INFO];
			$i++;
		}
		elsif ( $. == 2){
			$head1=join("\t",@a[@annot_ref]);
			$head3=join("\t",@a[@annot_sseq]);
			$head4=join("\t",@a[@sample]);
		}
		else	{
			my $id="$a[$CHR]\*$a[$START]\*$a[$STOP]\*$a[$ALT]\*$a[$BASE]";
			my $id_a="$a[$STOP]\*$a[$ALT]\*$a[$BASE]";
			##sseq values	
			my $snpeff_value=join("\t",@a[@annot_sseq]);
			push(@{$snpeff{$a[$CHR]}{$a[$START]}{$id_a}},$snpeff_value);
			if ($id ne $prev)	{
				my $value=$#{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}};
				$value=$value+2;
				my $sample_value;
				my $sam;
				if ($value == $i){
					$sample_value=join("\t",@a[@sample]);
					push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sample_value);
				}
				else {
					for(my $ll=$value; $ll <$i;$ll++){
						$sam="n/a\tn/a";
						push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sam);	
					}
					$sample_value=join("\t",@a[@sample]);
					push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sample_value);	
				}		
				##refernce value
				my $ref_value=join("\t",@a[@annot_ref]);
				$ref{$a[$CHR]}{$a[$START]}{$id_a}=$ref_value;
				$igv{$a[$CHR]}{$a[$START]}{$id_a}=$a[$IGV];
			}
			$prev=$id;
		}
	}
	$count++;
	close FILE;
}
close FH;
print OUT "-\t";
print OUT "\t" x $BASE;
print OUT join ("\t\t",@samples);
print OUT "\t\tSNPEFF Annotation";
print OUT "\n";
print OUT "IGV Link\tChr\tStart\tStop\t$head1\tAlt\tBase-Length\t";
for(my $i=0; $i <=$#samples; $i++)	{
	print OUT "$head4\t";
}
print OUT "$head3\n";	
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
print "Merging all teh sample files\n";
foreach my $c (sort {$chrvalue{$a}<=>$chrvalue{$b}} keys %ref)	{
	foreach my $p (sort {$a<=> $b} keys %{$ref{$c}})	{
		foreach my $a (keys %{$ref{$c}{$p}})	{	
			my @required_snpeff=uniq(@{$snpeff{$c}{$p}{$a}});
			my $num_snpeff=scalar(@required_snpeff);
			for(my $rows=0; $rows < $num_snpeff; $rows++)	{
				my $key = "$c\*$p\*$a";
				my @a=split(/\*/,$key);
				my $num_samples=$#samples;
				print OUT "$igv{$c}{$p}{$a}\t$a[0]\t$a[1]\t$a[2]\t";
				print OUT "$ref{$c}{$p}{$a}\t$a[3]\t$a[4]\t";
				my $values=$#{$sample_info{$c}{$p}{$a}};
				print OUT join ("\t", @{$sample_info{$c}{$p}{$a}});
				if ($values != $num_samples){
					for(my $k=$values+1; $k<$num_samples;$k++){
						print OUT "\tn/a\tn/a";
					}
					print OUT "\tn/a\tn/a";
				}
				print OUT "\t$snpeff{$c}{$p}{$a}[$rows]\t";	
				print OUT "\n";
			}
		}
	}		
}	
close OUT;
print "$out is generated\n";    
	    
	    
	    
