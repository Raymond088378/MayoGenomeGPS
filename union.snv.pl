#!/usr/local/biotools/perl/5.10.0/bin/perl

use strict;
use warnings;

die "Usage:$0 <file w/ list of files> <flag> <output file>\n" if ($#ARGV != 2); 
my $list=$ARGV[0];
my $flag=$ARGV[1];
my $out=$ARGV[2];
open FH, "$list" or die "can not open $list :$!\n";
open OUT, ">$out" or die "can not open $out :$!\n";
my $IGV=0;
my $CHR=1;
my $POS=2;
my $ALT=19;
my $dbSNP=3;
my $REF=18;
my ($START_INFO,$STOP_INFO,$SNPEFF_START,$SNPEFF_STOP,$SIFT_START,$SIFT_STOP);
my $num_col_samples;
my %chrvalue = ("chrX"=>23,"chrY"=>24,"chrM"=>25);
for (my $i=1; $i<23; $i++) {
	$chrvalue{"chr".$i} = $i;
}
my @samples;  ## store the sample names and print in the same order
my $i =0;
my @annot_ref=($dbSNP .. $REF);
my (@annot_sift,@annot_snpeff,@sample);
my (%igv,%ref, %sift, %snpeff, %sample_info);
my ($head1,$head2,$head3,$head4);
my $count=0;
my $a3;
my $prev_samples=0;
while(my $l = <FH>){
    chomp $l;
    open FILE, "$l" or die "can not open $l : $!\n";
	my $prev=0;
	print "Reading $l\n";
	$START_INFO=20;
	my $num_samples=0;
	while(my $k = <FILE>){
		chomp $k;
		my @a = split("\t",$k);
		my $last_col=$#a;
		if ( $. == 1){
			for (my $j=0; $j <= $last_col; $j++)	{
				if ($a[$j] =~ m/^SIFT/)	{
					last;
				}
				elsif ($a[$j] =~ m/^\S/ && $a[$j] !~ m/^Allele/ && $a[$j] !~ m/^-/)	{		
					$samples[$i]=$a[$j];
					$j++;
					$num_samples++;
					$i++;
				}
			}
		}
		elsif ( $. == 2){
			my $last_col=$#a;
			$STOP_INFO=$START_INFO+(7*$num_samples)-1;
			$num_col_samples=$STOP_INFO -$START_INFO +1;
			$SIFT_START=$STOP_INFO+1;
			$SIFT_STOP=$STOP_INFO+1+30;
			@annot_sift=($SIFT_START .. $SIFT_STOP);
			$SNPEFF_START=$SIFT_STOP+1;
			@annot_snpeff=($SNPEFF_START .. $last_col);
			$a3=$last_col-$SNPEFF_START;
			for (my $k=$START_INFO;$k<$STOP_INFO;)	{
				@sample=($k .. $k+6);
				$head4=join("\t",@a[@sample]);
				$k+=7;
			}	
			$head1=join("\t",@a[@annot_ref]);
			$head2=join("\t",@a[@annot_sift]);
			$head3=join("\t",@a[@annot_snpeff]);
			$prev_samples=$i-$num_samples;
		}
		else	{
			my $last_col=$#a;
			my $id="$a[$CHR],$a[$POS],$a[$ALT]";
			##snpeff values	
			my $snpeff_value=join("\t",@a[@annot_snpeff]);
			push(@{$snpeff{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$snpeff_value);
			if ($id ne $prev)	{
				my $value=$#{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}};
				#if ($flag eq "multi")	{$value=($value+2) * $prev_samples;}
				#else	{$value=($value+2);}	
				$value=$value+1;
				my $sample_value;
				my $sam;
				if ($value == $prev_samples){
					for (my $k=$START_INFO;$k<$STOP_INFO;)	{
						@sample=($k .. $k+6);
						$sample_value=join("\t",@a[@sample]);
						push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sample_value);
						$k+=7;
					}
				}
				else {
					$value=$prev_samples-$value;
					for(my $ll=0; $ll <$value;$ll++){
						$sam="n/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\t";  
						$sam =~ s/\s*$//;
						push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sam);	
					}
					for (my $k=$START_INFO;$k<$STOP_INFO;)	{
						@sample=($k .. $k+5);
						$sample_value=join("\t",@a[@sample]);
						push(@{$sample_info{$a[$CHR]}{$a[$POS]}{$a[$ALT]}},$sample_value);	
						$k+=7;
					}	
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
print OUT join ("\t\t\t\t\t\t\t",@samples);
print OUT "\t\t\t\t\t\t\tSIFT Annotation";
my $a2 = $SIFT_STOP - $SIFT_START +1;
print OUT "\t" x $a2;
print OUT "SNPEFF Annotation";
print OUT "\t" x $a3;
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
				my $val=join("\t",@{$sample_info{$c}{$p}{$a}});
				my @val1=split('\s+',$val);
				my $values=@val1;
				$values=$values/7;
				print OUT join ("\t", @{$sample_info{$c}{$p}{$a}});
				if ($values != $num_samples+1){
					for(my $k=$values; $k<=$num_samples;$k++){
						print OUT "\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a";
					}
				}
				print OUT "\t$sift{$c}{$p}{$a}\t$snpeff{$c}{$p}{$a}[$rows]";	
				print OUT "\n";
			}
		}
	}
}
close OUT;
print "$out is generated\n";	
    
	    
	    
	    
