#!/usr/local/biotools/perl/5.10.0/bin/perl

use strict;
use warnings;
die "Usage:$0 <file w/ list of files> <flag><output file>\n" if ($#ARGV != 2);
my $list=$ARGV[0];
open FH, "$list" or die "can not open $list : $!\n";
my $flag=$ARGV[1];
my $out=$ARGV[2];
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
my (@sample,@annot_snpeff);
my (%igv,%ref, %snpeff, %sample_info);
my ($head1,$head3,$head4);
my $count=0;
my $prev_samples=0;
my %chrvalue = ("chrX"=>23,"chrY"=>24,"chrM"=>25);
for (my $i=1; $i<23; $i++) {
	$chrvalue{"chr".$i} = $i;
}
my $a1;
while(my $l = <FH>){
   	chomp $l;
	open FILE, "$l" or die "";
	my $prev=0;
   	print "Reading $l\n";
	$START_INFO=12;
	my $num_samples=0;
	while(my $k = <FILE>){
		chomp $k;
		my @a = split("\t",$k);
		my $last_col=$#a;
		if ( $. == 1){
			for (my $j=0; $j <= $last_col; $j++)	{
				if ($a[$j] =~ m/^SNPEFF/)	{
					last;
				}
				elsif ($a[$j] =~ m/^\S/ && $a[$j] !~ m/^Allele/ && $a[$j] !~ m/^-/)	{		
					$samples[$i]=$a[$j];
					$num_samples++;
					$i++;
					$j++;
				}	
			}
		}
		elsif ( $. == 2){
			$last_col=$#a;
			$STOP_INFO=$START_INFO+(2*$num_samples)-1;
			$SNPEFF_START=$STOP_INFO+1;
			@annot_snpeff=($SNPEFF_START .. $last_col);
			$a1=$last_col-$SNPEFF_START;
			$head1=join("\t",@a[@annot_ref]);
			$head3=join("\t",@a[@annot_snpeff]);
			for (my $k=$START_INFO;$k<$STOP_INFO;)	{
				@sample=($k .. $k+1);
				$head4=join("\t",@a[@sample]);
				$k+=2;
			}	
			$prev_samples=$i-$num_samples;
		}
		else	{
			$last_col=$#a;
			my $id="$a[$CHR]\*$a[$START]\*$a[$STOP]\*$a[$ALT]\*$a[$BASE]";
			my $id_a="$a[$STOP]\*$a[$ALT]\*$a[$BASE]";
			##snpeff values	
			#if ($count == 2)	{
			#	<STDIN>;
			#}
			my $snpeff_value=join("\t",@a[@annot_snpeff]);
			push(@{$snpeff{$a[$CHR]}{$a[$START]}{$id_a}},$snpeff_value);
			if ($id ne $prev)	{
				my $value=$#{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}};
				#print "$a[$START]\tvalue=$value\tpre=$prev_samples\n";
				#<STDIN>;#$value=$prev_samples;
				$value=$value+1;	
				#if ($a[$START] == "1684347" || $a[$START] == "761957")	{print "value=$value\ti=$i\tpre=$prev_samples\n";}
				my $sample_value;
				my $sam;
				if ($value == $prev_samples){
					for (my $k=$START_INFO;$k<$STOP_INFO;)	{
						@sample=($k .. $k+1);
						$sample_value=join("\t",@a[@sample]);
						push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sample_value);
						$k+=2;
					}
				}
				else {
					$value=$prev_samples-$value;
					for(my $ll=0; $ll <$value;$ll++){
						$sam="n/a\tn/a\t";
						$sam =~ s/\s*$//;
						push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sam);	
					}
					for (my $k=$START_INFO;$k<$STOP_INFO;)	{
						@sample=($k .. $k+1);
						$sample_value=join("\t",@a[@sample]);
						push(@{$sample_info{$a[$CHR]}{$a[$START]}{$id_a}},$sample_value);	
						$k+=2;
					}	
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
	#<STDIN>;
}
close FH;
print OUT "-\t";
print OUT "\t" x $BASE;
print OUT join ("\t\t",@samples);
print OUT "\t\tSNPEFF Annotation";
print OUT "\t" x $a1;
print OUT "\n";
print OUT "IGV Link\tChr\tStart\tStop\t$head1\tAlt\tBase-Length\t";
for(my $i=0; $i <=$#samples; $i++)	{
	print OUT "$head4\t";
}
print OUT "$head3\n";	
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
print "Merging all the sample files\n";
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
				my $val=join("\t",@{$sample_info{$c}{$p}{$a}});
				my @val1=split('\s+',$val);
				my $values=@val1;
				$values=$values/2;
				print OUT join ("\t", @{$sample_info{$c}{$p}{$a}});
				if ($values != $num_samples+1){
					for(my $k=$values; $k<=$num_samples;$k++){
						print OUT "\tn/a\tn/a";
					}
				}
				print OUT "\t$snpeff{$c}{$p}{$a}[$rows]";	
				print OUT "\n";
			}
		}
	}		
}	
close OUT;
print "$out is generated\n";    
	    
	    
	    
