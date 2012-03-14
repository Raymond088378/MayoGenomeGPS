	## perl script to generate the tsv (comma seperaetd files for each exome run)
	
	use strict;
	#use warnings;
	use Getopt::Std;

	our ($opt_r,$opt_p);
	print "RAW paramters: @ARGV\n";
	getopt('rp');
	if ( (!defined $opt_r) && (!defined $opt_p) )	{
		die ("Usage: $0 \n\t -r [ un info file] \n\t -p [ ath to input folder] \n");
	}
	else	{
		my $run_info =  $opt_r;chomp $run_info;	
		my $path = $opt_p;chomp $path;
		my $dest=$path."/SampleStatistics.tsv";
		open OUT, ">$dest" or die "can not open $dest : $!\n";
		my @line;
		@line=split(/=/,`perl -ne "/^ANALYSIS/ && print" $run_info`);
		my $analysis=$line[$#line];chomp $analysis;
		@line=split(/=/,`perl -ne "/^SAMPLENAMES/ && print" $run_info`);
		my $sampleNames=$line[$#line];chomp $sampleNames;
		my @sampleArray = split(/:/,$sampleNames);
		@line=split(/=/,`perl -ne "/^VARIANT_TYPE/ && print" $run_info`);
		my $variant_type=$line[$#line];chomp $variant_type;
		@line=split(/=/,`perl -ne "/^SNV_CALLER/ && print" $run_info`);
		my $SNV_caller=$line[$#line];chomp $SNV_caller;
		# @line=split(/=/,`perl -ne "/^NUM_SAMPLES/ && print" $run_info`);
		# my $num_samples=$line[$#line];chomp $num_samples;
		my $num_samples=`echo $sampleNames | tr ":" "\n" | wc -l`;
		$analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`;chomp $analysis;
		$variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`;chomp $variant_type;
		@line=split(/=/,`perl -ne "/^TOOL_INFO/ && print" $run_info`);
		my $tool_info=$line[$#line];chomp $tool_info;
		@line=split(/=/,`perl -ne "/^dbSNP_SNV_rsIDs/ && print" $tool_info`);
		my $dbsnp_file=$line[$#line];chomp $dbsnp_file;
		$dbsnp_file =~ m/.+dbSNP(\d+)/;
		my $dbsnp_v = $1;
		my @To_find;
		print "ANALYSIS:$analysis\n";
		print "SAMPLES:$sampleNames\n";
		print "VARIANT:$variant_type\n";
		print "TOOL_INFO:$tool_info\n";
		# function for adding commas
		sub CommaFormatted
		{
			my $delimiter = ','; # replace comma if desired
			my($n,$d) = split /\./,shift,2;
			my @a = ();
			while($n =~ /\d\d\d\d/)
			{
				$n =~ s/(\d\d\d)$//;
				unshift @a,$1;
			}
			unshift @a,$n;
			$n = join $delimiter,@a;
			$n = "$n\.$d" if $d =~ /\d/;
			return $n;
		}
		# row header
		if ( $analysis eq 'alignment' )	{
			@To_find=("Total Reads","Mapped Reads","Mapped Reads (Q >= 20)");
		}		
		elsif ( ( $analysis eq 'mayo' ) || ( $analysis eq 'variant' ) ) {
			@To_find=("Total Reads","Mapped Reads","Percent duplication","Realigned Reads (Q >= 20)","Mapped Reads (Coding Region)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","Coding SNVs","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","Coding INDELs","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5","Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
		}
		my %sample_numbers=();
		my $uniq;
		# storing all the numbers in a Hash per sample (one hash)
		# print OUT "samples";
		for(my $k = 0; $k < $num_samples;$k++)	
		{
			print OUT "\t$sampleArray[$k]";
			my $file="$path/numbers/$sampleArray[$k].out";
			open SAMPLE, "<$file", or die "could not open $file : $!";
			print "reading numbers from $sampleArray[$k]\n";
			my $id=0;

			while(my $l = <SAMPLE>)	
			{			
				chomp $l;
				if ( $l !~ /^\d/)	{
					$uniq = $id;
					$id++;
				}	
				else	{
					push (@{$sample_numbers{$uniq}},$l);
				}
			}	
			close SAMPLE;
		}
		print OUT "\n";
	#printing the statistics for each sample
		foreach my $key (sort {$a <=> $b} keys %sample_numbers)	{
		print OUT "$To_find[$key]";
		if ( $key eq '1')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_mapped \%)";
				}
			}	
		if ( $key eq '39')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_deleted = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / 3200000000) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print $per_deleted \%)";
				}
			}	
		if ( $key eq '40')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_duplicated = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / 3200000000) * 100);
					# my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_duplicated \%)";
				}
			}	
		if ( ( $key eq '3' ) )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_mapped \%)";
				}
			}
		}		
		if ($analysis ne 'annotation' )	{
			for ( my $c=0; $c < $num_samples; $c++ )	{
				if ( ( $key eq '1') || ( $key eq '3' ) || ( $key eq '39' ) || ( $key eq '40' ) )	{
				}
				else	{
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print";	
				}
			}
		}
		print OUT "\n";
	}	
	undef %sample_numbers;
	close OUT;
	
}	
	
	
		
	
