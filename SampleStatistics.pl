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
		@line=split(/=/,`perl -ne "/^MULTISAMPLE/ && print" $run_info`);
		my $multi=$line[$#line];chomp $multi;
		@line=split(/=/,`perl -ne "/^TYPE/ && print" $run_info`);
		my $tool=$line[$#line];chomp $tool;
		$tool=lc($tool);
		@line=split(/=/,`perl -ne "/^LANEINDEX/ && print" $run_info`);
		my $lanes=$line[$#line];chomp $lanes;
		my @laneArray = split(/:/,$lanes);
        @line=split(/=/,`perl -ne "/^LABINDEXES/ && print" $run_info`);
		my $indexes=$line[$#line];chomp $indexes;
		my @IndexArray = split(/:/,$indexes);
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
		# header name
	$multi =~ s/\s+$//;
	$analysis =~ s/\s+$//;
	if ( ( $analysis ne 'alignment' ) && ( $analysis ne 'annotation' ) ) {
		if ($tool eq 'whole_genome')	{
			if ($analysis eq 'variant')	{
				@To_find=("Total Reads","Mapped Reads","Mapped Reads(CodingRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5","Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
			}
			else	{
				@To_find=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(CodingRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5","Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
			}
		}
		elsif ($tool eq 'exome')	{
			if ($analysis eq 'realignment' || $analysis eq 'external' || $analysis eq 'mayo' || $analysis eq 'realign-mayo' )	{
				@To_find=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(in CaptureRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","INDELs in CaptureRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
			}
			else	{
				@To_find=("Total Reads","Mapped Reads","Mapped Reads(in CaptureRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","INDELs in CaptureRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
			}
		}	
	}
	if ($analysis eq 'alignment')	{
		@To_find=("Total Reads","Mapped Reads","Percent duplication");
	}

	if ($analysis eq 'annotation')	{
		if ($variant_type eq 'BOTH' )	{
			@To_find=("Total SNVs","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs" ,"In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
		}
		elsif ($variant_type eq 'SNV')	{
			@To_find=("Total SNVs","Total SNVs(Known)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs(Novel)","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5");
		}	
		elsif ($variant_type eq 'INDEL')	{
			@To_find=("Total INDELs" ,"In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
		}	
	}
	
		my %sample_numbers=();
		my $uniq;
		# storing all the numbers in a Hash per sample (one hash)
		# print OUT "samples";
		print OUT "SampleNamesUsed/info";
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
		print OUT "lanes";
        for(my $k = 0; $k < $num_samples;$k++)  {
            print OUT "\t$laneArray[$k]";
        }
        print OUT "\n";
        print OUT "indexes";
        for(my $k = 0; $k < $num_samples;$k++)  {
            print OUT "\t$IndexArray[$k]";
        }
        print OUT "\n";
        
		
		
	#printing the statistics for each sample
		foreach my $key (sort {$a <=> $b} keys %sample_numbers)	{
		print OUT "$To_find[$key]";
		if ( $key eq '1' && $analysis ne 'annotation')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_mapped \%)";
				}
			}	
		if ($analysis eq 'variant' )	{
			if ( $key eq '37')	{
					for (my $c=0; $c < $num_samples;$c++)	{
						my $per_deleted = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{2}}[$c]) * 100);
						my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						#my $print=${$sample_numbers{$key}}[$c];
						print OUT "\t$print $per_deleted \%)";
					}
				}	
			if ( $key eq '38')	{
					for (my $c=0; $c < $num_samples;$c++)	{
						my $per_duplicated = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{2}}[$c]) * 100);
						# my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						my $print=${$sample_numbers{$key}}[$c];
						print OUT "\t$print ($per_duplicated \%)";
					}
				}	
		}
		if ($key eq '2' && $analysis ne 'annotation' && $analysis ne 'variant' )	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $print=sprintf("%.2f",$sample_numbers{$key}[$c]);
				print OUT "\t$print\%";
			}	
		}
		if ($key eq '2' && $analysis eq 'variant')	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print ($per_mapped \%)";
			}
		}	
		
		if ($analysis ne 'variant' )	{
			if ( $key eq '39')	{
					for (my $c=0; $c < $num_samples;$c++)	{
						my $per_deleted = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{4}}[$c]) * 100);
						my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						#my $print=${$sample_numbers{$key}}[$c];
						print OUT "\t$print $per_deleted \%)";
					}
				}	
			if ( $key eq '40')	{
					for (my $c=0; $c < $num_samples;$c++)	{
						my $per_duplicated = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{4}}[$c]) * 100);
						# my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						my $print=${$sample_numbers{$key}}[$c];
						print OUT "\t$print ($per_duplicated \%)";
					}
				}	
		}	
		if ( ( $key eq '3' ) && ($analysis ne 'variant') )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_mapped \%)";
				}
			}
		}	
		if ( ( $key eq '3' ) && ($analysis eq 'variant') )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					#my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ";
				}
			}
		}	
		if ( ( $key eq '4' ) && ($analysis ne 'variant') )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ($per_mapped \%)";
				}
			}
		}	
		if ( ( $key eq '4' ) && ($analysis eq 'variant') )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					#my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					#my $print=${$sample_numbers{$key}}[$c];
					print OUT "\t$print ";
				}
			}
		}			
		if ($analysis ne 'annotation' )	{
			for ( my $c=0; $c < $num_samples; $c++ )	{
				if ($analysis ne 'variant')	{
					if ( ( $key ne '1') && ( $key ne '3' ) && ($key ne '2') && ($key ne '4') && ( $key ne '39' ) && ( $key ne '40' ) )	{
						my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						print OUT "\t$print";	
					}
				}
				elsif ($analysis eq 'variant')	{
					if ( ( $key ne '1') && ( $key ne '3' ) && ($key ne '2') && ($key ne '4') && ( $key ne '37' ) && ( $key ne '38' ) )	{
						my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
						print OUT "\t$print";	
					}
				}
			}
		}
		else	{
			for ( my $c=0; $c < $num_samples; $c++ )	{
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print";	
			}
		}		
		print OUT "\n";
	}	
	undef %sample_numbers;
	close OUT;
	
}	
	
	
		
	
