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
		@line=split(/=/,`perl -ne "/^ALIGNER/ && print" $run_info`);
		my $Aligner=$line[$#line];chomp $Aligner;
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
		my (@align,@snv,@indel,@sv);
		my @names;
		if ($multi eq 'NO')	{
			if ($analysis eq 'realignment' || $analysis eq 'realign-mayo' || $analysis eq 'mayo' || $analysis eq 'external')	{
			if ($tool eq 'exome')	{
				@align=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(in CaptureRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","INDELs in CaptureRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@align,@snv,@indel);
			}	
			elsif ($tool eq 'whole_genome')	{
				@align=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(in CodingRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				
				@sv=("Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
				@names=(@align,@snv,@indel,@sv);
			}		
		}
		elsif ($analysis eq 'variant')	{
			if ($tool eq 'exome')	{
				@align=("Total Reads","Realigned Mapped Reads","Mapped Reads(in CaptureRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","INDELs in CaptureRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@align,@snv,@indel);
			}	
			elsif ($tool eq 'whole_genome')	{
				@align=("Total Reads","Realigned Mapped Reads","Mapped Reads(in CodingRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				
				@sv=("Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
				@names=(@align,@snv,@indel,@sv);
			}
		}
		elsif ($analysis eq 'alignment')	{
			@align=("Total Reads","Mapped Reads","Percent duplication");
			@names=(@align);
		}
		elsif ($analysis eq 'annotation')	{
			if ($variant_type eq 'BOTH' )	{
				@snv=("Total SNVs (${SNV_caller})","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@snv,@indel);
			}
			elsif ($variant_type eq 'SNV')	{
				@snv=("Total SNVs (${SNV_caller})","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@snv);
			}
			elsif ($variant_type eq 'INDEL')	{
				@indel=("Total INDELs (GATK)","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@indel);
			}
		}
		elsif ($analysis eq 'ontarget')	{
			if ($tool eq 'exome')	{
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","INDELs in CaptureRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@snv,@indel);
			}	
			elsif ($tool eq 'whole_genome')	{
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@snv,@indel);
			}
		}
		else	{
			print "incorrect analysis type \n";
			exit 1;
		}	
	}
	elsif ($multi eq 'YES')	{
		if ($analysis eq 'mayo' || $analysis eq 'external' || $analysis eq 'realign-mayo' || $analysis eq 'realignment')	{
			if ($tool eq 'exome')	{
				@align=("Total Reads","Mapped Reads","Percent duplication","Mapped Reads(in CaptureRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","INDELs in CaptureRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@align,@snv,@indel);
			}	
			elsif ($tool eq 'whole_genome')	{
				@align=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(in CodingRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				
				@sv=("Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
				@names=(@align,@snv,@indel,@sv);
			}	
		}
		elsif ($analysis == "variant")	{
			if ($tool eq 'exome')	{
				@align=("Realigned Mapped Reads","Mapped Reads(in CaptureRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME",
				"Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","INDELs in CaptureRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				@names=(@align,@snv,@indel);
			}	
			elsif ($tool eq 'whole_genome')	{
				@align=("Realigned Mapped Reads","Mapped Reads(in CodingRegion)");
				
				@snv=("Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME","Total SNVs","Ti/Tv Ratio","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","START_LOST","STOP_GAINED","STOP_LOST","RARE_AMINO_ACID","NON_SYNONYMOUS_CODING","SYNONYMOUS_START","NON_SYNONYMOUS_START","START_GAINED","SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_STOP","UTR_5_PRIME","UTR_3_PRIME");
				
				@indel=("Total INDELs (GATK)","Filtered INDELs (GATK)","INDELs in CodingRegion","EXON_DELETED","FRAME_SHIFT","CODON_CHANGE","UTR_5_DELETED","UTR_3_DELETED","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","UTR_5_PRIME","UTR_3_PRIME");
				
				@sv=("Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
				@names=(@align,@snv,@indel,@sv);
			}
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
		print OUT "$names[$key]";
		if ( $key eq '1' && $analysis ne 'annotation'  && $analysis ne 'ontarget' )	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print ($per_mapped \%)";
			}
		}
		if ($key eq '2' && $analysis ne 'annotation' && $analysis ne 'variant' && $analysis ne 'ontarget' )	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $print=sprintf("%.2f",$sample_numbers{$key}[$c]);
				print OUT "\t$print\%";
			}		
		}
		elsif ($key eq '2' && $analysis eq 'variant')	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print ($per_mapped \%)";	
			}
		}
		if ($key eq '3' && $analysis ne 'annotation' && $analysis ne 'variant' && $analysis ne 'ontarget' )	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print ($per_mapped \%)";	
			}
		}
		if ($key eq '4' && $analysis ne 'annotation' && $analysis ne 'variant' && $analysis ne 'ontarget' && $multi ne 'YES')	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
				my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
				print OUT "\t$print ($per_mapped \%)";	
			}
		}
		if ($analysis ne 'annotation' && $analysis ne 'ontarget')	{
			for ( my $c=0; $c < $num_samples; $c++ )	{
				if ($multi eq 'NO')	{
					if ($analysis ne 'variant')	{
						if ( ( $key ne '1') && ( $key ne '3' ) && ($key ne '2') && ($key ne '4') )	{
							my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
							print OUT "\t$print";	
						}
					}
					elsif ($analysis eq 'variant')	{
						if ( ( $key ne '1') && ( $key ne '2') )	{
							my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
							print OUT "\t$print";	
						}
					}
				}
				elsif ($multi eq 'YES')	{
					if ($analysis ne 'variant')	{
						if ( ( $key ne '1') && ( $key ne '3' ) && ($key ne '2'))	{
							my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
							print OUT "\t$print";	
						}
					}
					elsif ($analysis eq 'variant')	{
						if ( ( $key ne '1') && ( $key ne '2') )	{
							my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
							print OUT "\t$print";	
						}
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
	
	
		
	
