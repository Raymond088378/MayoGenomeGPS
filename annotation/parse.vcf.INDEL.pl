#!/usr/local/biotools/perl/5.10.0/bin/perl
	
	## parse VCF to create old GATK format file
	##06/13/2011
	## baheti.saurabh@mayo.edu
	## GTAK updated the format of the VCF so needs and update to parsing script
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sampleA
	#chr1    12811188        .       G       GC      252.31  .       AC=2;AF=1.00;AN=2;Dels=0.00;HRun=3;HaplotypeScore=0.0000;MQ=43.31;MQ0=0;QD=31.54;SB=-83.88
	#     GT:AD:DP:GQ:PL  1/1:1,7:8:21.07:294,21,0
	use strict;
	#use warnings;
	use Getopt::Std;

	our($opt_i, $opt_o, $opt_s, $opt_h);
	#print "INFO: script to parse the VCF indel to output old GATK format\n";
	#print "RAW paramters: @ARGV\n";
	getopt('iosh');
	if ( (!defined $opt_i) && (!defined $opt_o) && (!defined $opt_s) && (!defined $opt_h) )	{
		die ("Usage: $0 \n\t-i [nput vcf file] \n\t-o [utput gatk format] \n\t-s[ample name] \n\t -h [ header flag]\n");
	}		
	else	{
		my $source=$opt_i;
		my $output=$opt_o;
		my $sample=$opt_s;chomp $sample;
		open FH, "<$source" or die "can not open $source :$! \n";
		open OUT, ">$output" or die "can not open $output : $! \n";
		my $flag=$opt_h;
		if ($flag ==1)	{
			print OUT "\t" x 8 . "$sample" . "\t" . "\n";
			print OUT "Chr\tStart\tStop\tInCaptureKit\t#AlternateHits\tRef\tAlt\tBase-Length\tIndel-supportedRead\tReadDepth\n";
		}
		my ($header_len,$sample_col,$info_col,$format_col,$chr,$ReadDepth,$GenoType,$AllelicDepth,$format_len);
		my (@alt_reads,@format_data,@sample_data);
		while( my $l = <FH>)	{
			chomp $l;
			next if ( $l =~ /^##/ );   ## skipping header of the VCF format
			my @call = split(/\t/,$l);
			if( $l =~ /^#/)	{		## reading the header to get header information
				$header_len=scalar(@call);
				for( my $i = 0; $i < $header_len ; $i++ )	{	## to get the position of sample specific data	
					if ($call[$i] eq 'INFO')	{
						$info_col=$i;
					}
					if($call[$i] eq 'FORMAT')	{
						$format_col=$i;
					}	
					if($call[$i] eq $sample)	{
						$sample_col=$i;
					}	
				}		
				next;	
			}	
			next if ((length($call[3]) == 1) && (length($call[4]) == 1) );
            
			$chr=$call[0];
			@format_data = split(/:/,$call[$format_col]);
			@sample_data = split(/:/,$call[$sample_col]);
			$format_len=scalar(@format_data);
			for( my $i = 0; $i < $format_len ; $i++ )	{			## to get genotype and allelic depth columns
				if($format_data[$i] eq 'GT')	{
					$GenoType = $i;
				}
				if($format_data[$i] eq 'AD')	{
					$AllelicDepth = $i;
				}
			}	
			@alt_reads=split(/,/,$sample_data[$AllelicDepth]);
            
			if ($sample_data[$GenoType] ne './.')	{
				# to get INS or DEL information and bases
				my $ref_len=length($call[3]);
				my $alt_len=length($call[4]);
				my $INDEL;
				my $Bases;
				my $Start;
				my $Stop;
				if ( $ref_len < $alt_len )	{
					$INDEL = '+';
					$Bases = length(substr( $call[4],1,$alt_len ));
					$Start = $call[1];
					$Stop = $Start;
				}
				else	{
					$INDEL = '-';
					$Bases=length(substr( $call[3],1,$ref_len ));
					$Start = $call[1];
					$Stop = $Start + $Bases;
				}
				$ReadDepth=$alt_reads[0]+$alt_reads[$#alt_reads];
				my $capture_flag=1;
				if (grep /CAPTURE=0/,$call[$info_col])	{
					$capture_flag=0;
				}
				my $reg=$call[$info_col];
				my $multi="-";
				$multi=$2 if ($reg =~ /(\w*)ED=(\d*)/);
				if (length($multi) == 0)	{
					$multi=-1;
				}	
				if ($sample_data[$GenoType] eq '.')	{
					print OUT "$chr\t$Start\t$Stop\t$capture_flag\t$multi\t$call[3]\t$call[4]\t$Bases\t-\t-\n";		
				}
				else	{
					print OUT "$chr\t$Start\t$Stop\t$capture_flag\t$multi\t$call[3]\t$call[4]\t$Bases\t$alt_reads[$#alt_reads]\t$ReadDepth\n";
				}
			}
			
		}
		close FH;
		close OUT;
	}		
## END of script		
