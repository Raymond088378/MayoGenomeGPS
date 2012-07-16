#script to parse VCF format variant calls to readable format SNVs
#
#chr1	109	.	A	T	328.94	.	AC=1;AF=0.50;AN=2;DP=218;DS;Dels=0.00;HRun=0;HaplotypeScore=130.5110;MQ=19.84;MQ0=52;QD=1.51;SB=54.12;sumGLbyD=1.65	GT:AD:DP:GQ:PL	0/1:140,68:65:99:359,0,1222
#chr1	147	.	C	A	214.51	.	AC=1;AF=0.50;AN=2;DP=99;DS;Dels=0.00;HRun=2;HaplotypeScore=68.1006;MQ=18.83;MQ0=23;QD=2.17;SB=-54.86;sumGLbyD=2.47	GT:AD:DP:GQ:PL	0/1:76,23:27:99:244,0,471
#chr1	180	.	T	C	173.38	.	AC=1;AF=0.50;AN=2;DP=190;Dels=0.02;HRun=2;HaplotypeScore=130.1643;MQ=19.80;MQ0=38;QD=0.91;SB=-16.96;sumGLbyD=1.07	GT:AD:DP:GQ:PL	0/1:153,34:52:99:203,0,1141
# Assumptions: DP is read depth and AD is number of ref and alt call.

use strict;
use warnings;
use Getopt::Std;

our($opt_i, $opt_o, $opt_s, $opt_h);
#print "INFO: script to parse the VCF indel to output old GATK format\n";
#print "RAW paramters: @ARGV\n";
getopt('iosh');
if ( (!defined $opt_i) && (!defined $opt_o) && (!defined $opt_s)  && (!defined $opt_h) )	{
	die ("Usage: $0 \n\t-i [nput vcf file] \n\t-o [utput gatk format] \n\t-s[ample name] \n\t-h [ header boolean]\n");
}		
else	{
	my $source=$opt_i;
	my $output=$opt_o;
	my $sample=$opt_s;chomp $sample;
	my $flag=$opt_h;
	open FH, "<$source" or die "can not open $source :$! \n";
	open OUT, ">$output" or die "can not open $output : $! \n";		
	if ($flag == 1)	{
		print OUT "\t" x 6 . "$sample" . "\t" x 5 . "\n";
		print OUT "Chr\tPosition\tInCaptureKit\t#AlternateHits\tRef\tAlt\tGenotypeClass\tAlt-SupportedReads\tRef-SupportedReads\tReadDepth\tQuality\tCloseToIndel\n";
	}
	
	my ($header_len,$sample_col,$qual_col,$info_col,$format_col,$chr,$ReadDepth,$GenoType,$AllelicDepth,$format_len,$pos_col,$ref_col,$alt_col,$chr_col,$Genoqual);
	$Genoqual=0;
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
				if($call[$i] eq 'QUAL')	{
					$qual_col=$i;
				}	
				if($call[$i] eq 'POS')	{
					$pos_col=$i;
				}
				if($call[$i] eq 'REF')	{
					$ref_col=$i;
				}
				if($call[$i] eq 'ALT')	{
					$alt_col=$i;
				}
				if($call[$i] eq '#CHROM')	{
					$chr_col=$i;
				}
			}		
			next;	
		}	
		next if ($call[$sample_col] eq './.');
		next if (length($call[$ref_col]) > 1);
		next if (length($call[$alt_col]) > 1);
		@format_data = split(/:/,$call[$format_col]);
		@sample_data = split(/:/,$call[$sample_col]);
		$format_len=scalar(@format_data);
		my $quality;
		for( my $i = 0; $i < $format_len ; $i++ )	{			## to get genotype and allelic depth columns
			if($format_data[$i] eq 'GT')	{
				$GenoType = $i;
			}
			if($format_data[$i] eq 'AD')	{
				$AllelicDepth = $i;
			}
			if ($format_data[$i] eq 'GQ')	{
				$Genoqual = $i;
	    
			}
		}	
		if ($Genoqual == 0) {
			$quality=99;
		}
		else    {
			$quality=$sample_data[$Genoqual];
		}
		@alt_reads=split(/,/,$sample_data[$AllelicDepth]);
		$ReadDepth=$alt_reads[0]+$alt_reads[$#alt_reads];
		#CAPTURE=1
		my $close2indel_flag=0;
		if (grep /CLOSE2INDEL=1/,$call[$info_col])	{
			$close2indel_flag=1;	
		}
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
		if ($sample_data[$GenoType] eq '.' )  {
			print OUT "$call[$chr_col]\t$call[$pos_col]\t$capture_flag\t$multi\t$call[$ref_col]\t$call[$alt_col]\t-\t$alt_reads[$#alt_reads]\t$alt_reads[0]\t$ReadDepth\t$quality\t$close2indel_flag\n";
		}	
		if ($sample_data[$GenoType] eq '0/0' )  {
		        print OUT "$call[$chr_col]\t$call[$pos_col]\t$capture_flag\t$multi\t$call[$ref_col]\t$call[$alt_col]\t$call[$ref_col]$call[$ref_col]\t$alt_reads[$#alt_reads]\t$alt_reads[0]\t$ReadDepth\t$quality\t$close2indel_flag\n";
		}
		
		if ($sample_data[$GenoType] eq '0/1' )	{
			print OUT "$call[$chr_col]\t$call[$pos_col]\t$capture_flag\t$multi\t$call[$ref_col]\t$call[$alt_col]\t$call[$ref_col]$call[$alt_col]\t$alt_reads[$#alt_reads]\t$alt_reads[0]\t$ReadDepth\t$quality\t$close2indel_flag\n";	
		}
		if ($sample_data[$GenoType] eq '1/1' )	{
			print OUT "$call[$chr_col]\t$call[$pos_col]\t$capture_flag\t$multi\t$call[$ref_col]\t$call[$alt_col]\t$call[$alt_col]$call[$alt_col]\t$alt_reads[$#alt_reads]\t$alt_reads[0]\t$ReadDepth\t$quality\t$close2indel_flag\n";	
		}
	}	
	close FH;
	close OUT;
}		
## END of script		
	
