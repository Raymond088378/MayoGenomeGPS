#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   getSampleInfo.pl

=head1 SYNOPSIS

    USAGE: getSampleInfo.pl -i=/input_dir/primary -e=fastq.gz -o=/path/to/output_dir

=head1 OPTIONS


B<--input, -i>
	Input directory where all fastq or bam files are located.

B<--ext, -e>
	File extension eg: fastq.gz or bam

B<--output, -o>
	Output directory where all configuration files should go	

B<--help,-h>

=head1 DESCRIPTION
	Get info to use in sample_info.txt file and in run_info file realted to
	samples.

=head1 INPUT
	Input dir, file extension and output directory

=head1 OUTPUT
	Output run_info and sample_info files.

=head1 VERSION
	0.1.0

=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./getSampleInfo.pl -i=/input/dir/primary -e=fastq.gz -o=/path/to/outputdir

=cut

use strict;
use warnings;
#use Data::Dumper;
#use File::Basename;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'ext|e=s',
						  'output|o=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

opendir (DIR, $options{input}) or die "Could not open dir $!\n";

my $sample = "";
my $lane = "";
my $index = "";
my $sample_hash;

open S_INFO , ">$options{output}/sample_info.txt" or die "can not open $options{output}/sample_info.txt : $! \n";
open R_INFO , ">$options{output}/run_info.txt" or die "can not open $options{output}/run_info.txt : $! \n";
while (my $file = readdir(DIR)) {
	next if ($file !~ /$options{ext}$/);
	
	# expecting file in following format
	# some_unique_id.FLOWCELLID_LANE_INDEX.EXT
	my @bits = split(/\./, $file);
	# in case a sample sample was run multiple times
	# use id and lane number as sample name
	$bits[0] =~ s/^s_//;
	$bits[0] =~ s/L\d$//;

	my $key = "s_" . $bits[0];
	$key =~ s/[^a-zA-Z0-9_-]*//g;
	push @{$sample_hash->{$key}}, sampleArray($file); # .= $file ."\t";
}

#print "\n";
# print sample info for sample_info file
foreach my $key (keys %{$sample_hash}){
	$sample .= $key .":";
	$lane .= $sample_hash->{$key}[0]->{lane}.":";
	$index .= $sample_hash->{$key}[0]->{index}.":";

	my $name = "";

	## sort so R1 is always before R2
	my @sorted = sort{ $a->{file} cmp $b->{file} } @{$sample_hash->{$key}};
	foreach my $idx (@sorted) {
		$name .= $idx->{file} ."\t";
	}

	$name =~ s/\t$//;
	my $tag=$options{ext};
	
	$tag =~  s/.gz//g;
	print S_INFO uc($tag) .":";

	print S_INFO $key. "=" . $name."\n";
}

# remove all -/_ from sample
$sample =~ s/:$//; # remove last :
$lane =~ s/:$//;
$index =~ s/:$//;

#print sample info for run_info file
#print "\n";
my @parameters=split(/\//,$options{input});
my $delivery = $options{input} =~ s/primary/secondary/g;
print R_INFO "TOOL=GENOME_GPS\n" . "VERSION=1.2\n" . "TYPE=\n" . "DISEASE=NA\n" . "READLENGTH=\n" . "PAIRED=\n" . "ANALYSIS=\n" . "PI=$parameters[3]\n"
. "MULTISAMPLE=\n" . "INPUT_DIR=$options{input}\n" . "BASE_OUTPUT_DIR=/data2/bsi/secondary/\n" . "SAMPLENAMES=" .$sample. "\n" . "GROUPNAMES=\n" . "LANEINDEX=" .$lane. "\n".
"LABINDEXES=" .$index. "\n" . "CHRINDEX=1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:M\n" . "TOOL_INFO=\n" 
. "SAMPLE_INFO=$options{output}/sample_info.txt\n"  . "MEMORY_INFO=\n" . "OUTPUT_FOLDER=$parameters[4]\n" .
"GENOMEBUILD=hg19\nALIGNER=NOVOALIGN\nFASTQC=NO\nFOLDER_FASTQC=/data2/bsi/reports/$parameters[4]/fastqc\nVARIANT_TYPE=BOTH\nSNV_CALLER=GATK\nSOMATIC_CALLER=SOMATICSNIPER\n"
. "SAMPLEINFORMATION=\n" . "DELIVERY_FOLDER=$delivery\n"
. "TERTIARY_FOLDER=/data2/bsi/tertiary/$parameters[3]/<analsyis type>/$parameters[4]";

print "\nPlease fill these columns in run info file: TYPE,READLENGTH,PAIRED,ANALYSIS,MULTISAMPLE,TOOL_INFO,MEMORY_INFO,SAMPLEINFORMATION \nNOTE:\nFor Standard run user should copy the tool information and memory information file from the /path/to/config folder of the scripts\nValidate the configuration files again before running the workflow\n\n";
exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input ext output);

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}
}

#############################################################################
sub sampleArray {
	my $file = shift;

	my $obj = ();
	my $tag=$options{ext};
	my @bits = split(/\./, $file);
	my @flow = split(/_/, $bits[1]);

	$obj->{'file'} = $file;
	$obj->{'lane'} = substr($flow[1], 1);

	if ($tag =~ /bam/)	{
		if (scalar(@flow) == 4) {
			$obj->{'index'} = substr($flow[3], 1);
		} elsif(scalar(@flow) == 3) {
			$obj->{'index'} = substr($flow[2], 1);
		}
		else	{
			$obj->{'index'} = "-";
		}
	}else	{
		if (scalar(@flow) == 4) {
			$obj->{'index'} = substr($flow[3], 1);
		}
		else	{
			$obj->{'index'} = "-";
		}
	}
		
	return $obj;
}
