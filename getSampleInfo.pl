#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   getSampleInfo.pl

=head1 SYNOPSIS

    USAGE: getSampleInfo.pl -i=/input_dir/primary -e=fastq.gz -t=rna

=head1 OPTIONS


B<--input, -i>
	Input directory where all fastq or bam files are located.

B<--ext, -e>
	File extention eg: fastq.gz or bam

B<--help,-h>

=head1 DESCRIPTION
	Get info to use in sample_info.txt file and in run_info file realted to
	samples.

=head1 INPUT
	Input dir, file extention and type of workflow

=head1 OUTPUT

=head1 VERSION
	0.1.0

=head1  CONTACT
  bjaysheel@gmail.com


==head1 EXAMPLE
	./getSampleInfo.pl -i=/input/dir/primary -e=fastq.gz 

=cut

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'ext|e=s',
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

while (my $file = readdir(DIR)) {
	next if ($file !~ /$options{ext}$/);
	next if ($file !~ /$options{type}/i);

	# expecting file in following format
	# some_unique_id.FLOWCELLID_LANE_INDEX.EXT
	my @bits = split(/\./, $file);

	# in case a sample sample was run multiple times
	# use id and lane number as sample name
	$bits[0] =~ s/^s_//;
	$bits[0] =~ s/L\d$//;

	my $key = $bits[0];
	$key =~ s/_|-//g;

	push @{$sample_hash->{$key}}, sampleArray($file); # .= $file ."\t";
}

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

	if ($options{type} !~ /rna/i) {
		print uc($options{ext}) .":";
	}

	print $key. "=" . $name."\n";
}

# remove all -/_ from sample
$sample =~ s/:$//; # remove last :
$lane =~ s/:$//;
$index =~ s/:$//;

#print sample info for run_info file
print "SAMPLENAMES=" .$sample. "\n";
print "LANEINDEX=" .$lane. "\n";
print "LABINDEXES=" .$index. "\n";

exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = qw(input ext);

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

	my @bits = split(/\./, $file);
	my @flow = split(/_/, $bits[1]);

	$obj->{'file'} = $file;
	$obj->{'lane'} = substr($flow[1], 1);

	if (scalar(@flow) == 4) {
		$obj->{'index'} = substr($flow[3], 1);
	} else {
		$obj->{'index'} = substr($flow[2], 1);
	}

	return $obj;
}
