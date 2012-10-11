#!/usr/local/biotools/perl/5.14.2/bin/perl

=head1 NAME
   vcf_blat_verify.pl

=head1 SYNOPSIS

    USAGE: vcf_blat_verify.pl --input input_vcf_file  --output output_file --reference reference genome --window 50 [--minScore 70 --minidentity 90]

=head1 OPTIONS

B<--input,-i>
   VCF input file

B<--output, -o>
	Output file

B<--reference,-r>
   reference genome file

B<--window, -w>
	windows size to capture length of dna up and down stream from vcf location

B<--blat_path, -b>
	full path to the blat tool

B<--samtools_path, -sam>
	full path to the samtools tool

B<--blat_ref, -br>
	blat reference genome file

B<--minScore, -m>
	Optional sets minimum score.  This is twice the matches minus the
    mismatches minus some sort of gap penalty.  Default is 70

B<--minIdentity, -t>
	Optional Sets minimum sequence identity (in percent).  Default is 90

B<--help,-h>
   This help message

=head1  DESCRIPTION
    Identify uniqueness of vcf given location and window size that covers up and down stream

=head1  INPUT


=head1  OUTPUT
	input VCF with BLAT information added

=head1  CONTACT
  Saurabh.Baheti @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   vcf_blat_verify.pl --input /file/path/filename.vcf --output /file/path/output/filename.vcf --reference /file/path/hg19.fsa --window 50

=cut

use strict;
use warnings;
#use Data::Dumper;
#use Cwd;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use POSIX;

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'samtools_path|sam=s',
						  'blat_path|b=s',
						  'blat_ref|br=s',
						  'reference|r=s',
						  'window|w=s',
						  'minScore|m=s',
						  'minIdentity|t=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

print " Starting main program\n";
timer(); #call timer to see when process ended.

my $blat_ref=$options{blat_ref};
my $blat=$options{blat_path};
my $samtools = $options{samtools_path};

my $fsa = $options{input}. ".out.fsa";
my $psl = $options{input}. ".out.psl";
open (TMP, ">", "$fsa") or die "Could not write temp sequence file\n$!\n";

open(IN, "<", $options{input}) or die "can not open $options{input} : $! \n";

## read in vcf file
while (<IN>) {
	chomp $_;

	## skip header/comments
	next if ($_ =~ /^#/);

	## for each variant create a fasta seq to be blat(ed)
	my @data = split(/\t/,$_);
	my $start = $data[1]-$options{window};
	my $end = $data[1]+$options{window};

	## create temp file fsa file to blat
	print TMP `$samtools/samtools faidx $options{reference} $data[0]:$start-$end`;
}

## close file handler
close(TMP);
close(IN);

## run blat
my $cmd = "$blat/blat $blat_ref $fsa -noHead -minScore=$options{minScore} -minIdentity=$options{minIdentity} $psl";
system("$cmd");

## create hash of blat resutls, exlude self-hit
my $blat_hash = &create_blat_hash();

## create new vcf file with ED value
new_vcf_file($blat_hash);

## remove tmp and psl files
system("rm $fsa");
system("rm $psl");

timer(); #call timer to see when process ended.
exit();

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output", "reference", "window");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}


	$options{minScore} = 70 unless($options{minScore});
	$options{minIdentity} = 90 unless($options{minIdentity});

	$options{blat_ref} = "/data2/bsi/reference/db/hg19.2bit" unless ($options{blat_ref});

	## be sure that blat executable is 64-bit and not 32-bit
	## 32-bit could cause memory error.
	$options{blat_path} = "/projects/bsi/bictools/scripts/dev/jbhavsar/apps" unless ($options{blat_path});
	$options{samtools_path} = "/projects/bsi/bictools/apps/alignment/samtools/samtools-0.1.18" unless ($options{samtools_path});
}

#############################################################################
sub timer {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    print "Time now: " . $theTime."\n";
}

#############################################################################
sub create_blat_hash {
	## open blat output and create a hash of query and count to be later added to vcf file.
	## following columns are expected in order in blat output.
	##  1 - match
	##	2 - mis-match
	##	3 - rep. match
	##	4 - N's
	##	5 - Q gap count
	##	6 - Q gap bases
	##	7 - T gap count
	##	8 - T gap bases
	##	9 - strand
	##	10 - Q name --- has to be of form chr##:start-end
	##	11 - Q size
	## 	12 - Q start
	##	13 - Q end
	##	14 - T name
	##	15 - T size
	##	16 - T start
	##	17 - T end
	##	18 - block count
	##	19 - blockSizes
	##	20 - qStarts
	##	21 - tStarts

	open (PSL, "<", $psl) or die "Could not open PSL file\n$!\n";
	my $hash = {};

	while(<PSL>) {
		chomp $_;

		my @data = split(/\t/, $_);
		my $key = $data[9];

		$data[9] =~ s/:|-/\t/g;

		my ($chr, $qStart, $qEnd) = split(/\t/, $data[9]);

		## create blat hash,  eliminating self hits
		if ((abs($qStart - $data[15]) > 1) || (abs($qEnd - $data[16]) > 1)){
			if (exists $hash->{$key}) {
				$hash->{$key} += 1;
			} else {
				$hash->{$key} = 1;
			}
		}
	}

	return $hash;
}

#############################################################################
sub new_vcf_file {
	my $hash_ref = shift;

	## open original vcf file and update comments field with ED=? from the hash above
	## vcf file format v4.1 with following columns expected in order.
	## 	1 - CHROM
	##	2 -	POS
	##	3 - ID
	##	4 - REF
	##	5 - ALT
	## 	6 - QUAL
	##	7 - FILTER
	##	8 - INFO
	## 	9 -	FORMAT
	## 	10 - s_933236

	open(IN, "<", $options{input}) or die "can not open $options{input} : $! \n";
	open(OUT, ">", $options{output}) or die "can not open $options{output} : $! \n";
	while (<IN>) {
		chomp $_;

		if ($_ =~ /^##/) {
			## print original comments to new file
			print OUT $_."\n";

		} elsif ($_ =~ /^#/) {
			## print ED comments before printing column header
			print OUT "##INFO=<ID=ED,Number=1,Type=Integer,Description=\"Number of blat hits to reference genome, not counting self-hit\">\n";
			print OUT $_ ."\n";

		} else {
			## add ED value to each variant call.
			my @data = split(/\t/, $_);
			my $start = $data[1]-$options{window};
			my $end = $data[1]+$options{window};
			my $key = $data[0].":".$start."-".$end;

			if (exists $hash_ref->{$key}) {
				$data[7] .= ";ED=".$hash_ref->{$key};
			} else {
				$data[7] .= ";ED=-1";
			}

			print OUT join("\t", @data)."\n";
		}
	}
	close(IN);
	close(OUT);
}
