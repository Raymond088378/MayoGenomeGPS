#!/usr/bin/perl

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
  Jaysheel D. Bhavsar @ bjaysheel[at]gmail[dot]com


==head1 EXAMPLE
   vcf_blat_verify.pl --input /file/path/filename.vcf --output /file/path/output/filename.vcf --reference /file/path/hg19.fsa --window 50

=cut

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'blat_path|b=s',
						  'blat_ref|br=s',
						  'blat_server|bs=s',
						  'blat_port|bp=s',
						  'reference|r=s',
						  'window|w=s',
						  'minMatch|m=s',
						  'minMisMatch|s=s',
						  'minIdentity|t=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options);

#my $dir = getcwd;
my $input=$options{input};
my $tmp_fsa = "$input.tmp_seq.fsa";
my $tmp_psl = "$input.tmp_seq.psl";

my $blat_log = "$input.blat.log";
my $blat_ref=$options{blat_ref};
my $blat=$options{blat_path};
my $blat_server = $options{blat_server};
my $blat_port = $options{blat_port};

#check if server is already running
#make sure blat server is up before proceeding
my $cmd = "$blat/gfServer status $blat_server $blat_port | wc -l";
my $status = `$cmd`;

if ($status == 0){
	print STDERR "INFO: Initializing gfServer\n";

	#start blat server improve blat search response.
	system("$blat/gfServer start $blat_server $blat_port -log=$blat_log $blat_ref &");

	if ($? != 0){
		print STDERR "ERROR: Count not init BLAT gfServer\n$!\n";
		exit(-1);
	}

	print STDERR "INFO: Checking if server is ready\n";

	my $sec = "30";

	while (1) {
		$status = `$cmd`;
		print $?."\n";

		#server is up exit while loop
		unless ($status == 0){
			last;
		}

		sleep $sec;
		$sec += $sec;
	}

	print STDERR "INFO: Server ready\n";
}

#open input file
open (FHD, "<", $options{input}) or die "Could not open file $options{input}\n$!\n";
open (OUT, ">", $options{output}) or die "Could not create output file $options{output}\n$!\n";

# read input vcf file and create a temp multi fasta file to BLAT
while (<FHD>){
	#output comments and skip

	if ($_ =~ /^#/){
		#add info for blat homologs.
		if ($_ =~ /^#CHROM/){
			print OUT "##INFO=<ID=ED,Number=1,Type=Integer,Description=\"Number of blat hits to reference genome, not counting self-hit \">\n";
		}
		print OUT $_;
		next;
	}

	chomp $_;

	#split vcf data
	my @data = split(/\t/,$_);

	#output extracted seq to temp file.
	my $start = $data[1]-$options{window};
	my $end = $data[1]+$options{window};

	#create temp file fsa file to blat
	open (TMP, ">", "$tmp_fsa") or die "Could not write temp sequence file\n$!\n";
	print TMP `samtools faidx $options{reference} $data[0]:$start-$end`;
	close(TMP);

	#execute blat search.
	`$blat/gfClient $blat_server $blat_port -nohead -minScore=$options{minScore} -minIdentity=$options{minIdentity} / $tmp_fsa $tmp_psl`;

	if ($? == -1){
		print STDERR "ERROR: Problem executing gfClient\n$!\n";
		exit(-1);
	}

	#BLAT result count.
	my $b_count = `cat $tmp_psl | wc -l`;

	#do not count self hit, assumed self-hit is always in the output list.
	if ($b_count >= 1){
		$b_count -= 1;
	} else { $b_count = -1 };


	#append count to description.
	$data[7] .= ";ED=$b_count";

	print OUT join("\t",@data)."\n";
}

close(FHD);
close(OUT);

print STDERR "INFO: Clearing tmp files and shutting down gfServer\n";
`rm $tmp_psl`;
`rm $tmp_fsa`;
`rm $blat_log`;

#gfServer stop doesn't work.
#`gfServer stop $blat_server $blat_port`;
#$status = `ps | grep "gfServer" | cut -d ' ' -f1`;
#chomp $status;

#print "kill $status\n";

#`kill $status`;



exit(0);

#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output", "reference", "window", "blat_path", "blat_ref", "blat_server", "blat_port");

	foreach my $key (@required) {
		unless ($options{$key}) {
			print STDERR "ARG: $key is required\n";
			pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});
			exit(-1);
		}
	}

	unless($options{minScore}){
		$options{minScore} = 70;
	}

	unless($options{minIdentity}){
		$options{minIdentity} = 90;
	}
}
