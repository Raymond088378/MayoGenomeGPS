
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
	
B<--blat_server, -bs>
	name of the server to run the blat on Default is localhost

B<--threads, -th>
	number of threads to use ot run this Default is 1 
	
B<--blat_port, -bp>
	full path to the samtools tool	
	
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
use Data::Dumper;
use Cwd;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use threads;
use threads::shared;
use POSIX;
 
my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
						  'output|o=s',
						  'samtools_path|sam=s',
						  'blat_path|b=s',
						  'blat_ref|br=s',
						  'blat_server|bs=s',
						  'blat_port|bp=s',
						  'reference|r=s',
						  'window|w=s',
						  'minMatch|m=s',
						  'minMisMatch|s=s',
						  'minIdentity|t=s',
						  'threads|th=s',
						  'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
#############################################################################
## make sure everything passed was peachy
&check_parameters(\%options); 
 
print "Starting main program\n";
my @threads;
my $input=$options{input};
my $output = $options{output};
my $blat_log = "$input.blat.log";
my $blat_ref=$options{blat_ref};
my $blat=$options{blat_path};
my $blat_server = $options{blat_server};
my $blat_port = $options{blat_port};
my $samtools = $options{samtools_path};
my $num_threads = $options{threads};

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


open OUT , ">$output" or die "can not open $output : $! \n";
open IN, "$input" or die "can not open $input : $! \n";
my $head=<IN>;
my $skip=0;	
while($head =~ m/^##/)	{
	print OUT "$head";
	$head=<IN>;
	$skip++;
};
print OUT "$head";
print OUT "##INFO=<ID=ED,Number=1,Type=Integer,Description=\"Number of blat hits to reference genome, not counting self-hit \">\n";
$skip++;
my $len=`cat $input | awk '\$0 !~ /^#/' | wc -l`;

for ( my $count = 1; $count <= $num_threads; $count++) {
	my $start=$skip + ceil(($count-1)*($len/$num_threads)+1);
	my $end=$skip + ceil(($len/$num_threads)*$count);
	if ($end > $len+$skip)	{
		$end = $len+$skip;
	}
	my $t = threads->create(\&blat, $count, $input, $start, $end, $options{window}, $samtools, $options{reference}, $blat, $blat_server, $blat_port, $options{minScore}, $options{minIdentity});
	push(@threads,$t);
}
foreach (@threads) {
	my $num = $_->join;
	print "done with $num\n";
}
for ( my $count = 1; $count <= $num_threads; $count++) {
	my $out="$input.$count.out";
	open OUT1, "$out" or die "";
	while(<OUT1>)	{
		print OUT $_;
	}
	close OUT1;
	`rm $out`;
}
close OUT;	

print "End of main program\n";
 
sub blat {
	my $num = shift; 
	print "started thread $num\n";
	my $file = shift;
	my $start = shift;
	my $end = shift;
	my $window = shift;
	my $samtools = shift; 
	my $ref = shift ;
	my $blat = shift;
	my $blat_server = shift;
	my $blat_port = shift;
	my $minScore = shift;
	my $minIdentity = shift;

	open FH, "$file" or die "";
	my $dest="$file.$num.out";
	open OUT, ">$dest" or die "";
	my $fsa="$file.$num.out.fsa";
	my $psl="$file.$num.out.psl";
	while (<FH>)	{
		chomp $_;
		next if ( ($. > $end) || ($. < $start));
		next if ($_ =~ m/^#/);
		my @data = split(/\t/,$_);
		#output extracted seq to temp file.
		my $start = $data[1]-$window;
		my $end = $data[1]+$window;

		#create temp file fsa file to blat
		open (TMP, ">", "$fsa") or die "Could not write temp sequence file\n$!\n";
		print TMP `$samtools/samtools faidx $ref $data[0]:$start-$end`;
		close(TMP);
		### execute the blat
		`$blat/gfClient $blat_server $blat_port -nohead -minScore=$minScore -minIdentity=$minIdentity / $fsa $psl`;
		my $b_count = `cat $psl | wc -l`;
		#do not count self hit, assumed self-hit is always in the output list.
		if ($b_count >= 1){
			$b_count -= 1;
		} else { $b_count = -1 };
		#append count to description.
		$data[7] .= ";ED=$b_count";
		print OUT join("\t",@data)."\n";
	}
	close OUT;
	`rm $fsa $psl`;
	print "done with thread $num\n";
	return $num;
}	


#############################################################################
sub check_parameters {
    my $options = shift;

	my @required = ("input", "output", "reference", "window", "samtools_path" , "blat_path", "blat_ref", "blat_port");

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
	unless($options{threads}){
		$options{threads} = 1;
	}

	unless($options{blat_server}){
		$options{blat_server} = "localhost";
	}
}


