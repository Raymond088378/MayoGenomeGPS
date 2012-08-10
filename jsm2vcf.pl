use strict;
use warnings;
use Getopt::Long;

my($input, $output, $gzipped, $nsample, $tsample, $filter_depth, $filter_prob, $help);
GetOptions("in|i=s"	=> \$input,
	"out|o:s"	=> \$output,
	"gzipped|z"	=> \$gzipped,
	"normalsample|ns:s"	=> \$nsample,
    "tumorsample|ts:s"	=> \$tsample,
	"depth|d:i"	=> \$filter_depth,
	"prob|p:f"	=> \$filter_prob,
	"help|h|?|"	=> \&help);

if(not $input){
	print "Missing input file!\n";
	help();
	exit 1;
}
if($gzipped){
	open IN, "gunzip -c $input |" or die "opening gzipped $input\n";
}else{
	open IN, "<$input" or die "opening $input\n";
}
open OUT, ">$output" or die "opening $output\n" if defined $output;

# print VCF header
header($nsample,$tsample);
##chrom   position        ref_base        var_base        normal_counts_a normal_counts_b tumour_counts_a tumour_counts_b p_AA_AA p_AA_AB p_AA_BB p_AB_AA p_AB_AB p_AB_BB p_BB_AA p_BB_AB p_BB_BB post_processed_p_somatic	
while(my $row = <IN>){
    next if ($. == 1);
	chomp $row;
	my @line = split(/\t/,$row);
	my $chr = $line[0];
	my $pos = $line[1];
	my $ref = $line[2];
	my $alt = $line[3];
	my $normal_depth_ref = $line[4];
	my $normal_depth_alt = $line[5];
	my $normal_depth = $normal_depth_ref + $normal_depth_alt;
	my $tumor_depth_ref = $line[6];
	my $tumor_depth_alt = $line[7];
	my $tumor_depth = $tumor_depth_ref + $tumor_depth_alt;
	my $total_dp=$normal_depth + $tumor_depth;
	my $PGERM = $line[8] + $line[12] + $line[16];
	my $PLOH = $line[11] + $line[13];   
	my $PHETMUT = $line[9] + $line[16];
	my $PHOMMUT = $line[10] + $line[14];
	my $PSOM = $line[9] + $line[15] + $line[10] + $line[14];
	my $prob=0;
	my $id=0;
	my @ngeno=('0/0','0/0','0/0','0/1','0/1','0/1','1/1','1/1','1/1');
	my @tgeno=('0/0','0/1','1/1','0/0','0/1','1/1','0/0','0/1','1/1');
	my $c=0;
	for(my $i =8 ; $i <=14; $i++,$c++)	{
		if ($prob < $line[$i])	{
			$id=$c;
			$prob=$line[$i];
		}
	}
	print OUT "$chr\t$pos\t.\t$ref\t$alt\t.\tPASS\t" . "NS=2:DP=$total_dp:PGERM=$PGERM:PLOH=$PLOH:PHETMUT=$PHETMUT:PHOMMUT=$PHOMMUT:PSOM=$PSOM:PPS=$line[$#line]" . "\tGT:AD:DP:GQ\t" . "$ngeno[$id]:$normal_depth_ref,$normal_depth_alt:$normal_depth:99\t" . "$tgeno[$id]:$tumor_depth_ref,$tumor_depth_alt:$tumor_depth:99\n";	    
}
close IN;
close OUT;
sub header{
	my $nsm = $_[0];
    my $tsm = $_[1];
	my $date=&spGetCurDateTime();
    my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=jsm2vcf.pl
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth for This Sample">
##FORMAT=<ID=PGERM,Number=1,Type=Float,Description="probability of germ line call">
##FORMAT=<ID=PLOH,Number=1,Type=Float,Description="probability of Loss of heterozygosity">
##FORMAT=<ID=PHETMUT,Number=1,Type=Float,Description="probability of Hetero zygous mutation">
##FORMAT=<ID=PHOMMUT,Number=1,Type=Float,Description="probability of Homo zygous mutation">
##FORMAT=<ID=PSOM,Number=1,Type=Float,Description="probability of somatic mutation">
##FORMAT=<ID=PPS,Number=1,Type=Float,Description="post processed probability of somatic call ">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${nsm}\t${tsm}\n};
print OUT $header;
}

sub spGetCurDateTime {
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
    $year+1900, $mon+1, $mday, $hour, $min, $sec;
    return ($curDateTime);
}

sub help{

	print "DESCRIPTION:
	jsm2vcf.pl converts raw joint snvmix output to VCF format.

USAGE:
	jsm2vcf.pl -i input.gz -o sample.vcf -z

OPTIONS:
	--in,-i		Path to snvmix file. Required parameter. Input can be gzipped as long 
			as the -z flag is also used. Input files should not have a header line.

	--out,-o 	path to the output VCF output file.

	--gzipped,-z 	A flag used if the input file is gzipped.

	--normalsample,-ns 	normal sample name is defined.
        
	--tumorsample,-ts     tumor sample name is defined.           

	--help,-h,-?	Display this documentation.

";
}

