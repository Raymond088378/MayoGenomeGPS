use strict;
use warnings;
use Getopt::Long;

my($input, $output, $gzipped, $nsample, $tsample, $help);
GetOptions("in|i=s"	=> \$input,
	"out|o:s"	=> \$output,
	"gzipped|z"	=> \$gzipped,
	"normalsample|ns:s"	=> \$nsample,
    "tumorsample|ts:s"	=> \$tsample,
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
##contig	position	ref_allele	alt_allele	tumor_name	normal_name	score	dbsnp_site	covered	power	tumor_power	normal_power	total_pairs	improper_pairs	map_Q0_reads	t_lod_fstar	tumor_f	contaminant_fraction	contaminant_lod	t_ref_count	t_alt_count	t_ref_sum	t_alt_sum	t_ins_count	t_del_count	normal_best_gt	init_n_lod	n_ref_count	n_alt_count	n_ref_sum	n_alt_sum	judgement
	
while(my $row = <IN>){
    next if ($. == 1 || $. == 2);
	chomp $row;
	my @line = split(/\t/,$row);
	next if ($line[8] =~ m/^UNCOVERED/);
	my $chr = $line[0];
	my $pos = $line[1];
	my $ref = $line[2];
	my $alt = $line[3];
	my $power = $line[9];
	my $improper = $line[13];
	my $mq0 = $line[14];
	my $normal_depth_ref = $line[27];
	my $normal_depth_alt = $line[28];
	my $normal_depth = $normal_depth_ref + $normal_depth_alt;
	my $tumor_depth_ref = $line[19];
	my $tumor_depth_alt = $line[20];
	my $tumor_depth = $tumor_depth_ref + $tumor_depth_alt;
	my $total_dp=$normal_depth + $tumor_depth;
	my $insc = $line[23];
	my $delc = $line[24];
	my $mutlod = $line[26];
	my $somatic=0;
	if ( $line[$#line] eq 'KEEP')	{
		$somatic=1;
	}		
	my $n_gt = $line[25];
	my $n_gt_put;
	if ($n_gt eq "$ref$alt" || $n_gt eq "$alt$ref")	{
		$n_gt_put="0/1";
	}
	elsif ($n_gt eq "$alt$alt")	{
		$n_gt_put="1/1";
	}
	elsif ($n_gt eq "$ref$ref")	{
		$n_gt_put="0/0";
	}
	else	{
		$n_gt_put=".";
	}	
		
	if ($somatic == 1 ){
	print OUT "$chr\t$pos\t.\t$ref\t$alt\t.\tPASS\t" . "NS=2;DP=$total_dp;POW=$power;IMPAIR=$improper;MQ0=$mq0;MUTX_LOD=$mutlod;SOMATIC=$somatic" . "\tGT:AD:DP:GQ:INSC:DELC\t" . "$n_gt_put:$normal_depth_ref,$normal_depth_alt:$normal_depth:99:0:0\t" . "0/1:$tumor_depth_ref,$tumor_depth_alt:$tumor_depth:99:$insc:$delc\n";	    
}}
close IN;
close OUT;
sub header{
	my $nsm = $_[0];
    my $tsm = $_[1];
	my $date=&spGetCurDateTime();
    my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=mutect2vcf.pl
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth for This Sample">
##FORMAT=<ID=POW,Number=1,Type=Float,Description="given the tumor sequencing depth, what is the power to detect a mutation at 0.3 allelic fraction * given the normal sequencing depth, what power did we have to detect (and reject) this as a germline variant">
##FORMAT=<ID=IMPAIR,Number=1,Type=Integer,Description="number of reads which have abnormal pairing (orientation and distance)">
##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description="total number of mapping quality zero reads in the tumor and normal at this locus">
##FORMAT=<ID=MUTX_LOD,Number=1,Type=Float,Description="log likelihood of ( normal being reference / normal being altered )">
##FORMAT=<ID=SOMATIC,Number=1,Type=Integer,Description="keep/Reject (confident call)">
##FORMAT=<ID=INSC,Number=1,Type=Integer,Description="count of insertion events at this locus in tumor">
##FORMAT=<ID=DELC,Number=1,Type=Integer,Description="count of deletion events at this locus in tumor">
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
	mutect2vcf.pl converts raw mutect output to VCF format.

USAGE:
	mutect2vcf.pl -i input.gz -o sample.vcf -z

OPTIONS:
	--in,-i		Path to Mutect file. Required parameter. Input can be gzipped as long 
			as the -z flag is also used. Input files should not have a header line.

	--out,-o 	path to the output VCF output file.

	--gzipped,-z 	A flag used if the input file is gzipped.

	--normalsample,-ns 	normal sample name is defined.
        
	--tumorsample,-ts     tumor sample name is defined.           

	--help,-h,-?	Display this documentation.

";
}

