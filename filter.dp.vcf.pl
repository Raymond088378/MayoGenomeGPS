use strict;
use warnings;
use Getopt::Long;

my($vcf, $output, $gzipped, $sample, $filter_depth,$help);

GetOptions("in|i=s"     => \$vcf,
"out|o:s"       => \$output,
"gzipped|z"     => \$gzipped,
"sample|s:s"    => \$sample,
"depth|d:i"     => \$filter_depth,
"help|h|?|"     => \&help);

if (not $vcf)   {
    print "Missing input file!\n";
	help();  
}
if($gzipped){
	open VCF, "gunzip -c $vcf |" or die "opening gzipped $vcf\n";
}else{
	open VCF, "<$vcf" or die "opening $vcf\n";
}

open OUT, ">$output" or die "opening $output\n" if defined $output;

my ($format,$pos,$dp);
$dp=0;
while(my $l = <VCF>)	{
    if ($l =~ /^##/)    {
            print OUT "$l";
    }
    elsif ($l =~ /^#CHROM/)	{
		chomp $l;
		my @call=split(/\t/,$l);
		for (my $i=0; $i <=$#call;$i++)	{
			if ($call[$i] eq 'FORMAT')	{
				$format=$i;
			}
			if ($call[$i] eq "$sample")	{
				$pos=$i;
			}
		}
        print OUT "$l\n";
    }
    else    {
       my @call=split(/\t/,$l);
       my @parse=split(/:/,$call[$format]);
       for (my $j=0; $j <=$#parse;$j++)	{
           if ($parse[$j] eq 'DP')	{
               $dp=$j;
               last;
           }
       }
        my @data=split(/:/,$call[$pos]);
		if ( $data[$dp] >= $filter_depth)	{
            print OUT "$l\n";
		}
    }
}
close VCF;

sub help{
    
	print "DESCRIPTION:
	filter.dp.vcf.pl filter the vcf file and outputs VCF format. 
    
USAGE:
	filter.dp.vcf.pl -i input.vcf.gz -o input.filter.vcf -z
    
OPTIONS:
	--in,-i		path to input file. Required parameter. Input can be gzipped as long 
                as the -z flag is also used. 
    
	--out,-o 	path to VCF output file. 
    
	--gzipped,-z 	A flag used if the input file is gzipped.
        
    --sample,-s 	sample name can be defined. 
        
    --depth,-d	Optional parameter to filter variants by read depth. A variant will be
                skipped if the total read depth at the position of the variant is less 
                than this value.

                
    --help,-h,-?	Display this documentation.
                
                \n";
            exit 1;
}            



