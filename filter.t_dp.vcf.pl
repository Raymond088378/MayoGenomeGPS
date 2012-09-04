use strict;
use warnings;
use Getopt::Long;

my($vcf, $output, $gzipped, $sample, $filter_depth,$help);

GetOptions("in|i=s"     => \$vcf,
"out|o:s"       => \$output,
"gzipped|z"     => \$gzipped,
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

my ($info,$pos,$dp);
$dp=0;
while(<VCF>)	{
	if (/^\#\#/) {
		print OUT $_;
		next;
    }
	#Parse the header
	if (/^\#/) {
       chomp;
       my $line = $_;
	   $line =~ s/\#//;
       @nfields = split (/\t/, $line) ;
       @samples = @nfields[8..$#nfields];
       print OUT "#$line\n";
       next;
   }
   my %fields;
   my @values = split/\t/;

   #Parse the columns;
   foreach my $i (@nfields) {
	$fields{$i} = shift @values;
    }
   
   #Parse the values of the info field
   my @names_info = split (/;/, $fields{'INFO'});
   
   
   
}
close VCF;

sub help{
    
	print "DESCRIPTION:
	filter.t_dp.vcf.pl filter the vcf file and outputs VCF format. 
    
USAGE:
	filter.dp.vcf.pl -i input.vcf.gz -o input.filter.vcf -z
    
OPTIONS:
	--in,-i		path to input file. Required parameter. Input can be gzipped as long 
                as the -z flag is also used. 
    
	--out,-o 	path to VCF output file. 
    
	--gzipped,-z 	A flag used if the input file is gzipped.
        
    --depth,-d	Optional parameter to filter variants by read depth. A variant will be
                skipped if the total read depth at the position of the variant is less 
                than this value.

                
    --help,-h,-?	Display this documentation.
                
                \n";
            exit 1;
}            



