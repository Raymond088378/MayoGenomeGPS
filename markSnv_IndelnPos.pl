#!/usr/local/biotools/perl/5.10.0/bin/perl
use Getopt::Std;

## considering assumption as 2nd column in the file is position and works for per chromosome file
## contact : baheti.saurabh@mayo.edu
### added the condition to solve a problem with deletions as a delteion can be longer than 10 bases and then the flag will show up as 0

our($opt_s,$opt_i,$opt_n,$opt_o);
print "RAW paramters: @ARGV\n";
getopt('sino');

if( (!defined $opt_s) && (!defined $opt_i) && (!defined $opt_n) && (!defined $opt_o) ){
    die ("Usage: $0 \n\t -s [snv file] \n\t -i [indel file] \n\t -n [distance] \n\t -o [output snv file]\n" );
}
else{
    my $snv=$opt_s;
    my $indel=$opt_i;
    my $distance=$opt_n;
    my $pos=1;		## 2nd column is pos in vcf files
    my $ref=3;
	my $alt=4;
	my $output=$opt_o;
    open SNV , "$snv" or die "can not open $snv : $!\n";
    open INDEL, "$indel" or die "can not open $indel : $!\n";
    open OUT ,">$output" or die "can not open $output : $!\n";
    my $snv_l=<SNV>;
    while($snv_l =~ /^#/)	{
	print OUT "$snv_l";
	$snv_l=<SNV>;
    }
    
    my $indel_l=<INDEL>;
    while($indel_l =~ /^#/)	{
	$indel_l=<INDEL>;
    } 
    
    # disable the output buffer to ensure we don't use too much memory
    my $old_fh = select(OUT); $| = 1; select($old_fh);

    while(defined $snv_l || defined $indel_l){
		chomp $snv_l if defined $snv_l;
		chomp $indel_l if defined $indel_l;
		if (!defined $snv_l){
			last;
		}
		else{	
			if(!defined $indel_l){
			    print OUT "$snv_l\t0\n";
			    $snv_l=<SNV>;
			}
			else{
			    my @snv_array = split(/\t/,$snv_l);
			    my @indel_array = split(/\t/,$indel_l);
				my $ref_allele=$indel_array[$ref];
				my $alt_allele=$indel_array[$alt];
				if (length($ref_allele) > 1)	{	
					$start=$snv_array[$pos]-$distance-length($ref_allele);
					$stop=$snv_array[$pos]+$distance;
				}
				else	{
					$start=$snv_array[$pos]-$distance;
					$stop=$snv_array[$pos]+$distance;
				}	
				
				if( ($indel_array[$pos] > $start) && ($indel_array[$pos] < $stop) ){
					print OUT "$snv_l\t1\n";
					$snv_l=<SNV>;
					$indel_l=<INDEL>;
			    }
			    elsif($indel_array[$pos] < $snv_array[$pos]){
					$indel_l=<INDEL>;
			    } 
			    else{
					print OUT "$snv_l\t0\n";
					$snv_l=<SNV>;
			    }
			}
		}
	}
	close SNV;
	close INDEL;
	close OUT;
}
