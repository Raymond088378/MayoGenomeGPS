use strict;
use warnings;
use Getopt::Std;

## considering assumption as 2nd column in the file is position and works for per chromosome file

our($opt_s,$opt_i,$opt_n,$opt_p,$opt_o);
print "RAW paramters: @ARGV\n";
getopt('sinpo');

if( (!defined $opt_s) && (!defined $opt_i) && (!defined $opt_n) && (!defined $opt_p) && (!defined $opt_o) ){
    die ("Usage: $0 \n\t -s [snv file] \n\t -i [indel file] \n\t -n [distance] \n\t -p [col for position] \n\t -o [output snv file]\n" );
}
else{
    my $snv=$opt_s;
    my $indel=$opt_i;
    my $distance=$opt_n;
    my $pos=$opt_p;$pos=$pos-1;
    my $output=$opt_o;
    open SNV , "$snv" or die "can not open $snv : $!\n";
    open INDEL, "$indel" or die "can not open $indel : $!\n";
    open OUT ,">$output" or die "can not open $output : $!\n";
    my $snv_l=<SNV>;
    my $indel_l=<INDEL>;
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
				my $start=$snv_array[$pos]-$distance;
				my $stop=$snv_array[$pos]+$distance;
				if( ($indel_array[$pos] > $start) && ($indel_array[$pos] < $stop) ){
					print OUT "$snv_l\t1\n";
					$snv_l=<SNV>;
					$indel_l=<INDEL>;
				}
				elsif($indel_array[$pos] < $snv_array[$pos]){
					print OUT "$snv_l\t0\n";
					$snv_l=<SNV>;
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
