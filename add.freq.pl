#!/usr/local/biotools/perl/5.10.0/bin/perl
## Saurabh Baheti
## April 13th 2012
## contact baheti.saurabh@mayo.edu

## scrip to add rsid column to the input file

use strict;
use warnings;
use Getopt::Std;

our ($opt_i, $opt_h, $opt_k , $pot_p, $opt_o);
getopt('ihkpo');
if ( (!defined $opt_i) && (!defined $opt_h) && (!defined $opt_k) && (!defined $opt_p) && (!defined $opt_o) ){
    die ("Usage: $0 \n\t -i [input file] \n\t -h [ hapmap source file] \n\t -k [ kgenome source file ] \n\t -p [ population ] -o [ output file]\n");
}
else{
    my $infile = $opt_i;    
    my $hapmap = $opt_h;
	my $kgenome = $opt_k;
	my $population = $opt_p;
    my $outfile = $opt_o;
    ### 1 base means subtract one from dbsnp file pos
    open IN, "<$infile" or die "can not open file $infile: $!\n";
    open OUT, ">$outfile" or die "can not open file $outfile: $!\n";
    open HAP, "<$hapmap" or die " can not open file $hapmap: $!";
	open KGEN, "<$kgenome" or die " can not open file $kgenome: $!";
    my $header .= <IN>;
	chomp $header;
	$header .= "\tHapMap_".$population."_allele_freq\t1kgenome_".$population."_allele_freq";
	print OUT "$header\n";
    
    # disable the output buffer to ensure we don't use too much memory
    my $old_fh = select(OUT); $| = 1; select($old_fh);
    my $input_c=<IN>;
    my $in_hap=<HAP>;
	my $in_kgen=<KGEN>;
	
    while(defined $in_hap || defined $input_c || defined $in_kgen){
        chomp $input_c if defined $input_c;
        chomp $in_hap if defined $in_hap;
		chomp $in_kgen if defined $in_kgen;
		if (!defined $input_c)	{last;}
        elsif (!defined $snp)  {
            print OUT "$input_c\t";
            print OUT "-\t-\n";
            $input_c=<IN>;
        }
        else{
			
			my @in_array = split(/\t/,$input_c);
			my @snp_array = split(/\t/,$snp);
			if($in_array[0] eq $snp_array[1] && $in_array[1] == $snp_array[3]) {
				print OUT "$input_c\t";
				print OUT "$snp_array[4]\t$snp_array[9]\n";
				$snp=<SNP>;
				$input_c=<IN>;
			}
			elsif($in_array[0] eq $snp_array[1] && $in_array[1] > $snp_array[3]) {
				$snp=<SNP>;
			}
			elsif ($in_array[0] eq $snp_array[1] && $in_array[1] < $snp_array[3])   {
				print OUT "$input_c\t";
				print OUT "-\t-\n";
				$input_c=<IN>;
			}
		}
    }
}
close IN;
close SNP;
close OUT;

## end of the script

    
