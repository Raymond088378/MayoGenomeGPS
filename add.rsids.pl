#!/usr/local/biotools/perl/5.10.0/bin/perl
## Saurabh Baheti
## April 13th 2012
## contact baheti.saurabh@mayo.edu

## scrip to add rsid column to the input file

use strict;
use warnings;
use Getopt::Std;

our ($opt_i, $opt_s, $opt_o);
getopt('iso');
if ( (!defined $opt_i) && (!defined $opt_s) && (!defined $opt_o) ){
    die ("Usage: $0 \n\t -i [input file] \n\t -s [ dbsnp file] \n\t -o [ output file]\n");
}
else{
    my $infile = $opt_i;    
    my $dbsnp_file = $opt_s;
    $dbsnp_file =~ m/.+dbSNP(\d+)/;
    my $dbsnp_v = $1;
    my $outfile = $opt_o;
    ### 1 base means subtract one from dbsnp file pos
    open IN, "<$infile" or die "can not open file $infile: $!\n";
    open OUT, ">$outfile" or die "can not open file $outfile: $!\n";
    open SNP, "<$dbsnp_file" or die " can not open file $dbsnp_file: $!";
    my $header .= <IN>;
    chomp $header;
    $header .= "\tdbSNP".$dbsnp_v."\tdbSNP".$dbsnp_v."Alleles";
    print OUT "$header\n";
    
    # disable the output buffer to ensure we don't use too much memory
    ### input chr 0 pos 1
    ### dbsnp chr 1 pos 2 allele 9 stop 3
    my $old_fh = select(OUT); $| = 1; select($old_fh);
    my $input_c=<IN>;
    my $snp=<SNP>;
    while(defined $snp || defined $input_c){
        chomp $input_c if defined $input_c;
        chomp $snp if defined $snp;
        if (!defined $input_c)	{
            last;
		}
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

    
