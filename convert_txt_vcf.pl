use strict;
use warnings;

my $input = shift @ARGV;
my $sample = shift @ARGV;
open FH, "$input" or die " opening $input : $!\n";
## chr, pos, ref alt

header($sample);
while(my $l = <FH>) {
    chomp $l;
    my @a = split (/\t/,$l);
    print "$a[0]\t$a[1]\t.\t$a[2]\t$a[3]\t.\tPASS\tNS=1\tGT:AD:DP:GQ\t.:.,.:.:99\n";
}

close FH;

sub spGetCurDateTime {
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
    $year+1900, $mon+1, $mday, $hour, $min, $sec;
    return ($curDateTime);
}

sub header{
    my $sm = $_[0];
    my $date=&spGetCurDateTime();
    my $header = qq{##fileformat=VCFv4.1
##fileDate=$date
##source=convert_txt_vcf.pl
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth for This Sample">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sm}\n};
print $header;
}