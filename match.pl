use strict;
use warnings;
use Getopt::Std;

our ($opt_b, $opt_i,$opt_o);
print "RAW paramters: @ARGV\n";
getopt('bio');
if ( (!defined $opt_b) && (!defined $opt_i)  && (!defined $opt_o)) {
    die ("Usage: $0 \n\t-b [bed file] \n\t-i [interesect file]\n\t-o [out file]\n");
}
else    {
    my $bed=$opt_b;
    my $intersect=$opt_i;
    my $output= $opt_o;
    
    open BED, "$bed" or die "can not open $bed : $! \n";
    open IN, "$intersect" or die "can not open $intersect : $! \n";
    open OUT, ">$output" or die "can not open $output : $!\n";
    
    my $bed_l=<BED>;
    my $in_l=<IN>;
    # disable the output buffer to ensure we don't use too much memory
    my $old_fh = select(OUT); $| = 1; select($old_fh);
    
    while(defined $bed_l || defined $in_l){
        chomp $bed_l if defined $bed_l;
        chomp $in_l if defined $in_l;
        if(!defined $in_l){
            print OUT "-\n";
		$bed_l=<BED>;	
        }
        else    {
            my @bed_array=split(/\t/,$bed_l);
            my @in_array=split(/\t/,$in_l);
            if ($bed_array[1] == $in_array[1])  {
                print OUT "$in_array[$#in_array]\n";
                $in_l=<IN>;
                $bed_l=<BED>;  
            }
            elsif ($bed_array[1] > $in_array[1])    {
                print OUT "-\n";
                $in_l=<IN>;
            }
            elsif ($bed_array[1] < $in_array[1]){
                print OUT "-\n";
                $bed_l=<BED>;
            }
        }
    }
    close BED;
    close OUT;
    close IN;
}    
