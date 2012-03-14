use strict ; 
use warnings;

our $file = $ARGV[0]; 
our $reads = $ARGV[1];
my $num_lines=$reads*4;
our $min = 255 ; 
our $max = 0 ; 

#   $SIG{'INT'} = sub { print "overall MIN $min / MAX $max\t$file\n"; exit ; };
open FH, "$file" or die "";

while(<FH>){
    if ($. < $num_lines){
	if($_ =~ /^\+/)	{
	    my $line = <FH>; # get the quality line
	    chomp $line;
	    #print "$line\n";
	    my $fred_scores = $line; 
	    my $converted_fred_scores ; 
	
	    my $linemin = 255 ; 
	    my $linemax = 0 ; 
	
	    foreach my $this_score (split(//,$fred_scores)){
		$converted_fred_scores .= chr(ord($this_score) - 31) ; 
	
		my $score = ord($this_score); 
		if($score > $linemax){ $linemax = $score ; }
		if($score < $linemin){ $linemin = $score ; }
		
		if($linemax != 0 && $linemax > $max){ $max = $linemax ; }
		if($linemin != 255 && $linemin < $min){ $min = $linemin ; }
	
	    }
	}
    }
    else	{
	last;
    }	
}

print "$min";


