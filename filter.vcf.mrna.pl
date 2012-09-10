#!/usr/bin/perl -w

use strict;
my $ad;

while (<>) {
   if (/^\#\#/) {
	print $_;
	next;
    }
    
   #Parse the header
   if (/^\#/) {
       print "##INFO=<ID=RATIO,Number=1,Type=Float,Description=\"ratio of alt/total reads\">\n";
       chomp;
       my $line = $_;
       print "$line\n";
       next;
   }
   
   my $line = $_;
   chomp $line;
   my @values = split(/\t/,$line);
   my @format=split(/:/,$values[8]);
   for(my $i=0; $i <=$#format;$i++) {
        $ad=$i if ($format[$i] =~ /AD/) ;
   }
   if (defined $ad)	{
	   my @sample=split(/:/,$values[9]);
	   my ($ref,$alt)=split(/,/,$sample[$ad]);
	   if ($alt >=2)    {
			my $depth=$ref+$alt;
			my $ratio=sprintf("%.2f",$alt/$depth);
			if ($ratio >= 0.1)  {
				print join ("\t",@values[0..7]) . ";RATIO=$ratio\t"; 
				print join ("\t",@values[8..9]) . "\n";
			}
		}	
	}
	undef $ad;
}    
             