#!/usr/local/biotools/perl/5.10.0/bin/perl -w

use strict;
my @nfields;
my @samples;

while (<>) {
   if (/^\#\#/) {
	if ( /^##FORMAT=<ID=BCOUNT/)	{
		print "##FORMAT=<ID=BCOUNT,Number=.,Type=Integer,Description=\"Occurrence count for each base at this site (A,C,G,T)\">\n";	
	}
	else if(/^##FORMAT=<ID=DP4/){
		print "##FORMAT=<ID=DP4,Number=.,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n";
	}else	{
		print $_;
		next;
    }
}
   #Parse the header
   if (/^\#/) {
       print "##FORMAT=<ID=PGERM,Number=1,Type=Float,Description=\"probability of germ line call\">\n";
   print "##FORMAT=<ID=PLOH,Number=1,Type=Float,Description=\"probability of Loss of heterozygosity\">\n";
   print "##FORMAT=<ID=PHETMUT,Number=1,Type=Float,Description=\"probability of Hetero zygous mutation\">\n";
   print "##FORMAT=<ID=PHOMMUT,Number=1,Type=Float,Description=\"probability of Homo zygous mutation\">\n";
   print "##FORMAT=<ID=PSOM,Number=1,Type=Float,Description=\"probability of somatic mutation\">\n";
   print "##FORMAT=<ID=PPS,Number=1,Type=Float,Description=\"post processed probability of somatic call\">\n";
   print "##FORMAT=<ID=POW,Number=1,Type=Float,Description=\"given the tumor sequencing depth, what is the power to detect a mutation at 0.3 allelic fraction * given the normal sequencing depth, what power did we have to detect (and reject) this as a germline variant\">\n";
   print "##FORMAT=<ID=IMPAIR,Number=1,Type=Integer,Description=\"number of reads which have abnormal pairing (orientation and distance)\">\n";
   print "##FORMAT=<ID=MQ0,Number=1,Type=Integer,Description=\"total number of mapping quality zero reads in the tumor and normal at this locus\">\n";
   print "##FORMAT=<ID=MUTX_LOD,Number=1,Type=Float,Description=\"log likelihood of ( normal being reference / normal being altered )\">\n";
   print "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
   print "##FORMAT=<ID=SOMATIC,Number=1,Type=Integer,Description=\"keep/Reject (confident call)\">\n";
   print "##FORMAT=<ID=INSC,Number=1,Type=Integer,Description=\"count of insertion events at this locus in tumor\">\n";
	print "##FORMAT=<ID=DELC,Number=1,Type=Integer,Description=\"count of deletion events at this locus in tumor\">\n";
chomp;
       my $line = $_;
       $line =~ s/\#//;
       @nfields = split (/\t/, $line) ;
       @samples = @nfields[9..$#nfields];
       print "#$line\n";
       next;
   }

   chomp;
   my %fields;
   my @values = split/\t/;

   #Parse the columns;
   foreach my $i (@nfields) {
	$fields{$i} = shift @values;
    }
   
   #Parse the values of the samples
   my @names_format = split (/:/, $fields{'FORMAT'});
   my (%sample_values, %freq);
   my @bases=qw /A C G T/;
   
   
   foreach my $sample (@samples) {
       my @s_values = split (/:/, $fields{$sample});
       foreach my $i (@names_format) {
	   my $cvalue = shift (@s_values);
	   $sample_values{$sample}->{$i}=$cvalue;
       }
       my @bfreqs = split (/,/,$sample_values{$sample}->{BCOUNT});

       foreach my $base (@bases) {
	   my $cfreq = shift (@bfreqs);
	   $freq{$sample}->{$base}=$cfreq;
       }
       #Calculate alternate allele
       my $ad;
	   if (grep (/,/,$fields{'ALT'}))	{
			my @alter=split(/,/,$fields{'ALT'});
			foreach my $a (@alter)	{
			$ad .= "," . $freq{$sample}->{$a} ;
			}
			
		}
		else	{
			$ad = "," . $freq{$sample}->{$fields{'ALT'}};
		}		
	    
       my $ad1 = $freq{$sample}->{$fields{'REF'}};
      my @add = ("PGERM","PLOH","PHETMUT","PHOMMUT","PSOM","PPS","INSC","DELC","POW","IMPAIR","MQ0","MUTX_LOD","SOMATIC");
      for (my $i=0;$i<=$#add;$i++)  {
          $sample_values{$sample}->{$add[$i]}=".";   }
	$sample_values{$sample}->{AD}="${ad1}${ad}";
   }   
   #Change fields to reflect change
   my @add = ("PGERM","PLOH","PHETMUT","PHOMMUT","PSOM","PPS","INSC","DELC","POW","IMPAIR","MQ0","MUTX_LOD","SOMATIC","AD");
   for (my $i=0;$i<=$#add;$i++)  {
      push (@names_format,$add[$i]);
      
   }
   $fields{'FORMAT'} = join (':', @names_format);
   
   foreach my $sample (@samples) {
       my @values;
       foreach my $name (@names_format) {
	   push @values, $sample_values{$sample}->{$name};
       }
       $fields{$sample} = join (':', @values); 
   }

   #Print fields
   my @line_fields;
   foreach my $field (@nfields) {
       push (@line_fields, $fields{$field});
   }
   print join ("\t",@line_fields) . "\n";
}
