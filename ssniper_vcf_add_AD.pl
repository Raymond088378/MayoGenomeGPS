#!/usr/bin/perl -w

use strict;
my @nfields;
my @samples;

while (<>) {
   if (/^\#\#/) {
	print $_;
	next;
    }

   #Parse the header
   if (/^\#/) {
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
       my $ad = $freq{$sample}->{$fields{'ALT'}};
       $sample_values{$sample}->{AD}=$ad;
   }

   #Change fields to reflect change
   push (@names_format, "AD");
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
