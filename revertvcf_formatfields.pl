#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
#use warnings;
use Getopt::Long;

my($original, $input, $vcf, $help);
GetOptions("org|o=s"	=> \$original,
	"in|i:s"	=> \$input,
	"vcf|v=s"	=> \$vcf,
	"help|h|?|"	=> \&help);

if(not $input){
	print "Missing input file!\n";
	help();
	exit 1;
}
if(not $original){
	print "Missing original file!\n";
	help();
	exit 1;
}
if(not $vcf){
	print "Missing output vcf file!\n";
	help();
	exit 1;
}

my @nfields;
my @nfields_in;
my @samples;
my @samples_pos;


open OUT, ">$vcf" or die "can not open $vcf : $!\n";
open ORG, "$original" or die "can not open $original : $!\n";
open IN, "$input" or die "can not open $input : $!\n";
while (my $l = <ORG>)	{
	if ($l =~ /^\#\#/) {
		print OUT "$l";
		next;
	}
	if ($l =~ /^\#/) {
       $l =~ s/\#//;
       chomp $l;
	   @nfields = split (/\t/, $l) ;
	@samples = @nfields[9..$#nfields];
       print OUT "#$l\n";
       next;
	}
}	

while (my $l = <IN>)	{
	if ($l =~ /^\#\#/) {
		next;
	}
	if ($l =~ /^\#/) {
       $l =~ s/\#//;
	   chomp $l;
       @nfields_in = split (/\t/, $l) ;
		next;
	}
}	

close ORG;
close IN;
open IN, "$input" or die "can not open $input : $!\n";
open ORG, "$original" or die "can not open $original : $!\n";

my $calls=<ORG>;
while ($calls =~ /^\#/)	{
	$calls=<ORG>;
}	
my $computed_calls=<IN>;
while ($computed_calls =~ /^\#/)	{
	$computed_calls=<IN>;
}	
my %fields=();
my %fields_in=();
my %sample_value=();
my %sample_value_in=();
while(defined $calls || defined $computed_calls){
	chomp $calls if defined $calls;
	chomp $computed_calls if defined $computed_calls;
	
	my @data_calls=split(/\t/,$calls);
	foreach my $i (@nfields) {
		$fields{$i} = shift @data_calls;
	}
	my @names_format = split (/:/, $fields{'FORMAT'});
	
	my @data_computed_calls=split(/\t/,$computed_calls);
	foreach my $i (@nfields_in) {
		$fields_in{$i} = shift @data_computed_calls;
    }		
	my @names_format_in = split (/:/, $fields_in{'FORMAT'});
	
	foreach my $sample (@samples) {
		if ($fields{$sample} =~ /\.\/\./)	{
			my $val=".:" x $#names_format . ".";
			$fields{$sample}=$val;
		}		
		my @s_values = split (/:/, $fields{$sample});
		if ($fields_in{$sample} =~ /\.\/\./)	{
			my $val=".:" x $#names_format_in . ".";
			$fields_in{$sample}=$val;
		}	
		my @s_values_in = split (/:/, $fields_in{$sample});
		my @req=('AD','DP','DP4');
		foreach my $i (@names_format)	{
			my $cvalue = shift (@s_values);
			$sample_value{$sample}->{$i}=$cvalue;
		}
		foreach my $i (@names_format_in)	{
			my $cvalue = shift (@s_values_in);
			$sample_value_in{$sample}->{$i}=$cvalue;
		}
		foreach my $i (@req)	{
			$sample_value{$sample}->{$i}=$sample_value_in{$sample}->{$i};
		}	
	}
	$fields{'FORMAT'} = join (':', @names_format);
	foreach my $sample (@samples) {
       my @values;
       foreach my $name (@names_format) {
	   push @values, $sample_value{$sample}->{$name};
       }
       $fields{$sample} = join (':', @values); 
   }
	 my @line_fields;
   foreach my $field (@nfields) {
       push (@line_fields, $fields{$field});
   }
   print OUT join ("\t",@line_fields) . "\n";
	$calls=<ORG>;$computed_calls=<IN>;
}	
	
	



