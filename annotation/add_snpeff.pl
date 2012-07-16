

use strict;
use Getopt::Std;

our ($opt_i, $opt_s, $opt_o);
#print "INFO: script to add snpeff results to the variant report\n";
#print "RAW paramters: @ARGV\n";
getopt('iso');
if ( (!defined $opt_i) && (!defined $opt_s) && (!defined $opt_o) ) {
	die ("Usage: $0 \n\t-i [variant file] \n\t-s [snpeff] \n\t-o [output file] \n");
}
else    {
	my $source = $opt_i;
	my $eff = $opt_s;
	my $dest = $opt_o;
	my %hashreport=();
	my %hasheff=();
	
	open REPORT, "<$source" or die " can not open $source : $! \n";
	open OUT, ">$dest" or die "can not open $dest :$! \n";
	open SNPEFF, "<$eff" or die " can not open $eff :$! \n";
	my $len_header=0;
	my $eff_head=<SNPEFF>;chomp $eff_head;
	my @eff_head_array=split(/\t/,$eff_head);
	my $len_eff=$#eff_head_array;
	my $num_bock=$len_eff+1;
	my $num_tabs=$#eff_head_array-4;
	while(my $line = <REPORT>)	{
		chomp $line;
		if ($. == 1)	{
			print OUT "$line". "\t"x 5 . "SNPEFF Annotation" . "\t" x $num_tabs . "\n";
		}	
		elsif($. == 2)	{
			chomp $line;
            print OUT "$line\t$eff_head\n";
		}	
		else	{
			chomp $line;
			my @array = split(/\t/,$line);
			${$hashreport{$array[1]}{$array[17]}{$array[18]}}=join("\t",@array);
		}
	}	
	close REPORT;

	#Form a unique key with chr_pos from snpeff, push the duplicates per key into array
	while(my $line = <SNPEFF>)	{
		chomp $line;
		my @array = split(/\t/,$line);
		push( @{$hasheff{$array[1]}{$array[2]}{$array[3]}},join("\t",@array));
	}
	close SNPEFF;

	#Loop over unique key from %hashreport and compare with %hashsseq;
	
	foreach my $pos (sort { $a <=> $b } keys %hashreport)	{
            foreach my $ref (sort keys %{$hashreport{$pos}})    {
                foreach my $alt (sort keys %{$hashreport{$pos}{$ref}})  {
                    if(defined $hasheff{$pos}{$ref}{$alt} )	{
                        for (my $i=0; $i <= $#{$hasheff{$pos}{$ref}{$alt}}; $i++)   {
                            print OUT "${$hashreport{$pos}{$ref}{$alt}}" . "\t";
                            print OUT "${$hasheff{$pos}{$ref}{$alt}}[$i]" . "\n";
                        }
                    }
                    else    {
                        print OUT "${$hashreport{$pos}{$ref}{$alt}}" . "\t-" x $num_bock . "\n"; 
                    }
                }
            }
		}    
	undef %hashreport;
	undef %hasheff;
}
	
close OUT;
