

##chromosome      position        reference       Change  Homozygous      Bio_type        accession       Exon_ID Exon_Rank       functionGVS     aminoAcids      proteinPosition Codon_Degeneracy        geneList
##chr1    761958  *       +T      Hom     -       -       -       -       INTERGENIC      -       -       -       -
##chr1    978604  *       -CT     Het     mRNA    NM_198576       -       -       INTRON  -       -       -       AGRN

## Chr     Start   Stop    dbSNP132        dbSNP132Alleles COSMIC  InCaptureKit    #AlternateHits  Ref     Alt     Base-Length     Indel-supportedRead     ReadDepth
##chr1    761957  761957  rs70949521      -/T     -       0       1       A       AT      1       106     106
##chr1    978603  978605  -       -       -       1       0       CCT     C       2       40      81

use strict;
use Getopt::Std;

our ($opt_i, $opt_s, $opt_o);
print "INFO: script to add snpeff indel results to the variant report\n";
print "RAW paramters: @ARGV\n";
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
	my $allele;
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
			my @a = split(/\t/,$line);
			if (length($a[8]) eq 1)	{		##INS
				$allele=substr($a[9],1);
				${$hashreport{$a[1]}{$allele}{'+'}}=$line;
			#	print "$pos\t$allele\t+\n";
			#	<STDIN>;
			}
			if (length($a[8]) gt 1 )	{  ## DEL
				$allele=substr($a[8],1);
				${$hashreport{$a[1]}{$allele}{'-'}}=$line;
			#	print "$pos\t$allele\t-\n";
			#	<STDIN>;
			}	
		}
	}	
	close REPORT;

	#Form a unique key with chr_pos from snpeff, push the duplicates per key into array
	while(my $line = <SNPEFF>)	{
		chomp $line;
		my @a = split(/\t/,$line);
		my $ref_length=length($a[2]);
		my $alt_length=length($a[3]);
		my ($allele);
		if ($ref_length == 1)	{
			$allele = substr($a[3],1);
			push( @{$hasheff{$a[1]}{$allele}{'+'}},$line);
			#print "$a[1]\t$allele\t+\n";
			#<STDIN>;
		}
		else	{
			$allele = substr($a[2],1);
			push( @{$hasheff{$a[1]}{$allele}{'-'}},$line);
			#print "$a[1]\t$allele\t$-\n";
			#<STDIN>;
		}
		#print "$a[1]\t$allele\t$sign\n";
		#<STDIN>;
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
