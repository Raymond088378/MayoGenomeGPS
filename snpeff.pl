use strict;
use warnings;

my $eff=shift @ARGV;
### Chromo        Position        Reference       Change  Change_type     Homozygous      Quality Coverage        Warnings        Gene_ID Gene_name       Bio_type        Trancript_ID    Exon_ID Exon_Rank       Effect  #old_AA/new_AA   Old_codon/New_codon     Codon_Num(CDS)  Codon_Degeneracy        CDS_size        Codons_around   AAs_around      Custom_interval_ID

open FH, "$eff" or die "can not open $eff : $!\n";
my $cols;
while(my $l = <FH>)	{
	chomp $l;	
	my @a = split(/\t/,$l)  ;
	next if (( $. ==2) || ($. == 1));
	if ($. == 3)    {
		print "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tfunctionGVS\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n";
		$cols=$#a;
	}
	else	{
		my $lens=$#a;
		for(my $i = 0; $i <= $lens;$i++)        {
                        if (length($a[$i]) == 0)        {
                                $a[$i]='-';
                        }
                }
                if ($lens < $cols)      {
                        for (my $j=$lens+1;$j <=$cols; $j++)    {
                                $a[$j]='-';
                        }
                }
		my @cols=(0,1,2,3,5,11,12,13,14,15,16,18,19,10);
		print join ("\t",@a[@cols]);
		print "\n";
	}
}
close FH;
