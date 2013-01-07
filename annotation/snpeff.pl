#!/usr/local/biotools/perl/5.10.0/bin/perl

### Chromo        Position        Reference       Change  Change_type     Homozygous      Quality Coverage        Warnings        Gene_ID Gene_name       Bio_type        Trancript_ID    Exon_ID Exon_Rank       Effect  #old_AA/new_AA   Old_codon/New_codon     Codon_Num(CDS)  Codon_Degeneracy        CDS_size        Codons_around   AAs_around      Custom_interval_ID
##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )' ">
#NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|tGt/tTt|C195F|350|AADACL3||CODING|NM_001103170|NM_001103170.ex.4)
#chr1    12785494        G       T       -       coding  NM_001103170    NM_001103170.ex.4       4       NON_SYNONYMOUS_CODING   C/F     195     1       AADACL3

my $cols;
print "chromosome\tposition\treference\tChange\tHomozygous\tBio_type\taccession\tExon_ID\tExon_Rank\tEffect\taminoAcids\tproteinPosition\tCodon_Degeneracy\tgeneList\n";
while(<>)	{
	next if (/^#/);
	$l=$_;
	chomp $l;	
	my @a = split(/\t/,$l)  ;
	my @b = split(/;/,$a[7]);
	my @eff = split(/,/,$b[$#b]);
	for (my $i=0;$i<=$#eff;$i++)	{
		my @trans= split(/[( | )]/,$eff[$i]);	
		if (length($trans[9])>0)	{$accession=$trans[9]; }
		else	{$accession="-";}
		if (length($trans[8])>0)	{$bio_type=lc($trans[8]); }
		else	{$biotype="-";}
		if (length($trans[10])>0)	{$exon_id=$trans[10]; $exon_id=~ s/\)//g; }
		else	{$exon_id="-";}
		if (length($trans[0])>0)	{$effect=$trans[0]; $effect=~ s/EFF=//g; }
		else	{$effect="-";}
		if (length($trans[4])>0)	{
			$trans[4] =~ m/([A-Z]*)(\d+)([A-Z]*)/;
			if (length($1)>0 && length($3)>0)	{$aa=$1."/".$3;}
			elsif(length($1)>0)	{$aa=$1."/".$1;}
			else	{$aa="-/".$3;}
				$pp=$2;}
		else	{$aa="-";$pp="-";}
		if (length($trans[6])>0)	{$gene=$trans[6]; }
		else	{$gene="-";}	
			
		print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t-\t$bio_type\t$accession\t$exon_id\t-\t$effect\t$aa\t$pp\t-\t$gene\n";
	}
}
