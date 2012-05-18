use strict;
#use warnings;

my $file = shift @ARGV;
my $run_info = shift @ARGV;
my $type = shift @ARGV;

my @line=split(/=/,`perl -ne "/^TOOL_INFO/ && print" $run_info`);
my $tool_info=$line[$#line];chomp $tool_info;
@line=split(/=/,`perl -ne "/^HTTP_SERVER/ && print" $tool_info`);
my $server=$line[$#line];chomp $server;
@line=split(/=/,`perl -ne "/^GENOMEBUILD/ && print" $run_info`);
my $genome=$line[$#line];chomp $genome;

## IGV = =Hyperlink("http://localhost:60151/goto?locus=chr1:761732","chr1:761732")
## TIssue = =Hyperlink("http://bmidev2/reference/Tissue-Specific.Gene.Expression/html/148398.html","link")
 ## pathway = =Hyperlink("http://bmidev2/reference/metacore_kegg_pathway_files/9636_pathways.xls","link")
 ## gene card = =Hyperlink("http://www.genecards.org/cgi-bin/carddisp.pl?gene=SAMD11","link")
## kaviar = =Hyperlink("http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl?chr=chr1&pos=761732&frz=hg19&onebased=1","link")
open FH , "$file" or die " cannot open $file : $!\n";
my $header=<FH>;chomp $header;
print "\t$header" . "\t" x 4 . "\n";
$header=<FH>;chomp $header;
my ($chr,$pos,$gene,$gene_id);
my @head=split(/\t/,$header);
for(my $i = 0; $i <=$#head;$i++)    {
    if ( $head[$i] eq 'geneList')   {
        $gene=$i;
    }
    if ($head[$i] eq 'Entrez_id')   {
        $gene_id=$i;
    }
    if ($head[$i] eq 'Chr')   {
        $chr=$i;
    }
    if ($type eq 'SNV') {
        if ($head[$i] eq 'Position')   {
            $pos=$i;
        }
    }
    else    {
        if ($head[$i] eq 'Start')   {
            $pos=$i;
        }
    }
}
#print "$chr\t$pos\t$gene\t$gene_id\n";
#<STDIN>;
print "IGV Link\t$header\tTissue_specificity\tpathway\tGeneCards\tKaviar_Variants\n";
while(my $l = <FH>) {
    chomp $l;
    my @a = split(/\t/,$l);
    #print "$l\n";
    my ($igv,$tissue,$pathway,$genecard,$kaviar);
    chomp $a[$gene_id];
    chomp $a[$gene];
    chomp $a[$chr];
    chomp $a[$pos];
    #print "$a[$gene_id]\t$a[$gene]\t$a[$chr]\t$a[$pos]\n";
    my ($c,$p,$g,$gid);
    $c=$a[$chr];
    $p=$a[$pos];
    $g=$a[$gene];
    $gid=$a[$gene_id];
    $igv = "=Hyperlink(\"http://localhost:60151/goto?locus=$c:$p\",\"$c:$p\")";
    if ($a[$gene_id] ne '-' )   {
        $tissue = "=Hyperlink(\"http://$server/reference/Tissue-Specific.Gene.Expression/html/${gid}.html\",\"link\")";
        $pathway = "=Hyperlink(\"http://$server/reference/metacore_kegg_pathway_files/${gid}_pathways.xls\",\"link\")";
        $genecard = "=Hyperlink(\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=$g\",\"link\")";
    }
    else    {
        $tissue = '-';
        $pathway = '-';
        $genecard = '-';
    }
    $kaviar = "=Hyperlink(\"http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl?chr=$c&pos=$p&frz=${genome}&onebased=1\",\"link\")";
    
    print "$igv\t$l\t$tissue\t$pathway\t$genecard\t$kaviar\n";
    #<STDIN>;
}
close FH;



