#!/usr/local/biotools/perl/5.10.0/bin/perl


use strict;
#use warnings;

my $uniprot=shift @ARGV;
my $prediction=shift @ARGV;

open PRED, "$prediction" or die "can not open $prediction : $! \n";
my $head=<PRED>;
## 1,2,3,4,11
#o_acc, o_pos, o_aa1, o_aa2, prediction
my @header=split('\s+',$head);

my ($acc,$pos,$aa1,$aa2,$pre);
for(my $i=0;$i <=$#header;$i++) {
    if ($header[$i] eq "#o_acc")    {
        $acc=$i;
    }
    if ($header[$i] eq 'o_pos')    {
        $pos=$i;
    }
    if ($header[$i] eq 'o_aa1')    {
        $aa1=$i;
    }
    if ($header[$i] eq 'o_aa2')    {
        $aa2=$i;
    }
    if ($header[$i] eq 'prediction')    {
        $pre=$i;
    }
}
#print "$acc,$pos,$aa1,$aa2,$pre\n";
#<STDIN>;
my %hash_p=();
while(my $l = <PRED>)   {
    chomp $l;
    my @a=split(/\t/,$l);
    my $lens=$#a;
    for(my $i = 0; $i <= $lens;$i++)        {
        if (length($a[$i]) == 0)        {
            $a[$i]='-';
        }
    }
    if ($lens < $#header)      {
        for (my $j=$lens+1;$j <=$#header; $j++)    {
            $a[$j]='-';
        }
    }
    $a[$acc] =~ m/(\S+)/;
    my $accession=$1;
    $a[$pos] =~ m/(\S+)/;
    my $position=$1;
    $a[$aa1] =~ m/(\S+)/;
    my $amino1=$1;
    $a[$aa2] =~ m/(\S+)/;
    my $amino2=$1;
    my $id=$accession."_".$position."_".$amino1."_".$amino2;
    #print "$id\n";
    #<STDIN>;
    $a[$pre] =~ m/(\D+)/; 
    my $prediction=$1;
    $hash_p{$id}=$prediction;
}
close PRED;

open UNI, "$uniprot" or die "can not open $uniprot : $!\n";
$head=<UNI>;
my $chrpos;
@header=split('\s+',$head);
for(my $i=0;$i <=$#header;$i++) {
    if ($header[$i] eq 'spacc')    {
        $acc=$i;
    }
    if ($header[$i] eq 'cdnpos')    {
        $pos=$i;
    }
    if ($header[$i] eq 'aa1')    {
        $aa1=$i;
    }
    if ($header[$i] eq 'aa2')    {
        $aa2=$i;
    }
    if ($header[$i] eq '#snp_pos')    {
        $chrpos=$i;
    }
}
while (my $l = <UNI>)   {
    chomp $l;
    my @a=split(/\t/,$l);
    my $lens=$#a;
    for(my $i = 0; $i <= $lens;$i++)        {
        if (length($a[$i]) == 0)        {
            $a[$i]='-';
        }
    }
    if ($lens < $#header)      {
        for (my $j=$lens+1;$j <=$#header; $j++)    {
            $a[$j]='-';
        }
    }
    $a[$acc] =~ m/(\S+)/;
    my $accession=$1;
    $a[$pos] =~ m/(\S+)/;
    my $position=$1;
    $a[$aa1] =~ m/(\S+)/;
    my $amino1=$1;
    $a[$aa2] =~ m/(\S+)/;
    my $amino2=$1;
    my $id=$accession."_".$position."_".$amino1."_".$amino2;
    $a[$chrpos] =~ m/(\S+)/;
    print "$1\t$a[$acc]\t$hash_p{$id}\n" if exists $hash_p{$id};
}
close UNI;




