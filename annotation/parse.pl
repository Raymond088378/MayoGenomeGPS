#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
#chr22   17590755        1       .       AGT     A       36,19   55      47,15   63
#chr22   17640399        1       .       G       GGGC    45,32   75      73,37   107

$l=$_;
chomp $l;
@a= split (/\t/,$l);
for ($i=0 ; $i <= $#a ; $i++)	{
	if ($a[$i] eq ".")	{
		$a[$i] = "n/a";
	}
}	

$stop=0;
$bases=0;
$start=$a[1];
if (length($a[4]) < length($a[5]))	{
	$stop=$start;
	$bases=length(substr($a[5],1));
}
else	{
	$bases=length(substr($a[4],1));
	$stop=$start+$bases;
}	
$num_samples=($#a-6+1)/2;

	
print "$a[0]\t$start\t$stop\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$bases";
my $init=5;
for (my $i=1; $i <=$num_samples; $i++)	{
	$read=$a[$init+1];
	if ($read eq "n/a")	{
		$read="n/a";
		$depth="n/a";
	}
	else{	
	($ref,$alt)=split(/,/,$read);
	$depth=$ref+$alt;
	}
	print "\t$alt\t$depth";
$init=$init+2;
	}
print "\n";
}	