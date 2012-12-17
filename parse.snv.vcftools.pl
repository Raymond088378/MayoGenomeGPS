#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
#chr22   16287541        0       .       G       T       G/G     12,0    7       99      0       G/G     24,1    15      99      0


$l=$_;
chomp $l;
@a= split (/\t/,$l);
for ($i=0 ; $i <= $#a ; $i++)	{
	if ($a[$i] eq ".")	{
		$a[$i] = "n/a";
	}
}	
$num_samples=($#a-6+1)/5;	
print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]";
my $init=5;
for (my $i=1; $i <=$num_samples; $i++)	{
	$geno=$a[$init+1];
	$geno =~ s/\///g if ($geno ne "n/a");
	$read=$a[$init+2];
	if ($read eq "n/a")	{
		$ref="n/a";
		$alt="n/a";
		$depth="n/a";
	}
	else{	
	($ref,$alt)=split(/,/,$read);
	$depth=$ref+$alt;
	}
	$qual=$a[$init+4];
	$c2i=$a[$init+5];
	print "\t$geno\t$ref\t$alt\t$depth\t$qual\t$c2i";
$init=$init+5;
	}
print "\n";
}	