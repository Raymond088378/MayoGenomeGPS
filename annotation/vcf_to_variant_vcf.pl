#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;
use Getopt::Long;


my($input, $snv, $indel, $sample, $chr, $filter, $type, $help);
GetOptions("in|i=s"	=> \$input,
	"snv|v:s"	=> \$snv,
	"indel|l:s"	=> \$indel,
	"sample|s:s"	=> \$sample,
	"chr|c:s"	=> \$chr,
	"filter|f"	=> \$filter,
	"type|t:s" => \$type,
	"help|h|?|"	=> \&help);

if(not $input){
	print "Missing input file!\n";
	help();
	exit 1;
}
if (not $type) {
	print "missing the type fied! \n";
	help();
	exit 1;
}	

open OUT, ">$snv" or die "opening $snv : $! \n" if ($type eq "both" || $type eq "snv");
open OUT1, ">$indel" or die "opening $indel : $! \n" if ($type eq "both" || $type eq "indel");


open FH, "$input" or die " opening $input : $! \n";
my $s;

while (my $l = <FH>)	{
	chomp $l;
	if ($l =~ /^##/)	{
		print OUT "$l\n" if ($type eq "both" || $type eq "snv");
		print OUT1 "$l\n" if ($type eq "both" || $type eq "indel");
	}
	elsif ($l =~ /^#CHROM/)	{
		my @a = split(/\t/,$l);
		if (defined $sample)	{
			for(my $i =0 ; $i <=$#a; $i++)	{
				if ($a[$i] eq "$sample")	{
					$s=$i;
					last;
				}
			}
			my @cols=(0,1,2,3,4,5,6,7,8,$s);
			print OUT join ("\t",@a[@cols]) . "\n" if ($type eq "both" || $type eq "snv");
			print OUT1 join ("\t",@a[@cols]) . "\n" if ($type eq "both" || $type eq "indel");
		}
		else	{
			print OUT "$l\n" if ($type eq "both" || $type eq "snv");	
			print OUT1 "$l\n" if ($type eq "both" || $type eq "indel");
		}	

	}
	else	{
		my @a = split(/\t/,$l);
		
		my @cols=(0,1,2,3,4,5,6,7,8,$s);
		if (defined $chr)	{
			if ($a[0] eq "$chr")	{
				### ref=3 alt=4	
				my $ref=length($a[3]);
				my @alt=split(/,/,$a[4]);
				my $len=1;
				for(my $i=0; $i <= $#alt; $i++)	{
					$len = length($alt[$i]);
				}
				if ($ref == 1 && $len == 1)	{
					if ($type eq "both" || $type eq "snv")	{
						if (defined $sample)	{
							if (defined $filter)	{
								if ($a[6] eq 'PASS')	{
									next if ($a[$s] eq "\.");
									print OUT join ("\t",@a[@cols]) . "\n";
								}
							}
							else	{
								next if ($a[$s] eq "\.");
								print OUT join ("\t",@a[@cols]) . "\n";
							}	
						}
						else	{
							if (defined $filter)	{
								if ($a[6] eq 'PASS')	{	
									print OUT "$l\n";
								}
							}
							else	{
								print OUT "$l\n";
							}
						}
					}
				}
				else	{
					if ($type eq "both" || $type eq "indel")	{
						if (defined $sample)    {
							if (defined $filter)	{
								if ($a[6] eq 'PASS')	{
									next if ($a[$s] eq "\.");
									print OUT1 join ("\t",@a[@cols]) . "\n";
								}
							}
							else	{
								next if ($a[$s] eq "\.");
								print OUT1 join ("\t",@a[@cols]) . "\n";	
							}
						}
						else	{
							if (defined $filter)	{
								if ($a[6] eq 'PASS')	{
									print OUT1 "$l\n";
								}
							}
							else	{
								print OUT1 "$l\n";	
							}
						}	
					}
				}			
			}	
		}
		else	{
			### ref=3 alt=4 
			my @alt=split(/,/,$a[4]);
			my $ref=length($a[3]);
			my $len=1;
			for(my $i=0; $i <= $#alt; $i++) {
				$len = length($alt[$i]);
			}
			if ($ref == 1 && $len == 1)  {
				if ($type eq "both" || $type eq "snv")	{
					if (defined $sample)    {
						if (defined $filter)	{
							if ($a[6] eq 'PASS')	{
								next if ($a[$s] eq "\.");
								print OUT join ("\t",@a[@cols]) . "\n";
							}
						}
						else	{
							next if ($a[$s] eq "\.");
							print OUT join ("\t",@a[@cols]) . "\n";
						}
					}
					else    {
						if (defined $filter)	{
							if ($a[6] eq 'PASS')	{
								print OUT "$l\n";
							}
						}
						else	{
							print OUT "$l\n";
						}
					}
				}
			}
			else    {
				if ($type eq "both" || $type eq "indel")	{
					if (defined $sample)    {
						if (defined $filter)	{
							if ($a[6] eq 'PASS')	{
								next if ($a[$s] eq "\.");
								print OUT1 join ("\t",@a[@cols]) . "\n";
							}
						}
						else	{
							next if ($a[$s] eq "\.");
							print OUT1 join ("\t",@a[@cols]) . "\n";
						}	
					}
					else    {
						if (defined $filter)	{
							if ($a[6] eq 'PASS')	{
								print OUT1 "$l\n";
							}
						}
						else	{
							print OUT1 "$l\n";
						}
					}
				}
			}
		}
	}
}
close FH;		
close OUT if ($type eq "both" || $type eq "snv");
close OUT1 if ($type eq "both" || $type eq "indel");

sub help{
	print "DESCRIPTION:
	vcf_to_variant_vcf.pl splits vcf into snvs and indels

USAGE:\n\tvcf_to_variant_vcf.pl -i input.vcf -v snv.vcf -i indel.vcf -t both
\nOPTIONS:\n\t--in,-i\tPath to vcf file. Required parameter. 
\n\t--snv,-v\tpath to the output SNV VCF output file.
\n\t--indel,-l\tpath to the output INDEL VCF output file.
\n\t--sample,-s\tsample name (optional)
\n\t--chr,-c\tchromsome (optional)
\n\t--type,-t\ttype of variant (both,snv,indel). Required paramter
\n\t--filter,-f\tif -f is there then only output PASS filter in vcf
\n\t--help,-h,-?\tDisplay this documentation.

";
exit 1;
}




