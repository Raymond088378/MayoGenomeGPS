#!/usr/local/biotools/perl/5.10.0/bin/perl
use strict;
use warnings;
use Getopt::Long;


my($input, $snv, $indel, $sample, $chr, $filter, $help);
GetOptions("in|i=s"	=> \$input,
	"snv|v:s"	=> \$snv,
	"indel|l:s"	=> \$indel,
	"sample|s:s"	=> \$sample,
	"chr|c:s"	=> \$chr,
	"filter|f"	=> \$filter,
	"help|h|?|"	=> \&help);

if(not $input){
	print "Missing input file!\n";
	help();
	exit 1;
}


open OUT, ">$snv" or die "opening $snv : $! \n";
open OUT1, ">$indel" or die "opening $indel : $! \n";

open FH, "$input" or die " opening $input : $! \n";
my $s;

while (my $l = <FH>)	{
	chomp $l;
	if ($l =~ /^##/)	{
		print OUT "$l\n";
		print OUT1 "$l\n";
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
			print OUT join ("\t",@a[@cols]) . "\n";
			print OUT1 join ("\t",@a[@cols]) . "\n";
		}
		else	{
			print OUT "$l\n";		
			print OUT1 "$l\n";
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
				else	{
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
		else	{
			### ref=3 alt=4 
			my @alt=split(/,/,$a[4]);
			my $ref=length($a[3]);
			my $len=1;
                        for(my $i=0; $i <= $#alt; $i++) {
                                $len = length($alt[$i]);
                        }
                        if ($ref == 1 && $len == 1)  {
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
                        else    {
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
close FH;		
close OUT;
close OUT1;

sub help{
	print "DESCRIPTION:
	vcf_to_variant_vcf.pl splits vcf into snvs and indels

USAGE:
	vcf_to_variant_vcf.pl -i input.vcf -v snv.vcf -i indel.vcf

OPTIONS:
	--in,-i		Path to snvmix file. Required parameter. Input can be gzipped as long 
			as the -z flag is also used. Input files should not have a header line.

	--snv,-v 	path to the output SNV VCF output file.

	--indel,-l	path to the output INDEL VCF output file.

	--sample,-s 	sample name (optional)
        
	--chr,-c    	chromsome (optional)
	
	--filter,-f	if -f is there then only output PASS filter in vcf

	--help,-h,-?	Display this documentation.

";
}




