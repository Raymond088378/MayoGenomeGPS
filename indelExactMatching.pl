#! /usr/bin/perl
### Raymond Moore
### Jan 2nd 2013

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

=head1 DESCRIPTION
	Documentation Standard: http://perldoc.perl.org/perlpod.html
	This script reads the TempReport/INDEL file, and kGenome VCF
	To identify exact matches in InDels to report out to annotation.
=cut 



## Setup Input Variables
our ($opt_i, $opt_d, $opt_o);

getopt('ido');
if ( (!defined $opt_i) && (!defined $opt_d) && (!defined $opt_o)){
    die ("Usage: $0 \n\t-i [input file] \n\t-d [1000 Genomes file] \n\t-o [output filtered file]\n");
}
else{
	open(DB, "<", $opt_d)||die$!;
	open(IN, "<", $opt_i)||die$!;
	open(OUT, ">", $opt_o)||die$!;

	unless($opt_o =~ /[\\\/]/){ print STDERR "Please Provide the Full Path\n"; exit 1;}
	
	my $h1=<IN>;
	print OUT $h1;
	$h1=<IN>;
	$h1 =~ s/\n//;
	print OUT $h1."\tkGenomeIndelMatch\n";

	my $baseLn = <DB>;
	my $addtoLn = <IN>;
	while( defined $baseLn || defined $addtoLn ){
		if( !defined($addtoLn) ){last;}
		elsif( defined($addtoLn) && !defined($baseLn) ){
			print OUT $addtoLn."\t.\n";
			$addtoLn = <IN>; next;
		}

		#### PASS THROUGH HEADER! ####
		if( $baseLn =~ /^#/){
			$baseLn = <DB>;
			next;
		}
		if( $addtoLn =~ /^#/){
			print OUT $addtoLn;
			$addtoLn = <IN>;
			next;
		}
		
		chomp $baseLn;
		chomp $addtoLn;	
		

		my @baseCol=split(/\t/, $baseLn);
		my @addCol=split(/\t/, $addtoLn);
		## Do they Match? If so...REPORT!
		

		if($baseCol[0] eq $addCol[1] && $baseCol[1] == $addCol[2] && $baseCol[3] eq $addCol[9] && $baseCol[4] eq $addCol[10]){
			push(@addCol, $baseCol[7]);
			print OUT join("\t", @addCol)."\n";
			$baseLn=<DB>;
			$addtoLn=<IN>;
			next;
		}
		### If they don't match....who gets incremented?
		else{
			my $dbChr = numerize($baseCol[0]);
			my $addChr = numerize($addCol[0]);
			if($dbChr < $addChr){
				$baseLn = <DB>; next;
			}
			elsif($addChr < $dbChr){
				$addtoLn = <IN>; next;
			}
			else{
				if($baseCol[1] <= $addCol[2]){
					$baseLn = <DB>; next;
				}
				else{ 
					print OUT $addtoLn."\t.\n";
					$addtoLn = <IN>; next; 
				}
			}
			print OUT join("\t", @addCol)."\t.\n";
		}

	}

	close DB; close IN; close OUT;
}

sub numerize{
	my $str = $_[0];
	my $num=0;
	if($str =~ /chr(\d+|[MXY])/){
		if($1 =~ /X/i){$num=22;}
		elsif($1 =~ /Y/i){$num=24;}
		elsif($1 =~ /M/i){$num=25;}
		else{$num=$1;}
	}
}