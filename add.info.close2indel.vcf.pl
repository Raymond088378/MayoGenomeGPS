#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
	if (/^\#\#/) {
		print $_;
		next;
    }

	#Parse the header
	if (/^\#/) {
       print "##FORMAT=<ID=C2I,Number=1,Type=Integer,Description=\"if a snp is close to indel for this sample\">\n";

	chomp;
       my $line = $_;
       $line =~ s/\#//;
       @nfields = split (/\t/, $line) ;
       @samples = (9..$#nfields);
       print "#$line\n";
       next;
	}
	
	my $line = $_;
	chomp $line;
	@fields = split (/\t/,$line);
	@first=(0..7);
	if ( $fields[$#fields] ==1 )	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[8]" . ":C2I" ;
		for ($i=9; $i <= $#nfields; $i++)	{
			if ($fields[$i] =~ m/^(.:)+/){
				print "\t$fields[$i]:.";
			}	
			elsif ($fields[$i] ne "./." && $fields[$i] ne "." )	{
				print "\t$fields[$i]:1";  
			}
			else	{
				print "\t$fields[$i]";
			}	
		}
	}
	elsif ($fields[$#fields] == 0 && $fields[$#fields] =~ /^[+-]?\d+$/)	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[8]" . ":C2I" ; 
		for ($i=9; $i <= $#nfields; $i++)	{
			if ($fields[$i] =~ m/^(.:)+/){
				print "\t$fields[$i]:.";
			}	
			elsif ($fields[$i] ne "./." && $fields[$i] ne ".")	{
			print "\t$fields[$i]:0";  
			}
			else	{
				print "\t$fields[$i]";
			}	
		}
	}	
	else	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[8]" . ":C2I" ; 
		for ($i=9; $i <= $#nfields; $i++)	{
			if ($fields[$i] =~ m/^(.:)+/){
				print "\t$fields[$i]:.";
			}	
			elsif ($fields[$i] ne "./." && $fields[$i] ne ".")	{
			print "\t$fields[$i]:0";  
			}
			else	{
				print "\t$fields[$i]";
			}	
		}
	}
	print "\n";
}
