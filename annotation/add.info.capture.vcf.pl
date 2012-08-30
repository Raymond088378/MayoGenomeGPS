#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
	if (/^\#\#/) {
		print $_;
		next;
	}

	#Parse the header
	if (/^\#/) {
		print "##INFO=<ID=CAPTURE,Number=1,Type=Integer,Description=\"variant in capture kit ot not\">\n";
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
	@first=(0..6);
	if ( $fields[$#fields] > 0)	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[7]" . ";CAPTURE=1\t$fields[8]\t" . join ("\t",@fields[@samples]) . "\n";  
	}
	elsif ($fields[$#fields] == 0 && $fields[$#fields] =~ /^[+-]?\d+$/)	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[7]" . ";CAPTURE=0\t$fields[8]\t" . join ("\t",@fields[@samples]) . "\n"; 
	}	
	else	{
		print join ("\t",@fields[@first]) . "\t" . "$fields[7]" . ";CAPTURE=1\t$fields[8]\t" . join ("\t",@fields[@samples]) . "\n"; 
	}
}


