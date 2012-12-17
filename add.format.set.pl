#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
	if (/^\#\#/) {
		print $_;
		next;
	}

	#Parse the header
	if (/^\#/) {
	   print "##FORMAT=<ID=SET,Number=1,Type=String,Description=\"Source VCF for the merged record in CombineVariants\">\n";

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
	my @format=split(/;/,$fields[7]);
	my ($form,$value);
	for (my $i=0;$i <=$#format;$i++)	{
		if ($format[$i] =~ /^set/)	{
			($form,$value)=split(/=/,$format[$i]);
			last;
		}
	}	
	print join ("\t",@fields[@first]) . "\t" . "$fields[8]" . ":SET" ;
		for ($i=9; $i <= $#nfields; $i++)	{
				if ($fields[$i] ne "./.")	{
					print "\t$fields[$i]:$value";  
				}
				else	{
					print "\t$fields[$i]";
				}
		}

	print "\n";
}
