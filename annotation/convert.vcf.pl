#!/usr/local/biotools/perl/5.10.0/bin/perl

while(<>)	{
	if (/^\#\#/) {
		print $_;
		next;
	}

	#Parse the header
	if (/^\#/) {
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
	@first=(0..8);
	print join ("\t",@fields[@first]);
	@format= split (/:/,$fields[8]);
	$num_format=$#format;
	for ($i=9;$i <= $#nfields; $i++)	{
		if ($fields[$i] eq "./.")	{
			print "\t." . ":." x $num_format; 
		}
		else {
			@a=split(/:/,$fields[$i]);
			if ($#a < $num_format)	{
				$val=$num_format-$#a;
				print "\t$fields[$i]" . ":." x $val;
			}	
			else	{
				print "\t$fields[$i]";
			}
		}
	}
	print "\n";	
}