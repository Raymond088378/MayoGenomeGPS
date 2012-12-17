#! usr/bin/perl

%ucscTracks;
open(IDX, "<", $ARGV[0])||die$!;
while(<IDX>){
	chomp;
	@line=split(/\t/, $_);
	if($line[0] =~ /\_/){ $line[0] =~ s/\_.+$//;}
	$ucscTracks->{$line[0]}->{$.}->{start}=$line[1];
	$ucscTracks->{$line[0]}->{$.}->{end}=$line[2];
	@gene=split(/\#/,$line[3]);
	$ucscTracks->{$line[0]}->{$.}->{gene}=$gene[0];
}

$input=$ARGV[1];
open(IN, "<", $input)|| die $!;
$input =~ s/\.[A-Za-z]{3}$//;
open(OUT, ">", $input."_GENE");
if($input =~ /\_1(\_)*/){
	$h=<IN>;
	$h=~s/\n//g;
	print OUT $h."\tGene\n";
}

while(<IN>){
	chomp;
	@row=split(/\t/, $_);
	while ( my ($key, $value) = each(%{$ucscTracks->{$row[0]}}) ) {
		if( $value->{start} < $row[1] && $row[1] < $value->{end} ){
			push(@genes, $value->{gene});
			#print "$row[1] => ". $value->{gene} ."\n";
		}		
    }
	$geneAnnot = join("|", uniq(@genes) );
	
	print OUT join("\t", (@row, $geneAnnot))."\n";
	undef(@genes);
}



sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
