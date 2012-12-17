#! /usr/bin/perl
use Spreadsheet::ParseExcel;

$input=$ARGV[0];
$features=$ARGV[1];
$tmpDir=$ARGV[2];
$splitTot=$ARGV[3];

##### Obtain desired Features #####
my $oBook = Spreadsheet::ParseExcel::Workbook->Parse($features);
my $sh = ${$oBook->{Worksheet}}[0];
if(!defined $sh->{MaxRow}){die "EMPTY WORKSHEET!!\n"};

%treat;
for ($row = 0; $row <= $sh->{MaxRow}; $row++) {
	if(defined($sh->{Cells}[$row][2])){
		next if( $sh->{Cells}[$row][2]->Value =~ /Require(d)*/i );
	}
	if( $sh->{Cells}[$row][0]->Value =~ /Y(es)*/i ){
		if(defined($sh->{Cells}[$row][2])){
			$treat{$sh->{Cells}[$row][1]->Value}=$sh->{Cells}[$row][2]->Value;
		}
		else{ $treat{$sh->{Cells}[$row][1]->Value}="."}
	#	print $sh->{Cells}[$row][1]->Value."\n";
	}
}	

#### Read & Parse Input Tabbed File ####
open(IN, "<", $input)|| die $!;
@head1=split(/\t/, <IN>);
@header=split(/\t/, <IN>);
@requiredPositions; @treatPositions; @treatTitle;
for($i=0;$i<$#header;$i++){
	## grab required columns --- input file changes...needs to be position independant.
	if($header[$i] =~ /^chr(omosome)?$/i){$requiredPositions[0]=$i;}
	elsif($header[$i] =~ /^pos(ition)?$/i){$requiredPositions[1]=$i;}
	elsif($header[$i] =~ /^(dbSNP\d+|rsid)$/i){$requiredPositions[2]=$i;}
	elsif($header[$i] =~ /^ref$/i){$requiredPositions[3]=$i;}
	elsif($header[$i] =~ /^alt$/i){$requiredPositions[4]=$i;}
	
	### Grab columns selected in feature selection by user
	elsif( exists($treat{$header[$i]}) ){ 
		push(@treatPositions, $i); 
		push(@treatTitle, $header[$i]." (".mainHeader($i).")"); 
	}
}
$varfile=1;
open(OUT, ">", $tmpDir."/VARFILE_".$varfile)|| die $!;
print OUT "#Chr\tPosition\tdbSNP\tRef\tAlt\n";
open(OUT_T, ">", $tmpDir."/VARFILE_".$varfile."_TREAT")|| die $!;
print OUT_T join("\t", @treatTitle)."\n";

$l=1;
while(<IN>){
	chomp;
	@line=split(/\t/, $_);
	print OUT $line[$requiredPositions[0]]."\t".$line[$requiredPositions[1]]."\t".$line[$requiredPositions[2]]."\t";
	print OUT $line[$requiredPositions[3]]."\t".$line[$requiredPositions[4]]."\n";
	
	foreach $nt (@treatPositions){ push(@treatOut, $line[$nt]) }
	print OUT_T join("\t", @treatOut)."\n";
	undef(@treatOut);
	
	$l++;	
	if($l % $splitTot == 0){
##		print "Line $l\n"; ### Print out the number of lines per file
		close(OUT); close(OUT_T);
		$varfile++;
		open(OUT, ">", $tmpDir."/VARFILE_".$varfile)|| die $!;
		open(OUT_T, ">", $tmpDir."/VARFILE_".$varfile."_TREAT")|| die $!;
	}
	
}

close(OUT); close(OUT_T);

### Reverse look up last title in Super Header
sub mainHeader{
	for($k=$_[0];$k>=0;$k--){
		if($head1[$k] ne ""){
			return $head1[$k];
		}
	}
	return ""
}
