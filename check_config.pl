#!/usr/bin/perl -w 
my $usage = "check_config.pl runinfo \n";

my $runinfo = shift or die $usage;
my ($runinfovars, $runinfomsg) = read_file_var($runinfo);

# Check for duplicated lines in run_info
if ($runinfomsg ne "") {
    print $runinfomsg."\n";
    exit 1;
}

my %runinfofmt = (
   READLENGTH => 'int',
   PI => 'string',
   TOOL => 'string',
   VERSION => 'real',
   TYPE => 'string',
   MULTISAMPLE => 'boolean',
   INPUT_DIR => 'dir',
   BASE_OUTPUT_DIR => 'string',
   EMAIL => 'string',
   ANALYSIS => 'string',
   SAMPLENAMES => 'csep',
   GROUPNAMES => 'csep',
   CHRINDEX => 'csep',
   LANEINDEX => 'csep',
   LABINDEXES => 'csep',
   TOOL_INFO => 'file',
   SAMPLE_INFO => 'file',
   OUTPUT_FOLDER => 'string',
   QUEUE => 'string',
   CENTER => 'string',
   PLATFORM => 'string',
   GENOMEBUILD => 'string',
   MARKDUP => 'boolean',
   FASTQC => 'boolean',
   FOLDER_FASTQC => 'dir',
   UPLOAD_TABLEBROWSER => 'boolean',
   VARIANT_TYPE => 'string',
   SNV_CALLER => 'string',
   SAMPLEINFORMATION => 'string',
   DELIVERY_FOLDER => 'dir',
   TERTIARY_FOLDER => 'dir'
);

$runinfomsg = check_variables ($runinfovars, \%runinfofmt, $runinfo);

if ($runinfomsg ne "") {
    print $runinfomsg."\n";
    exit 1;
}
my $analysis=$runinfovars->{ANALYSIS};

if ( ($analysis ne 'mayo') && ($analysis ne 'external') && ($analysis ne 'variant') && ($analysis ne 'alignment') && ($analysis ne 'realignment') && ($analysis ne 'realign-mayo') && ($analysis ne 'ontarget') && ($analysis ne 'annotation') ) {
    print "$runinfo: TYPE=$analysis should be mayo, external, variant , realignment , realign-mayo, ontarget, annotation or alignment\n";
    exit 1;
}

my $tool=$runinfovars->{TYPE};

if ( ($tool ne 'whole_genome') && ($tool ne 'exome'))	{
    print "$runinfo : TOOL= $tool should be whole_genome or exome\n";
}		

my $input=$runinfovars->{INPUT_DIR};
my @samples= split(/:/,$runinfovars->{SAMPLENAMES});
my @groups = split (/:/,$runinfovars->{GROUPNAMES});
## min samples needed per group is 2
if (($runinfovars->{MULTISAMPLE} eq 'YES') &&(@samples == 1))	{
        print "$runinfo: MULTISAMPLE=YES and SAMPLENAMES has only one element\n";
}

if (($runinfovars->{MULTISAMPLE} eq 'YES') &&(@groups == 0))	{
        print "$runinfo: MULTISAMPLE=YES and no groups defined\n";
}

my $toolinfo=$runinfovars->{TOOL_INFO};
my $sampleinfo=$runinfovars->{SAMPLE_INFO};
die "Tool info file is empty\n\n" if (-z $toolinfo || ! -e $toolinfo) ;
die "sample info file is empty\n\n" if (-z $sampleinfo || ! -e $sampleinfo) ;

my ($toolinfovars, $toolinfomsg) = read_file_var($toolinfo);
if ($toolinfomsg ne "") {
    print $toolinfomsg."\n";
    exit 1;
}

my ($sampleinfovars, $sampleinfomsg) = read_files_var($sampleinfo);

if ($sampleinfomsg ne "") {
    print $sampleinfomsg."\n";
    exit 1;
}

my ($samplenames, $samplenamesmsg)= read_sample_names($sampleinfo);
if ($samplenamesmsg ne "") {
    print $samplenamesmsg."\n";
    exit 1;
}

my ($pairsamplenames,$pairnamesmsg) = read_sample_pair_names($sampleinfo);
if ($pairnamesmsg ne "") {
    print $pairnamesmsg."\n";
    exit 1;
}



my %toolinfofmt = (
	REF_GENOME=>"file",
	SPLIT_GENOME=>"dir",
	GAP_GENOME=>"file",
	BWA_REF=>"file",
	BLAT_REF=>"file",
	BIT_DIR=>"dir",
	KG_INDELS_VCF=>"file",
	DBSNP_VCF=>"file",
	HAPMAP_VCF=>"file",
	OMNI_VCF=>"file",
	NOVO_REF=>"file",
	dbSNP_REF=>"file",
	KGENOME_REF=>"file",
	dbSNP_SNV_rsIDs=>"file",
	MASTER_GENE_FILE=>"file",
	MASTER_ENTREZ_FILE=>"file",
	CAPTUREKIT=>"file",
	ONTARGET=>"file",
	BLACKLISTED=>"file",
	MAPABILITY=>"file",
	REPEATREGIONS=>"file",
	miRbase=>"file",
	SNP_SR=>"file",
	SNP_CS=>"file",
	BLACKLIST_SV=>"file",
	METACORE_PATHWAY=>"dir",
	TISSUE_SPECIFIC=>"dir",
	KGENOME=>"dir",
	HAPMAP=>"dir",
	CODON_REF=>"file",
	WHOLEGENOME_PATH=>"dir",
	NOVOALIGN=>"dir",
	BWA=>"dir",
	JAVA=>"dir",
	SAMTOOLS=>"dir",
	GATK=>"dir",
	SOMATIC_SNIPER=>"dir",
	BEDTOOLS=>"dir",
	SIFT=>"dir",
	VCFTOOLS=>"dir",
	TABIX=>"dir",
	CIRCOS=>"dir",
	PICARD=>"dir",
	FASTQC=>"dir",
	CNVNATOR=>"dir",
	ROOTLIB=>"dir",
	CREST=>"dir",
	CAP3=>"dir",
	BLAT=>"dir",
	BREAKDANCER=>"dir",
	PERL_BREAKDANCER=>"dir",
	PERLLIB_BREAKDANCER=>"dir",
	EMIT_ALL_SITES=> 'boolean',
	TARGETTED=> 'boolean',
	MATLAB=>"dir",
	SEGSEQ=>"dir",
	SCRIPT_PATH=>"dir",
	CNVNATOR_BINSIZE=>"int",
	PCT_READS_SEGSEQ=>"real",
	MINFOLD=>"real",
	MAXFOLD=>"real",
	DISTGAP=>"int",
	BLAT_PORT=>"int",
	BLAT_SERVER=>"string",
	HTTP_SERVER=>"string",
	PERLLIB=>"csep",
	PERLLIB_VCF=>"csep",
	PERL_CIRCOS=>"csep"
);

$toolinfomsg = check_variables ($toolinfovars, \%toolinfofmt, $toolinfo);

$sampleinfomsg = check_files ($sampleinfovars, $sampleinfo, $input );


if ($toolinfomsg ne "") {
    print $toolinfomsg."\n";
    exit 1;
}

if ($sampleinfomsg ne "") {
    print $sampleinfomsg."\n";
    exit 1;
}

$sampleinfomsg = check_samples($samplenames,@samples,@groups);
if ($sampleinfomsg ne "") {
    print $sampleinfomsg."\n";
    exit 1;
}
$sampleinfomsg = check_pairs($pairsamplenames,@samples);
if ($sampleinfomsg ne "") {
    print $sampleinfomsg."\n";
    exit 1;
}


sub check_pairs{
	my ($svars,@samples) = @_; 
    my $errmsg="";
	foreach my $key (keys %{$svars})	{
	my @sam=split(/\t/,$svars->{$key}); 
	foreach (@sam)	{
		if (! grep /$_/,@samples) {
			$errmsg .= "$_ : does not match run_info for pairs \n";
	    }
    }
	}
    return $errmsg;	
}
sub check_samples{
    my ($svars,@samples,@groups) = @_; 
    my $errmsg="";
	foreach my $key (keys %{$svars})	{
	my @sam=split(/\t/,$svars->{$key}); 
	foreach (@sam)	{	
		if (! grep /$_/,@samples || ! grep /$_/,@groups) {
			$errmsg .= "$_: does not match run_info\n";
	    }
    }
	}
    return $errmsg;	
}


sub check_files	{
    my ($rvars,$fname,$input) = @_;
    my $errmsg="";
    foreach my $key (keys %{$rvars})	{
	my @files=split(/\t/,$rvars->{$key});
	foreach (@files)	{
	    my $name="$input/$_";
	    if (! -e $name)	{
		$errmsg .= "$name: does not exist\n";
	    }		
	}	
    }
    return $errmsg;	
}		


sub check_variables {
    my ($rvars, $rformat, $fname) = @_;
    my $errmsg="";

    foreach my $var (keys %{$rformat}) {

	if (! exists $rvars->{$var}) {
	    $errmsg .= '$fname: $var is not defined\n';
	}
	else {
	    my $value = $rvars->{$var};
	    my $type = $rformat->{$var};

	    if ( ($type eq 'file') && (!-e $rvars->{$var})) {
		$errmsg .= "$fname: $var, file $value does not exist\n";
	    }

	    if ( ($type eq 'dir') && (!-d $rvars->{$var})) {
		$errmsg .= "$fname: $var, directory $value does not exist\n";
	    }

	    if ( ($type eq 'boolean') && ( ($value ne 'YES') && ($value ne 'NO'))) {
		$errmsg .= "$fname: $var, $value should be YES or NO\n";
	    }

	    if ( $type eq 'csep') {
		chomp $value;
		my @fields = split (/:/,$value);
		foreach my $field (@fields) {
		    $errmsg .= "$fname: $var, $value should not have whitespace (\\t,\" \", \\n), it is colon separated \n" 
			if ($field =~ /.*\s+.*/);
		}
	    }
	}
    }
}

## check for sample names (misspled)

sub read_sample_names{
    my ($filename) = @_;
    open INFILE, "<$filename" or die "$filename $!\n";
    my %variables;
    my %names;
    my $errmsg="";
 
    LINE:while (<INFILE>) {
		next LINE if /^#|^\s/;
		chomp;

		my ($var, $value) = split (/=/,$_);
		if ($var =~ /:/)	{
			my ($type, $sam) = split(/:/,$var);	
			if (exists $variables{$var}) {
				$errmsg.="$var in $filename defined twice\n";
			}
			$variables{$var}=$value;
			$names{$sam}=$sam;
		}
    }
    return (\%names, $errmsg)	
}

sub read_sample_pair_names{
    my ($filename) = @_;
    open INFILE, "<$filename" or die "$filename $!\n";
    my %variables;
    my %names;
    my $errmsg="";
 
    LINE:while (<INFILE>) {
		next LINE if /^#|^\s/;
		chomp;

		my ($var, $value) = split (/=/,$_);
		if ($var !~ /:/)	{
			my @sample= split(/\t/,$value);	
			foreach my $sam (@sample)	{
			    if (exists $names{$sam} ) {
				$errmsg.="$filename:$sam in group $var is defined twice\n";
			    }
			    $names{$sam}=$sam;
			}
		}
    }
	return (\%names, $errmsg)	
}


sub read_files_var{
    my ($filename) = @_;
    open INFILE, "<$filename" or die "$filename $!\n";
    my %variables;
    my $errmsg="";
 
    LINE:while (<INFILE>) {
	next LINE if /^#|^\s/;
	chomp;

	my ($var, $value) = split (/=/,$_);
	if ($var =~ /:/)	{
	    if (exists $variables{$var}) {
		$errmsg.="$var in $filename defined twice\n";
	    }
	    $variables{$var}=$value;
	}
    }
    return (\%variables, $errmsg)
}

sub read_file_var {
    my ($filename) = @_;
    open INFILE, "<$filename" or die "$filename $!\n";
    my %variables;
    my $errmsg="";
 
    LINE:while (<INFILE>) {
	next LINE if /^#|^\s/;
	chomp;

	my ($var, $value) = split (/=/,$_);
	if (exists $variables{$var}) {
	    $errmsg.="$var in $filename defined twice\n";
	}
	$variables{$var}=$value;
    }
    return (\%variables, $errmsg)
}


