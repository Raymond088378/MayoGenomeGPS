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
		   TYPE => 'string',
		   MULTISAMPLE => 'boolean',
		   INPUT_DIR => 'dir',
		   BASE_OUTPUT_DIR => 'string',
		   EMAIL => 'string',
		   SAMPLENAMES => 'csep',
		   CHRINDEX => 'csep',
		   TOOL_INFO => 'file',
		   SAMPLE_INFO => 'file',
		   OUTPUT_FOLDER => 'string',
		   QUEUE => 'string',
		   LQUEUE => 'string',
		   CENTER => 'string',
		   PLATFORM => 'string',
		   GENOMEBUILD => 'string',
		   MARKDUP => 'boolean'
		   );

$runinfomsg = check_variables ($runinfovars, \%runinfofmt, $runinfo);

if ($runinfomsg ne "") {
    print $runinfomsg."\n";
    exit 1;
}
my $analysis=$runinfovars->{TYPE};

if ( ($analysis ne 'whole_genome') && ($analysis ne 'external') && ($analysis ne 'variant') && ($analysis ne 'alignment') && ($analysis ne 'exome') ) {
    print "$runinfo: TYPE=$analysis should be whole_genome, exome, external, variant or alignment\n";
    exit 1;
}

my @samples = split(/:/,$runinfovars->{SAMPLENAMES});

if ( ($runinfovars->{MULTISAMPLE}) && (@samples < 2)) {
    print "$runinfo: MULTISAMPLE=YES and SAMPLENAMES has only one element\n";
    exit 1;
}

my $toolinfo=$runinfovars->{TOOL_INFO};
my ($toolinfovars, $toolinfomsg) = read_file_var($toolinfo);


if ($toolinfomsg ne "") {
    print $toolinfomsg."\n";
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
		   WHOLEGENOME_PATH=>"dir",
		   NOVOALIGN=>"dir",
		   BWA=>"dir",
		   JAVA=>"dir",
		   SAMTOOLS=>"dir",
		   GATK=>"dir",
		   SOMATIC_SNIPER=>"dir",
		   BEDTOOLS=>"dir",
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
		   PERLLIB=>"csep" 
		   );

$toolinfomsg = check_variables ($toolinfovars, \%toolinfofmt, $toolinfo);

if ($toolinfomsg ne "") {
    print $toolinfomsg."\n";
    exit 1;
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


sub read_file_var {
    my ($filename) = @_;
    open INFILE, "<$filename" or die $!;
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
