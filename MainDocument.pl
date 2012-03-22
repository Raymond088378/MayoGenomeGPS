
use strict;
# use warnings;
use Getopt::Std;
our ($opt_r,$opt_p);
print "RAW paramters: @ARGV\n";
getopt('rp');
if ( (!defined $opt_r) && (!defined $opt_p) ) {
        die ("Usage: $0 \n\t-r [ un info file ] \n\t-p [ ath output folder ] \n");
}
else    {
	my ($run_info) = $opt_r;chomp $run_info;
	my ($path) = $opt_p;chomp $path;
	my $desc = $path."/StatisticsDescription.html";
	open DESC, ">$desc" or die "can not open $desc : $!\n";
	my @line=split(/=/,`perl -ne "/^DISEASE/ && print" $run_info`);
	my $disease=$line[$#line];chomp $disease;
	@line=split(/=/,`perl -ne "/^DATE/ && print" $run_info`);
	my $date=$line[$#line];chomp $date;
	@line=split(/=/,`perl -ne "/^SAMPLENAMES/ && print" $run_info`);
	my $sampleNames=$line[$#line];chomp $sampleNames;
	my @sampleArray = split(/:/,$sampleNames);
	@line=split(/=/,`perl -ne "/^LANEINDEX/ && print" $run_info`);
	my $laneNumbers=$line[$#line];chomp $laneNumbers;
	my @laneArray = split(/:/,$laneNumbers);
	@line=split(/=/,`perl -ne "/^PAIRED/ && print" $run_info`);
	my $paired=$line[$#line];chomp $paired;
	@line=split(/=/,`perl -ne "/^GENOMEBUILD/ && print" $run_info`);
	my $GenomeBuild=$line[$#line];chomp $GenomeBuild;
	@line=split(/=/,`perl -ne "/^TYPE/ && print" $run_info`);
	my $tool=$line[$#line];chomp $tool;
	$tool=lc($tool);
	@line=split(/=/,`perl -ne "/^ANALYSIS/ && print" $run_info`);
	my $analysis=$line[$#line];chomp $analysis;
	@line=split(/=/,`perl -ne "/^SAMPLEINFORMATION/ && print" $run_info`);
	my $sampleinfo=$line[$#line];chomp $sampleinfo;
	@line=split(/=/,`perl -ne "/^OUTPUT_FOLDER/ && print" $run_info`);
	my $run_num=$line[$#line];chomp $run_num;
	@line=split(/=/,`perl -ne "/^TOOL_INFO/ && print" $run_info`);
	my $tool_info=$line[$#line];chomp $tool_info;
	@line=split(/=/,`perl -ne "/^LABINDEXES/ && print" $run_info`);
    my $labindex=$line[$#line];chomp $labindex;
	my @IndexArray = split(/:/,$labindex);
	 @line=split(/=/,`perl -ne "/^dbSNP_SNV_rsIDs/ && print" $tool_info`);
	my $dbsnp_file=$line[$#line];chomp $dbsnp_file;
	$dbsnp_file =~ m/.+dbSNP(\d+)/;
	my $dbsnp_v = $1;
	@line=split(/=/,`perl -ne "/^WHOLEGENOME_PATH/ && print" $tool_info`);
	my $script_path=$line[$#line];chomp $script_path;
	my ($read_length, $variant_type, $target_region, $SNV_caller, $Aligner, $ontarget, $fastqc, $fastqc_path, $server, $upload_tb );
	@line=split(/=/,`perl -ne "/^UPLOAD_TABLEBROWSER/ && print" $run_info`);
	my $upload_tb=$line[$#line];chomp $upload_tb;
	@line=split(/=/,`perl -ne "/^TABLEBROWSER_PORT/ && print" $run_info`);
	my $port=$line[$#line];chomp $port;
	@line=split(/=/,`perl -ne "/^TABLEBROWSER_HOST/ && print" $run_info`);
	my $host=$line[$#line];chomp $host;
	@line=split(/=/,`perl -ne "/^READLENGTH/ && print" $run_info`);
	$read_length=$line[$#line];chomp $read_length;
	@line=split(/=/,`perl -ne "/^MULTISAMPLE/ && print" $run_info`);
	my $multi=$line[$#line];chomp $multi;
	
	if ( ( $analysis eq 'external' ) || ( $analysis eq 'variant' ) || ($analysis eq 'mayo') || ($analysis eq 'realignment') )	{
		@line=split(/=/,`perl -ne "/^ALIGNER/ && print" $run_info`);
		$Aligner=$line[$#line];chomp $Aligner;
		@line=split(/=/,`perl -ne "/^VARIANT_TYPE/ && print" $run_info`);
		$variant_type=$line[$#line];chomp $variant_type;
		@line=split(/=/,`perl -ne "/^SNV_CALLER/ && print" $run_info`);
		$SNV_caller=$line[$#line];chomp $SNV_caller;
		if ($tool eq 'exome')	{
			@line=split(/=/,`perl -ne "/^CAPTUREKIT/ && print" $tool_info`);
			$ontarget=$line[$#line];chomp $ontarget;
			$target_region=`awk '{sum+=\$3-\$2+1; print sum}' $ontarget | tail -1`;chomp $target_region;
		}
		elsif ($tool eq 'whole_genome'){
			$ontarget=$path."/bed_file.bed";
			$target_region=`awk '{sum+=\$3-\$2+1; print sum}' $ontarget | tail -1`;chomp $target_region;
		}
		@line=split(/=/,`perl -ne "/^FASTQC/ && print" $run_info`);
		$fastqc=$line[$#line];chomp $fastqc;
		@line=split(/=/,`perl -ne "/^HTTP_SERVER/ && print" $tool_info`);
		$server=$line[$#line];chomp $server;	
	}
	print "Run : $run_num \n";
	print "Result Folder : $path \n";
	my $numbers=$path."/numbers";
	print "Disease: $disease\n";
	print "Date: $date\n";
	print "Samples: $sampleNames\n";
	print "lanes: $laneNumbers\n";
	print "Region: $target_region\n";
	print "PairedEnd: $paired\n";
	print "Type of Tool: $tool\n";
	print "ReadLength: $read_length\n";
	print "AnalysisType: $analysis\n";
	
	sub CommaFormatted
	{
		my $delimiter = ','; # replace comma if desired
		my($n,$d) = split /\./,shift,2;
		my @a = ();
		while($n =~ /\d\d\d\d/)
		{
			$n =~ s/(\d\d\d)$//;
			unshift @a,$1;
		}
		unshift @a,$n;
		$n = join $delimiter,@a;
		$n = "$n\.$d" if $d =~ /\d/;
		return $n;
	}
	
	my $output = "$path/Main_Document.html";
	print "Generating the Document... \n";
	open (OUT,">$output");
	print OUT "<html>"; 
	print OUT "<head>"; 
	print OUT "<title>${tool} Analysis Main Document for $run_num</title>"; 
	print OUT "<style>";
	print OUT "table.helpT	{ text-align: center;font-family: Verdana;font-weight: normal;font-size: 11px;color: #404040;width: auto;
	background-color: #fafafa;border: 1px #6699CC solid;border-collapse: collapse;border-spacing: 0px; }
	td.helpHed	{ border-bottom: 2px solid #6699CC;border-left: 1px solid #6699CC;background-color: #BEC8D1;text-align: left;
	text-indent: 5px;font-family: Verdana;font-weight: bold;font-size: 11px;color: #404040; }
	td.helpBod	{ border-bottom: 1px solid #9CF;border-top: 0px;border-left: 1px solid #9CF;border-right: 0px;text-align: left;
	text-indent: 10px;font-family: Verdana, sans-serif, Arial;font-weight: normal;font-size: 11px;color: #404040;
	background-color: #fafafa; }
	table.sofT	{ text-align: center;font-family: Verdana;font-weight: normal;font-size: 11px;color: #404040;width: 580px;
	background-color: #fafafa;border: 1px #6699CC solid;border-collapse: collapse;border-spacing: 0px; }"; 
	print OUT "</style>";
	print OUT "</head>";
	print OUT "<body>";
	print DESC "<html>"; 
	print DESC "<head>"; 
	print OUT "<title>${tool} Column Description for Statistics per sample</title>"; 
	print DESC "<style>";
	print DESC "table.helpT	{ text-align: center;font-family: Verdana;font-weight: normal;font-size: 11px;color: #404040;width: auto;
	background-color: #fafafa;border: 1px #6699CC solid;border-collapse: collapse;border-spacing: 0px; }
	td.helpHed	{ border-bottom: 2px solid #6699CC;border-left: 1px solid #6699CC;background-color: #BEC8D1;text-align: left;
	text-indent: 5px;font-family: Verdana;font-weight: bold;font-size: 11px;color: #404040; }
	td.helpBod	{ border-bottom: 1px solid #9CF;border-top: 0px;border-left: 1px solid #9CF;border-right: 0px;text-align: left;
	text-indent: 10px;font-family: Verdana, sans-serif, Arial;font-weight: normal;font-size: 11px;color: #404040;
	background-color: #fafafa; }
	table.sofT	{ text-align: center;font-family: Verdana;font-weight: normal;font-size: 11px;color: #404040;width: 580px;
	background-color: #fafafa;border: 1px #6699CC solid;border-collapse: collapse;border-spacing: 0px; }"; 
	print DESC "</style>";
	print DESC "</head>";
	print DESC "<body>";
	
	##########################################################
	
	print OUT "<p align='center'> §§§ <b>Mayo BIC PI Support</b> §§§   </p>";
	# making the index for the document
	print OUT "<a name=\"top\"></a>";
	print OUT "
	<table id=\"toc\" class=\"toc\" ><tr><td><div id=\"toctitle\"><h2>Contents</h2></div>
	<ul>
	<li class=\"toclevel-1\"><a href=\"#Project Title\"><span class=\"tocnumber\">1</span>
	<span class=\"toctext\">Project Title</span></a></li>
	<li class=\"toclevel-1\"><a href=\"#Project Description\"><span class=\"tocnumber\">2</span>
	<span class=\"toctext\">Project Description</span></a></li>
	<ul>
	<li class=\"toclevel-2\"><a href=\"#Background\"><span class=\"tocnumber\">2.1</span>
	<span class=\"toctext\">Background</span></a></li>
	<li class=\"toclevel-2\"><a href=\"#Study design\"><span class=\"tocnumber\">2.2</span>
	<span class=\"toctext\">Study design</span></a></li>
	</ul>
	<li class=\"toclevel-1\"><a href=\"#Analysis plan\"><span class=\"tocnumber\">3</span>
	<span class=\"toctext\">Analysis plan</span></a></li>
	<li class=\"toclevel-1\"><a href=\"#Received Data\"><span class=\"tocnumber\">4</span>
	<span class=\"toctext\">Received Data</span></a></li>
	<ul>
	<li class=\"toclevel-2\"><a href=\"#Sample Summary\"><span class=\"tocnumber\">4.1</span>
	<span class=\"toctext\">Sample Summary</span></a></li>
	</ul>
	<li class=\"toclevel-1\"><a href=\"#Results Summary\"><span class=\"tocnumber\">5</span>
	<span class=\"toctext\">Results Summary</span></a></li>
	<ul>
	<li class=\"toclevel-2\"><a href=\"#Statistics based on per sample analysis\"><span class=\"tocnumber\">5.1</span>
	<span class=\"toctext\">Statistics based on per sample analysis</span></a></li>";
	if ( ( $analysis eq 'variant' ) || ( $analysis eq 'external' ) || ( $analysis eq 'mayo' ) || ($analysis eq 'realignment') )
	{
		if ($tool eq 'exome' )	{
			print OUT "
			<li class=\"toclevel-2\"><a href=\"#Percent coverage of CaptureRegion\"><span class=\"tocnumber\">5.2</span>
			<span class=\"toctext\">Percent coverage of CaptureRegion</span></a></li>";
		}
		elsif ($tool eq 'whole_genome')	{
			print OUT "
			<li class=\"toclevel-2\"><a href=\"#Percent coverage of CodingRegion\"><span class=\"tocnumber\">5.2</span>
			<span class=\"toctext\">Percent coverage of CodingRegion</span></a></li>";
		}	
	}
	print OUT "
	</ul>
	<li class=\"toclevel-1\"><a href=\"#Results and Conclusions\"><span class=\"tocnumber\">6</span>
	<span class=\"toctext\">Results and Conclusions</span></a></li>
	<li class=\"toclevel-1\"><a href=\"#Results Delivered\"><span class=\"tocnumber\">7</span>
	<span class=\"toctext\">Results Delivered</span></a></li>
	<li class=\"toclevel-1\"><a href=\"#Useful Links\"><span class=\"tocnumber\">8</span>
	<span class=\"toctext\">Useful Links</span></a></li></ul><br>
	</td></tr></table>
	</script>";
	print OUT "<a name=\"Project Title\" id=\"Project Title\"></a><p align='left'><b><u> I. Project Title : </p></b></u>\n";
	print OUT "<ul><table cellspacing=\"0\" class=\"sofT\" > <tr> <td class=\"helpHed\">NGS Bioinformatics for ${tool} sequencing</td> </tr> </table> <br></ul>\n";
	my $read_call;	if($paired == 1)	{	$read_call = 'PE';	}	else	{	$read_call = 'SR';	}	
	my $num_samples=scalar(@sampleArray);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;   ## to pull todays date
	$year += 1900;$mon++;
	print OUT "<a name=\"Project Description\" id=\"Project Description\"></a><p align='left'><b><u> II. Project Description</p></b></u>";
	print OUT "<ul><a name=\"Background\" id=\"Background\"></a><p align='left'> 1. Background</p>";
	print OUT "<table cellspacing=\"0\" class=\"sofT\"><tr>
	<td class=\"helpHed\">Item</td>
	<td class=\"helpHed\">Description</td></tr>
	<td class=\"helpBod\">Disease Type</td><td class=\"helpBod\">$disease</td></tr>
	<td class=\"helpBod\">Number of Samples</td><td class=\"helpBod\">$num_samples</td></tr>";
	if ($analysis ne 'annotation')	{
		print OUT "<td class=\"helpBod\">Read Length</td><td class=\"helpBod\">$read_length</td></tr>
		<td class=\"helpBod\">PairedEnd(PE)/SingleRead(SR)</td><td class=\"helpBod\">$read_call</td></tr>";
	}
	print OUT "<td class=\"helpBod\">Genome Build (hg18/hg19)</td><td class=\"helpBod\">$GenomeBuild</td></tr>
	<td class=\"helpBod\">StartDate</td><td class=\"helpBod\">$date</td></tr>
	<td class=\"helpBod\">EndDate</td><td class=\"helpBod\">$mon/$mday/$year</td></tr>
	</table>";
	print OUT "Note: Further raw NGS data will be used for statistical analysis<br>\n";
	my $loc=$path;
	$loc =~ s/\//\\/g;
	my @WG_ver=split(/\//,$script_path);
	print OUT "Location:: <b><u>\\\\rcfcluster-cifs$loc</b></u> <br>";
	print OUT "(Data is available for 60 Days from the Delivered Date)<br>";
	print OUT "<a name=\"Study design\" id=\"Study design\"></a><p align='left'> 2. Study design</p>";
	print OUT "<ul>
	<li><b> What are the samples? </b><br>
	${sampleinfo}
	<li><b> Goals of the project</b><br>";
	if ($analysis eq 'alignment')	{
		print OUT "Aligning sequencing samples using ${Aligner}";
	}
	else	{	
		if ($tool eq 'whole_genome')	{
			print OUT "${tool} sequencing analysis on the given set of samples, to align, realign, recalibrate, identify SNVs, INDELs, CNVs and SVs and annotate them.";
		}
		elsif ($tool eq 'exome')	{
			print OUT "${tool} sequencing analysis on the given set of samples, to align, realign, recalibrate, identify SNVs, INDELs and annotate them";
		}	
	}	
	print OUT "</ul></ul>\n";
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	print OUT "<a name=\"Analysis Plan\" id=\"Analysis Plan\"></a><p align='left'><b><u> III. Analysis Plan</p></b></u><br>\n";
	print OUT "<b><P ALIGN=\"CENTER\">What's new: New features and updates </b> <a href= \"http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/TREAT\"target=\"_blank\">Genome_GPS</a></P>";
	print OUT "<a href= \"${tool}_workflow.png\"target=\"_blank\"><P ALIGN=\"CENTER\"><img border=\"0\" src=\"${tool}_workflow.png\" width=\"700\" height=\"478\"><b><u><caption align=\"bottom\">Genome_GPS_0.1</caption></b></u></p>";
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	print OUT "<a name=\"Received Data\" id=\"Received Data\"></a><p align='left'><b><u> IV. Received Data</p></b></u> 
	<ul>
	1. Run Name<br>
	<br><table cellspacing=\"0\" class=\"sofT\"><tr><td class=\"helpHed\">Run #</td><td class=\"helpBod\">$run_num</td> </table><br>\n";
	# printing the samples information table
	
	print OUT "<a name=\"Sample Summary\" id=\"Sample Summary\"></a> 2. Sample Summary<br>
        <br><table cellspacing=\"0\" class=\"sofT\"><tr>
        <td class=\"helpHed\">Lane</td>
        <td class=\"helpHed\">Index</td>
        <td class=\"helpHed\">Sample Names</td></tr>";
        for (my $i = 0; $i < $num_samples; $i++)        {
                my @LaneNumbers= split(/,/,$laneArray[$i]);
                my @Indexes=split(/,/,$IndexArray[$i]);
                my $len_lanes = scalar(@LaneNumbers);
                for (my $j =0; $j < $len_lanes; $j++)   {
                        print OUT "
                        <td class=\"helpBod\">$LaneNumbers[$j]</td>
                        <td class=\"helpBod\">$Indexes[$j]</td>
                        <td class=\"helpBod\">$sampleArray[$i]</td></tr>
                        \n";
                }
        }
        print OUT "<td class=\"helpHed\"></td><td class=\"helpHed\"></td><td class=\"helpHed\"></td></tr>";
        print OUT "</table>";
        print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";

	if ( $analysis eq 'mayo' )	{
		print OUT "<b><u>NOTE:</u></b>Meta Data Information available for samples in the form of Lab Tracking Report(<u><a href=\"LTR.xls\"target=\"_blank\">LTR</a></u>)<br>";
	}	
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	print OUT "</ul>
	<a name=\"Results Summary\" id=\"Results Summary\"></a><p align='left'><u><b> V.  Results Summary:</p></u></b>\n";
	if ( ($analysis eq "mayo") || ($analysis eq "external" ) || ($analysis eq "alignment") )	{
		if($fastqc eq "YES")	{
			print OUT "
			<ul>
			<li><a name=\"QC steps\" id=\"QC steps\"></a> QC steps - FastQC-report
			<ul>
			FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.
			<u> <a href= \"fastqc\"target=\"_blank\">ClickMe</a></u>
			<br><b><u>NOTE:</b></u> FastQC runs a series of tests and will flag potential problems with your data<br>\n";
			print OUT "</ul></ul>";
		}
		else	{
			print OUT "
			<ul>
			<li><a name=\"QC steps\" id=\"QC steps\"></a> QC steps - FastQC-report
			<ul>
			FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.
			<u> <a href= \"http://${server}/reports/$run_num/fastqc/ \"target=\"_blank\"><br>ClickMe</a></u>
			<br><b><u>NOTE:</b></u> FastQC runs a series of tests and will flag up and potential problems with your data<br>\n";
		print OUT "</ul></ul>";
		}	
	}
	print OUT"
	<ul><li><a name=\"Statistics based on per sample analysis\" id=\"Statistics based on per sample analysis\"></a>Statistics based on per Sample Analysis(<u><a href=\"StatiticsDesciption.html\"target=\"_blank\">ColumnDescription</a></u>)<br>\n";
	print OUT "<br><table cellspacing=\"0\" class=\"sofT\"><tr><td class=\"helpHed\"><p align='center'></td>";
	my %sample_numbers=();
	my $pairs=();
	my $uniq;
	# storing all the numbers in a Hash per sample (opne hash)
	for(my $k = 0; $k < $num_samples;$k++)	
	{
		print OUT "<td class=\"helpHed\"><p align='center'>$sampleArray[$k]</td>";
		my $file="$path/numbers/$sampleArray[$k].out";
		open SAMPLE, "<$file", or die "could not open $file : $!";
		print"reading numbers from $sampleArray[$k]\n";
		my $id=0;
		while(my $l = <SAMPLE>)	
		{			
			chomp $l;
			if ( $l !~ /^\d/)	{
				$uniq = $id;
				$id++;
			}	
			else	{
				push (@{$sample_numbers{$uniq}},$l);
			}
		}	
		close SAMPLE;
	}
	
	if ($multi eq 'YES')	{
		
	}	
	
	print OUT "</tr>";
	my @To_find;
	my ( $avg_per, $avg_mapped_reads, $per_mapped, $per_mapped_reads );
	my @what;
	# header description
	if ( ( $analysis ne 'annotation' ) && ( $analysis ne 'alignment' ) ) {		
		if ($tool eq 'whole_genome')	{
			if ($analysis eq 'variant')	{
				@what=("Total number of reads obtained","Number of reads mapped to the reference using ${Aligner}","Number of mapped reads overlapping with coding regions in UCSC RefFlat","Total number of SNVs obtained using ${SNV_caller}","Filtered SNVs obtained using ${SNV_caller} recommendations","Number of SNVs observed in  UCSC RefFlat coding regions","Total number of SNVs in dbSNP$dbsnp_v or 1000 Genomes", "Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Number of SeattleSeq annotated Nonsense mutations","Number of SeattleSeq annotated Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total number of  SNVs not in dbSNP$dbsnp_v or 1000 Genomes", "Novel Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Novel number of Nonsense mutations","Novel number of Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total  number of INDELs obtained using ${SNV_caller}","Filtered INDELs obtained using ${SNV_caller} recommendations","Number of INDELs observed in  UCSC RefFlat coding region", "Number of SeattleSeq annotated INDELs in coding regions with inserted or deleted bases being multiples of 3","Number of SeattleSeq annotated Frameshifft INDELs with inserted or deleted bases which are not multiples of 3","Number of SeattleSeq annotated Splice-3 INDELs","Number of SeattleSeq annotated Splice-5 INDELs","Number of SeattleSeq annotated UTR-3 INDELs","Number of SeattleSeq annotated UTR-5 INDELs","Total number of CNVs","Number of CNVs observed in UCSC RefFlat coding regions","Number of deletions in UCSC RefFlat coding regions","Number of duplications in UCSC RefFlat coding regions","total number of SVs","Number of SVs observed in UCSC RefFlat coding regions","Number of ITX in UCSC RefFlat coding regions","Number of INV in UCSC RefFlat coding regions","Number of DEL in UCSC RefFlat coding regions","Number of INS in UCSC RefFlat coding regions","Number of CTX in UCSC RefFlat coding regions");
			}
			else	{
				@what=("Total number of reads obtained","Number of reads mapped to the reference using ${Aligner}","Percent of duplicated reads flagged in BAM","Number of reads after recalibration and realignment using GATK ","Number of mapped reads overlapping with coding regions in UCSC RefFlat","Total number of SNVs obtained using ${SNV_caller}","Filtered SNVs obtained using ${SNV_caller} recommendations","Number of SNVs observed in  UCSC RefFlat coding regions","Total number of SNVs in dbSNP$dbsnp_v or 1000 Genomes", "Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Number of SeattleSeq annotated Nonsense mutations","Number of SeattleSeq annotated Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total number of  SNVs not in dbSNP$dbsnp_v or 1000 Genomes", "Novel Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Novel number of Nonsense mutations","Novel number of Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total  number of INDELs obtained using ${SNV_caller}","Filtered INDELs obtained using ${SNV_caller} recommendations","Number of INDELs observed in  UCSC RefFlat coding region", "Number of SeattleSeq annotated INDELs in coding regions with inserted or deleted bases being multiples of 3","Number of SeattleSeq annotated Frameshifft INDELs with inserted or deleted bases which are not multiples of 3","Number of SeattleSeq annotated Splice-3 INDELs","Number of SeattleSeq annotated Splice-5 INDELs","Number of SeattleSeq annotated UTR-3 INDELs","Number of SeattleSeq annotated UTR-5 INDELs","Total number of CNVs","Number of CNVs observed in UCSC RefFlat coding regions","Number of deletions in UCSC RefFlat coding regions","Number of duplications in UCSC RefFlat coding regions","total number of SVs","Number of SVs observed in UCSC RefFlat coding regions","Number of ITX in UCSC RefFlat coding regions","Number of INV in UCSC RefFlat coding regions","Number of DEL in UCSC RefFlat coding regions","Number of INS in UCSC RefFlat coding regions","Number of CTX in UCSC RefFlat coding regions");
			}
		}
		elsif ($tool eq 'exome')	{
			if ($analysis eq 'variant')	{
				@what=("Total number of reads obtained","Number of reads mapped to the reference using ${Aligner}","Number of mapped reads overlapping with coding regions in UCSC RefFlat","Total number of SNVs obtained using ${SNV_caller}","Filtered SNVs obtained using ${SNV_caller} recommendations","Number of SNVs observed in  UCSC RefFlat coding regions","Number of SNVs observed in  Capture Region","Total number of SNVs in dbSNP$dbsnp_v or 1000 Genomes", "Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Number of SeattleSeq annotated Nonsense mutations","Number of SeattleSeq annotated Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total number of  SNVs not in dbSNP$dbsnp_v or 1000 Genomes", "Novel Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Novel number of Nonsense mutations","Novel number of Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total  number of INDELs obtained using ${SNV_caller}","Filtered INDELs obtained using ${SNV_caller} recommendations","Number of INDELs observed in  UCSC RefFlat coding region","Number of INDELs observed in  Capture Region", "Number of SeattleSeq annotated INDELs in coding regions with inserted or deleted bases being multiples of 3","Number of SeattleSeq annotated Frameshifft INDELs with inserted or deleted bases which are not multiples of 3","Number of SeattleSeq annotated Splice-3 INDELs","Number of SeattleSeq annotated Splice-5 INDELs","Number of SeattleSeq annotated UTR-3 INDELs","Number of SeattleSeq annotated UTR-5 INDELs","Total number of CNVs","Number of CNVs observed in UCSC RefFlat coding regions","Number of deletions in UCSC RefFlat coding regions","Number of duplications in UCSC RefFlat coding regions","total number of SVs","Number of SVs observed in UCSC RefFlat coding regions","Number of ITX in UCSC RefFlat coding regions","Number of INV in UCSC RefFlat coding regions","Number of DEL in UCSC RefFlat coding regions","Number of INS in UCSC RefFlat coding regions","Number of CTX in UCSC RefFlat coding regions");
			}
			else	{
				@what=("Total number of reads obtained","Number of reads mapped to the reference using ${Aligner}","Percent of duplicated reads flagged in BAM","Number of reads after recalibration and realignment using GATK ","Number of mapped reads overlapping with coding regions in UCSC RefFlat","Total number of SNVs obtained using ${SNV_caller}","Filtered SNVs obtained using ${SNV_caller} recommendations","Number of SNVs observed in  UCSC RefFlat coding regions","Number of SNVs observed in  Capture Region","Total number of SNVs in dbSNP$dbsnp_v or 1000 Genomes", "Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Number of SeattleSeq annotated Nonsense mutations","Number of SeattleSeq annotated Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total number of  SNVs not in dbSNP$dbsnp_v or 1000 Genomes", "Novel Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Novel number of Nonsense mutations","Novel number of Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total  number of INDELs obtained using ${SNV_caller}","Filtered INDELs obtained using ${SNV_caller} recommendations","Number of INDELs observed in  UCSC RefFlat coding region","Number of INDELs observed in  Capture Region", "Number of SeattleSeq annotated INDELs in coding regions with inserted or deleted bases being multiples of 3","Number of SeattleSeq annotated Frameshifft INDELs with inserted or deleted bases which are not multiples of 3","Number of SeattleSeq annotated Splice-3 INDELs","Number of SeattleSeq annotated Splice-5 INDELs","Number of SeattleSeq annotated UTR-3 INDELs","Number of SeattleSeq annotated UTR-5 INDELs","Total number of CNVs","Number of CNVs observed in UCSC RefFlat coding regions","Number of deletions in UCSC RefFlat coding regions","Number of duplications in UCSC RefFlat coding regions","total number of SVs","Number of SVs observed in UCSC RefFlat coding regions","Number of ITX in UCSC RefFlat coding regions","Number of INV in UCSC RefFlat coding regions","Number of DEL in UCSC RefFlat coding regions","Number of INS in UCSC RefFlat coding regions","Number of CTX in UCSC RefFlat coding regions");
			}
		}
	}
	
	# header name
	if ( ( $analysis ne 'alignment' ) && ( $analysis ne 'annotation' ) ) {
		if ($tool eq 'whole_genome')	{
			if ($analysis eq 'variant')	{
				@To_find=("Total Reads","Mapped Reads","Mapped Reads(CodingRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5","Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
			}
			else	{
				@To_find=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(CodingRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5","Total CNVs","Coding CNVs","Coding Deletions","Coding Duplications","Total SVs","Coding SVs","Intra-chr translocations","Inversions","Deletions","Insertions","Inter-chr translocations");
			}
		}
		elsif ($tool eq 'exome')	{
			if ($analysis eq 'variant')	{
				@To_find=("Total Reads","Mapped Reads","Percent duplication","Realigned Mapped Reads","Mapped Reads(in CaptureRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","INDELs in CaptureRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
			}
			else	{
				@To_find=("Total Reads","Mapped Reads","Mapped Reads(in CaptureRegion)","Total SNVs (${SNV_caller})","Filtered SNVs (${SNV_caller})","SNVs in CodingRegion","SNVs in CaptureRegion","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs (${SNV_caller})","Filtered INDELs (${SNV_caller})","INDELs in CodingRegion","INDELs in CaptureRegion","In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
			}
		}	
	}
	
	if ($analysis eq 'alignment')	{
		@To_find=("Total Reads","Mapped Reads","Percent duplication");
		@what=("Total number of reads obtained","Number of reads mapped to the reference using ${Aligner}","Percent of duplicated reads flagged in BAM");
	}	
	
	if ($analysis eq 'annotation')	{
		@what=("Total SNVs","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total SNVs","Transition to Trasnversion Ratio","Nonsense","Missense","Coding-synonymous","Coding-notMod3","Splice-3","Splice-5","UTR-3","UTR-5","Total INDELs" ,"In Coding","Leading to Frameshift","Splice-3","Splice-5","UTR-3","UTR-5");
		@To_find=("Total number of SNVs obtained using","Total number of SNVs in dbSNP$dbsnp_v or 1000 Genomes", "Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Number of SeattleSeq annotated Nonsense mutations","Number of SeattleSeq annotated Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total number of  SNVs not in dbSNP$dbsnp_v or 1000 Genomes", "Novel Transition to Transversion Ratio. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)","Novel number of Nonsense mutations","Novel number of Missense mutations", "Number of SeattleSeq annotated coding-synonymous mutations","Number of SeattleSeq annotated coding not-Mod3 mutations","Number of SeattleSeq annotated splice-3 mutations","Number of SeattleSeq annotated splice-5 mutations","Number of SeattleSeq annotated utr-3 mutations","Number of SeattleSeq annotated utr-5 mutations","Total  number of INDELs obtained using" , "Number of SeattleSeq annotated INDELs in coding regions with inserted or deleted bases being multiples of 3","Number of SeattleSeq annotated Frameshifft INDELs with inserted or deleted bases which are not multiples of 3","Number of SeattleSeq annotated Splice-3 INDELs","Number of SeattleSeq annotated Splice-5 INDELs","Number of SeattleSeq annotated UTR-3 INDELs","Number of SeattleSeq annotated UTR-5 INDELs");
	}	
	
	foreach my $key (sort {$a <=> $b} keys %sample_numbers)	{
		print OUT "<td class=\"helpHed\"><p align='left'><a href=\"#$To_find[$key]\" title=\"$what[$key]\">$To_find[$key]</a></td>";
		if ( $key eq '1')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c] / ${$sample_numbers{0}}[$c]) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\">$print <br> <b>($per_mapped \%) <b></td>";
					$avg_mapped_reads=$avg_mapped_reads + ${$sample_numbers{$key}}[$c];
					$per_mapped_reads = $per_mapped_reads + $per_mapped; 
				}
				$avg_mapped_reads=$avg_mapped_reads/$num_samples;$avg_mapped_reads=int($avg_mapped_reads/1000000);
				$per_mapped_reads =$per_mapped_reads/$num_samples;
				print OUT "</tr>\n";
		}		
		if ($tool eq 'whole_genome')	{
			if ( $key eq '39')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_deleted = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / 3200000000) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\">$print <br> <b>($per_deleted \%) <b></td>";
				}
			}	
			if ( $key eq '40')	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_duplicated = sprintf("%.5f",(${$sample_numbers{$key}}[$c] / 3200000000) * 100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\"> $print<br><b>($per_duplicated \%) <b></td>";
				}
			}	
		}
		if ($key eq '2')	{
			for (my $c=0; $c < $num_samples;$c++)	{
				my $print=sprintf("%.2f",$sample_numbers{$key}[$c]);
				print OUT "<td class=\"helpBod\"> $print\%</td>";
			}	
			print OUT "</tr>\n";
		}	
		if ( ( $key eq '3' ) )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\">$print <br><b>($per_mapped \%)<b></td>";
				}
				print OUT "</tr>\n";
			}
		}	
		if ( ( $key eq '4' ) )	{
			if ( $analysis ne 'annotation' )	{
				for (my $c=0; $c < $num_samples;$c++)	{
					my $per_mapped_t = sprintf("%.1f",(${$sample_numbers{$key}}[$c]/${$sample_numbers{0}}[$c])*100);
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\">$print <br><b>($per_mapped_t \%)<b></td>";
				}
				print OUT "</tr>\n";
			}
		}		
		if ($analysis ne 'annotation' )	{
			for ( my $c=0; $c < $num_samples; $c++ )	{
				if ( ( $key eq '1') || ( $key eq '3' ) || ($key eq '2') || ($key eq '4') || ( $key eq '39' ) || ( $key eq '40' ) )	{
				}
				else	{
					my $print=CommaFormatted(${$sample_numbers{$key}}[$c]);
					print OUT "<td class=\"helpBod\">$print</td>\n";	
				}
			}
		}
		print OUT "</tr>";
		# key 1 is for mapped reads
		# key 2 is for mapped reads ontarget
		if ( ( $analysis eq 'external' ) || ( $analysis eq 'variant' ) || ($analysis eq 'realignment') || ($analysis eq 'mayo') ) {	
			if ($key eq '4' )	{
				print OUT "<th class=\"helpBod\">Single Nucleotide Variants (SNVs)</th>";
				for (my $c=0; $c < $num_samples;$c++)	{
					print OUT "<td class=\"helpHed\"></td>";
				}		
				print OUT "</tr>\n";
			}	
			if ($key eq '7' )	{
				print OUT "<td class=\"helpBod\"><b>Known SNVs<b></td>";
				for (my $c=0; $c < $num_samples;$c++)	{
					print OUT "<td class=\"helpHed\"></td>";
				}		
				print OUT "</tr>\n";
			}	
			if ($key eq '17' )	{
				print OUT "<td class=\"helpBod\"><b>Novel SNVs<b></td>";
				for (my $c=0; $c < $num_samples;$c++)	{
					print OUT "<td class=\"helpHed\"></td>";
				}		
				print OUT "</tr>\n";
			}	

			if ($key eq '27' )	{
				print OUT "<th class=\"helpBod\">INsertions DELetions (INDELs)</th>";
				for (my $c=0; $c < $num_samples;$c++)	{
					print OUT "<td class=\"helpHed\"></td>";
				}		
				print OUT "</tr>\n";
			}
			if ($tool eq 'whole_genome')	{
				if ($key eq '36' )	{
					print OUT "<th class=\"helpBod\">Copy Number Variants (CNVs)</th>";
					for (my $c=0; $c < $num_samples;$c++)	{
						print OUT "<td class=\"helpHed\"></td>";
					}		
					print OUT "</tr>\n";
				}
				if ($key eq '40' )	{
					print OUT "<th class=\"helpBod\">Structural Variants (SVs)</th>";
					for (my $c=0; $c < $num_samples;$c++)	{
						print OUT "<td class=\"helpHed\"></td>";
					}		
					print OUT "</tr>\n";
				}
			}
		}
	# print OUT "</tr>";
	}	
	undef %sample_numbers;
	
	print OUT "</table>";
	print DESC "
	<p align='left'><b> Row description for Statistics Table:</p></b>
	<table cellspacing=\"0\" class=\"sofT\"><tr><td class=\"helpHed\">Column</td><td class=\"helpHed\">Description</td></tr>";
	if ($analysis eq 'alignment' )	{
		print DESC"
		<td class=\"helpBod\">Total Reads</td><td class=\"helpBod\">Total number of reads obtained</td></tr>
		<td class=\"helpBod\">Mapped Reads</td><td class=\"helpBod\">Number and percentage of reads mapped to reference genome(${GenomeBuild})</td></tr>
		<td class=\"helpBod\">Percent duplication</td><td class=\"helpBod\">Percentage of duplicated reads as identified by(${SNV_caller})</td></tr>
		<td class=\"helpBod\">Realigned Reads (Q >= 20)</td><td class=\"helpBod\">Number and percenatge of realigned reads with a quality >= 20 mapped to reference genome(${GenomeBuild})</td></tr>
		</table></ul>";	
	}	
	
	if ( ($analysis eq "external") || ($analysis eq "variant" ) || ($analysis eq 'realignment') || ($analysis eq "mayo")  )	{
		print DESC"
		<td class=\"helpBod\">Total Reads</td><td class=\"helpBod\">Total number of reads obtained</td></tr>
		<td class=\"helpBod\">Mapped Reads</td><td class=\"helpBod\">Number and percentage of reads mapped to reference genome(${GenomeBuild})</td></tr>
		<td class=\"helpBod\">Percent duplication</td><td class=\"helpBod\">Percentage of duplicated reads as identified by(${SNV_caller}) and flagged in the BAM files</td></tr>
		<td class=\"helpBod\">Realigned Reads (Q >= 20)</td><td class=\"helpBod\">Number and percenatge of realigned reads with a quality >= 20 mapped to reference genome(${GenomeBuild})</td></tr>
		<td class=\"helpBod\">Total SNVs ${SNV_caller}</td><td class=\"helpBod\">Total number of SNVs obtained using ${SNV_caller} </td></tr>
		<td class=\"helpBod\">Filtered SNVs ${SNV_caller}</td><td class=\"helpBod\">Number of SNVs obtained after passing through ${SNV_caller} VQSR filtering</td></tr>
		<td class=\"helpBod\">Coding SNVs</td><td class=\"helpBod\">Number of SNVs observed to occur in UCSC RefFlat genomic coordinates</td></tr>
		<td class=\"helpBod\">Total SNVs </td><td class=\"helpBod\">Total number of known SNVs, either in dbSNP$dbsnp_v or 1000 Genomes</td></tr>
		<td class=\"helpBod\">Transition to Transversion ratio</td><td class=\"helpBod\">Number of transtitions over number of transversions. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)</td></tr>
		<td class=\"helpBod\">Nonsense</td><td class=\"helpBod\">Number of Known SNVs that lead to a stop codon</td></tr>
		<td class=\"helpBod\">Missense</td><td class=\"helpBod\">Number of Known SNVs that lead to codons coding for different amino acids</td></tr>
		<td class=\"helpBod\">Coding-synonymous</td><td class=\"helpBod\">Number of Known SNVs that lead to codon change without changing the amino acid</td></tr>
		<td class=\"helpBod\">Coding-notMod3</td><td class=\"helpBod\">Number of Known SNVs that lead to codon change, but number of coding bases is not a multiple of 3</td></tr>
		<td class=\"helpBod\">Splice-3</td><td class=\"helpBod\">Number of Known SNVs present in the splice-3 sites</td></tr>
		<td class=\"helpBod\">Splice-5</td><td class=\"helpBod\">Number of Known SNVs present in the splice-5 sites</td></tr>
		<td class=\"helpBod\">UTR-3</td><td class=\"helpBod\">Number of Known SNVs present in the utr-3 sites</td></tr>
		<td class=\"helpBod\">UTR-5</td><td class=\"helpBod\">Number of Known SNVs present in the utr-5 sites</td></tr>
		<td class=\"helpBod\">Total SNVs </td><td class=\"helpBod\">Total number of Novel SNVs, either in dbSNP or 1000 Genomes</td></tr>
		<td class=\"helpBod\">Transition to Transversion ratio</td><td class=\"helpBod\">Number of transtitions over number of transversions. Transition is defined as a change among purines or pyrimidines (A to G, G to A, C to T, and T to C). Transversion is defined as a change from purine to pyrimidine or vice versa (A to C, A to T, G to C, G to T, C to A, C to G, T to A, and T to G)</td></tr>
		<td class=\"helpBod\">Nonsense</td><td class=\"helpBod\">Number of Novel SNVs that lead to a stop codon</td></tr>
		<td class=\"helpBod\">Missense</td><td class=\"helpBod\">Number of Novel SNVs that lead to codons coding for different amino acids</td></tr>
		<td class=\"helpBod\">Coding-synonymous</td><td class=\"helpBod\">Number of Novel SNVs that lead to codon change without changing the amino acid</td></tr>
		<td class=\"helpBod\">Coding-notMod3</td><td class=\"helpBod\">Number of Novel SNVs that lead to codon change, but number of coding bases is not a multiple of 3</td></tr>
		<td class=\"helpBod\">Splice-3</td><td class=\"helpBod\">Number of Novel SNVs present in the splice-3 sites</td></tr>
		<td class=\"helpBod\">Splice-5</td><td class=\"helpBod\">Number of Novel SNVs present in the splice-5 sites</td></tr>
		<td class=\"helpBod\">UTR-3</td><td class=\"helpBod\">Number of Novel SNVs present in the utr-3 sites</td></tr>
		<td class=\"helpBod\">UTR-5</td><td class=\"helpBod\">Number of Novel SNVs present in the utr-5 sites</td></tr>
		<td class=\"helpBod\">Total INDELs ${SNV_caller}</td><td class=\"helpBod\">Total number of INDELs identified by ${SNV_caller} </td></tr>
		<td class=\"helpBod\">Filtered INDELs ${SNV_caller}</td><td class=\"helpBod\">Number of INDELs obtained after passing through the VQSR filtering by ${SNV_caller}</td></tr>
		<td class=\"helpBod\">Coding INDELs</td><td class=\"helpBod\">Number of INDELs observed to occur in genomic coordinates</td></tr>
		<td class=\"helpBod\">In coding</td><td class=\"helpBod\">Number of SeattleSeq annotated INDELs present in coding regions of the genome where inserted or deleted bases are multiples of 3</td></tr>
		<td class=\"helpBod\">Leading to Frameshift</td><td class=\"helpBod\">Number of INDELs leading to frameshift mutations where inserted or deleted bases ae not multiples of 3</td></tr>
		<td class=\"helpBod\">Splice-3</td><td class=\"helpBod\">Number of INDELs present in the splice-3 sites</td></tr>
		<td class=\"helpBod\">Splice-5</td><td class=\"helpBod\">Number of INDELs present in the splice-5 sites</td></tr>
		<td class=\"helpBod\">UTR-3</td><td class=\"helpBod\">Number of INDELs present in the utr-3 sites</td></tr>
		<td class=\"helpBod\">UTR-5</td><td class=\"helpBod\">Number of INDELs present in the utr-5 sites</td></tr>";
		
		if ($tool eq 'whole_genome')	{
			print DESC"
			<td class=\"helpBod\">Total CNVs</td><td class=\"helpBod\">Total number of CNVs identified</td></tr>
			<td class=\"helpBod\">Coding CNVs</td><td class=\"helpBod\">Number of CNVs observed to occur in genomic coordinates</td></tr>
			<td class=\"helpBod\">Deletions</td><td class=\"helpBod\">Number of deletions found in coding regions</td></tr>
			<td class=\"helpBod\">Duplications</td><td class=\"helpBod\">Number of duplications found in coding regions</td></tr>
			<td class=\"helpBod\">Total SVs</td><td class=\"helpBod\">Ttoal number of SVs identified</td></tr>
			<td class=\"helpBod\">Coding SVs</td><td class=\"helpBod\">Number of SVs observed to occur in genomic coordinates</td></tr>
			<td class=\"helpBod\">Intra-chr translocations</td><td class=\"helpBod\">Number of intra-chromosomal translocations</td></tr>
			<td class=\"helpBod\">Inversions</td><td class=\"helpBod\">Number of inversions</td></tr>
			<td class=\"helpBod\">Deletions</td><td class=\"helpBod\">Number of deletions</td></tr>
			<td class=\"helpBod\">Insertions</td><td class=\"helpBod\">Number of insertions</td></tr>
			<td class=\"helpBod\">Inter-chr translocations</td><td class=\"helpBod\">Number of inter-chromosomal translocations</td></tr>";
		}	
		print DESC "</table></ul>";
	}
	print DESC "</ul>";
	print DESC "</body>\n"; 
	print DESC "</html>\n"; 
	close DESC;
	print OUT "</ul>";
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	my $re;
	if ($tool eq 'exome')	{
		$re="CaptureRegion";
	}
	elsif($tool eq 'whole_genome')	{
		$re="CodingRegion";
	}	
	if ( ( $analysis eq "external" ) || ( $analysis eq "variant") || ($analysis eq "realignment") || ($analysis eq "mayo")   )	{
		my $target=$target_region/1000000;$target=sprintf("%.2f",$target);
		print OUT "
		<ul><li><a name=\"Percent coverage of $re\" id=\"Percent coverage of $re\"></a>Percent coverage of $re
		<ul>
		The Probes Captured a region of $target Mbp. The figure below lists the percentage of that target region which is covered by at least N depth of coverage.
		</ul></ul>";
		print OUT "<a href= \"Coverage.JPG\"target=\"_blank\"><P ALIGN=\"CENTER\"><img border=\"0\" src=\"Coverage.JPG\" width=\"704\" height=\"628\">";
		print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";		
	}
	print OUT "<a name=\"Results and Conclusions\" id=\"Results and Conclusions\"></a><p align='left'><u><b> VI.  Results and Conclusions:</p></u></b>";
	$per_mapped_reads=sprintf("%.2f",$per_mapped_reads);
	if  ($analysis ne "annotation") 	{
		print OUT"
		<ul>
		<li> $per_mapped_reads % of the data has been mapped to the genome.
		<li> A high throughput with approximately $avg_mapped_reads million reads passed filtering.
		";
	}	
	if ($analysis ne "alignment")	{
		print OUT "</ul><b><u><ul> IGV Visualization</b></u><br>
		The SNV and INDEL annotation reports (both standard and filtered) include visualization links to IGV to enable a realistic view of the variants. Please follow steps in the following link to setup IGV (takes less than 5 minutes) and utilize this feature.<br>
		<a href= \"IGV_Setup.doc\"target=\"_blank\">IGV setup for variant visualization</a><br><br>";
		print OUT "<b><u>Synonymous codon change</b></u><br>
		Genetic code associates a set of sibling codons to the same amino acid, some codons occur more frequently than others in gene sequences. Biased codon usage results from a diversity of factors such as:
		<ul>
		<li>GC-content
		<li>preference for codons with G or C at the third nucleotide position
		<li>leading strand richer in G+T than lagging strand
		<li>translational bias frequently observed in fast growing organisms
		</ul>";
		print OUT "<br><b><u>BGI-Danish200 Exome Dataset</b></u><br>
		The data is from Exome resequencing of 200 individuals of Danish nationality, a collaborative project between BGI and Danish researchers. 
		<ul>
		<li>The BGI200exome column shows the two alleles and their frequencies (minor allele/major allele)
		</ul>";
		print OUT "<br><b><u>COSMIC Dataset</b></u><br>
		COSMIC (Catalogue of Somatic Mutations In Cancer) is a comprehensive source of somatic mutations and associated data found in human cancer. It's developed and maintained by the Cancer Genome Project (CGP) at the Wellcome Trust Sanger Institute. The curated data come from scientific papers in the literature and large scale experimental screens from CGP. 
			<ul>
			<li>The COSMIC column shows the COSMIC mutation ID; Mutation CDS; Mutation AA; strand
		</ul>";
	}	
	print OUT "</ul>";
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	print OUT "<a name=\"Results Delivered\" id=\"Results Delivered\"></a><p align='left'><u><b> VIII. Results Delivered</p></u></b>";
	if (( $analysis eq "external") || ($analysis eq "variant") || ($analysis eq 'realignment') || ($analysis eq "mayo")  )	{
		if ($analysis ne 'alignment')	{
		 print OUT "<ul>
                <li>Merged ${SNV_caller} results along with SIFT and SeattleSeq SNP annotation and Indel Annotation from Seattle Seq for all samples(<u><a href=\"ColumnDescription_Reports.xls\"target=\"_blank\">Column Description for Reports</a></u>)<br>
                <u> <a href= \"Reports/SNV.cleaned_annot.xls\"target=\"_blank\">SNV Report</a></u> <br>
                <u> <a href= \"Reports/INDEL.cleaned_annot.xls\"target=\"_blank\">INDEL Report</a></u> <br><br></ul>";

		print OUT "<ul>
		<li>Per sample SNV and INDEL files are available here comprising of filtered and unfiltered reports. The filtered reports comprise of a much smaller list of variants for investigation and are based on the following criteria:
		<ul>
		- dbSNP$dbsnp_v column does not have an rs ID (novel), OR<br>
		- functionGVS column having 'missense', 'nonsense', 'splice-3', 'splice-5', 'coding-notMod3', 'utr-3' or 'utr-5' (intron, intergenic and coding-synonymous removed using SeattleSeq annotation, OR<br>
		- any variant reported within +/-2bp of an exon edge using 'distance' report for variants
		</ul><br>
		<u> <a href= \"Reports_per_Sample\"target=\"_blank\">Per Sample Reports</a></u></ul>";
		print OUT "<ul>
		<li>Per sample Gene Summary files are available here<br>
		<u> <a href= \"Reports_per_Sample\"target=\"_blank\">Per Sample Reports</a></u></ul>";
		}
	}
	if ( ($analysis eq 'external') || ($analysis eq 'variant' ) || ($analysis eq 'realignment') || ($analysis eq "mayo")  )	{
		print OUT "<ul>
			<li>The variant distance for SNPs and INDELs (the distance to closest Exon) are recorded in two files<br>
			<u> <a href= \"Reports/variantLocation_SNVs\"target=\"_blank\">SNVs</a></u> <br>
            <u> <a href= \"Reports/variantLocation_INDELs\"target=\"_blank\">INDELs</a></u> <br><br></ul>";
	}
	print OUT "</ul>";
	if ($analysis eq 'alignment')	{
		print OUT "<ul>
		<li>Alignment and sorted BAM<br>
		<u> <a href= \"Alignment\"target=\"_blank\">Aligned Bam</a></u></ul>";
	}	
	print OUT "<ul>
		<li>Statistics based on per Sample Analysis are recorded in the tab delimited file<br>
		<u> <a href= \"SampleStatistics.tsv\"target=\"_blank\">SampleStatistics</a></u></ul>";
	
	
	if ($tool eq 'whole_genome' && $analysis ne 'alignment' && $analysis ne 'annotation')	{
		print OUT "<ul> <li>Per sample Circos plots are displayed (and linked) below. Copy Number Variants and Sructural Variants are plotted for each sample. <br>
	<u> <a href= \"circos\"target=\"_blank\">Per sample CIRCOS plots</a></u></ul></ul>";
		for(my $k = 0; $k < $num_samples;$k++)	
		{
			print OUT "<a href= \"circos/${sampleArray[$k]}.sv_cnv.png\"target=\"_blank\"><P ALIGN=\"CENTER\"><img border=\"0\" src=\"circos/${sampleArray[$k]}.sv_cnv.png\" width=\"50\" height=\"50\"><u><caption align=\"bottom\">${sampleArray[$k]}</caption></u></p>";
		}
	}	
	if ($upload_tb eq 'YES')	{
		print OUT "<ul>
		<li>Variant calls can be visualized and filtered using TableBrowser<br>
		<u> <a href= \"http://${host}:${port}/TREATTableBrowser/\"target=\"_blank\">TableBrowser</a></u></ul>";
	}
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
	print OUT "<a name=\"Useful Links\" id=\"Useful Links\"></a><p align='left'><b><u> VII. Useful Links</p></b></u>
	<ul>
	<li><a href= \"http://sift.jcvi.org/www/SIFT_help.html\"target=\"_blank\">SIFT</a>
	<li><a href= \"http://gvs.gs.washington.edu/GVS/HelpSNPSummary.jsp\"target=\"_blank\">SeattleSeq </a>
	<li><a href= \"http://circos.ca/\"target=\"_blank\">Circos </a>
	<li><a href= \"http://www.broadinstitute.org/software/igv/\"target=\"_blank\">Integrative Genomics Viewer</a>
	<li><a href= \"http://www.genecards.org/\"target=\"_blank\">GeneCards</a>
	<li><a href= \"http://www.ncbi.nlm.nih.gov/geo/\"target=\"_blank\">Gene Expression Omnibus</a>
	<li><a href= \"http://genome.ucsc.edu/\"target=\"_blank\">UCSC Genome Browser</a>
	</ul><br>";
	print OUT "
	<b><u>Authorship Consideration</u></b>: Advancing scientific research is a primary motivation of all bioinformaticians and acknowledgment of contribution through authorship on manuscripts arising from this analysis is one way our work is assessed and attributed. We request to be considered for authorship on any manuscripts using the analysis results provided if you believe we have made substantive intellectual contributions to the study.  
	";
	print OUT "
	<p align='center'> §§§ <b>Powered by Mayo BIC PI Support</b> §§§  </p> ";
	print OUT "<p align='right'><a href=\"#top\">-top-</a></p>";
 	print OUT "</body>\n"; 
	print OUT "</html>\n"; 
	close OUT;
	print "Document is generated with path as $output.......... \n";
}
