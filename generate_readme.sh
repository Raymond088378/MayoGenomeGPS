#!/bin/bash

if [ $# != 2 ]
then
	echo "Usage: </path/to/output folder><path/to/run info file>"
else
	echo `date`
	output=$1
	run_info=$2
	
	##Reports Reports_per_Sample 
	multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")
	groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" " ")
	
	##
	echo -e "###File Descriptions\n"
	echo -e "* Session file to visualize alignment BAMs, coverage and UCSC tracks using IGV: igv_session.xml"
	echo -e "* Help guide to set up and visualize in IGV: IGV_Setup.doc"
	echo -e "* Summary report on analysis performed and statistics obtained: Main_Document.html"
	echo -e "* Description of statistics table in Main_Document.html: StatisticsDescription.html"
	if [ $multi == "YES" ]
	then
		echo -e "* Plain text version of statistics table in Main_Document.html for single samples: SampleStatistics.pair.tsv"
	fi
	echo -e "* Plain text version of statistics table in Main_Document.html for single samples: SampleStatistics.tsv"
	echo -e "* Description of columns in SNV and INDEL reports: ColumnDescription_Reports.xls"
	if [  $tool == "whole_genome" ]
	then
		echo -e "* Whole genome workflow diagram: whole_genome_workflow.png"
	else
		echo -e "* Exome workflow diagram: exome_workflow.png"	
	fi	
	if [[ $analysis != "annotation" && $analysis != "ontarget" && $analysis != "alignment" ]]
	then
		echo -e "* Per sample coverage plot: Coverage.JPG"
	fi	
	echo -e "Compressed file with the statistical numbers in GENOME_GPS for all the samples : numbers.tar.gz"
  	
	echo -e "\n* Configuration files for GENOME_GPS:\n"
  	echo -e "-run_info.txt\n-sample_info.txt\n-tool_info.txt"
 	echo -e "\n#### Annotation Reports\n"
	echo -e "Directory with merged SNV and INDEL reports: : Reports/"
	echo -e "**ANNOTATION files (filtered - single annotation line per variant call(most impacting transcript)"
	echo -e "-SNVs\nSNV.xls\nSNV.filtered.xls\n-INDELs:\nINDEL.xls\nINDEL.filtered.xls
	if [ $multi == "YES" ]
	then
		echo -e "-SNVs\nTUMOR.SNV.xls\nTUMOR.SNV.filtered.xls\n-INDELs:\nTUMOR.INDEL.xls\nTUMOR.INDEL.filtered.xls
	fi	
	
	echo -e "Directory with per sample SNV, INDEL reports: Reports_per_Sample/"
	echo -e "-SNVs"
	if [ $multi == "NO" ]
	then
		for sam in $samples
		do
			echo -e " $sam.SNV.xls - There could be multiple annotation lines for each variant for $sam\n$sam.SNV.filtered.xls - just one annotation for each variant for $sam (most impacting one)"
		done
	else
		for pair in $groups
		do
			echo -e " $pair.SNV.xls - There could be multiple annotation lines for each variant for $pair\n$pair.SNV.filtered.xls - just one annotation for each variant for $pair (most impacting one)"
			echo -e " TUMOR.$pair.SNV.xls - There could be multiple annotation lines for each somatic variant for $pair\nTUMOR.$pair.SNV.filtered.xls - just one annotation for each somatic variant for $pair (most impacting one)"
		done
	fi	
	echo -e "-INDELs"
	if [ $multi == "NO" ]
	then
		for sam in $samples
		do
			echo -e " $sam.INDEL.xls - There could be multiple annotation lines for each variant for $sam\n$sam.INDEL.filtered.xls - just one annotation for each variant for $sam (most impacting one)"
		done
	else
		for pair in $groups
		do
			echo -e " $pair.INDEL.xls - There could be multiple annotation lines for each variant for $pair\n$pair.INDEL.filtered.xls - just one annotation for each variant for $pair (most impacting one)"
			echo -e " TUMOR.$pair.INDEL.xls - There could be multiple annotation lines for each somatic variant for $pair\nTUMOR.$pair.INDEL.filtered.xls - just one annotation for each somatic variant for $pair (most impacting one)"
		done
	fi	
	
	echo -e "**VCF files (Variant Calling format files. Could be opened in IGV)"
	if [ $multi == "YES " ]
	then
		for pair in $groups
		do
			echo -e "filtered somatic variant file for $pair : $pair.somatic.variants.filter.vcf"
			echo -e "unfiltered multi allelic somatic variant file for $pair : $pair.somatic.variants.raw.multi.vcf"
			echo -e "filtered multi sample variant file for $pair : $pair.variants.filter.vcf"
			echo -e "unfiltered multi allelic multi sample variant file for $pair : $pair.variants.raw.multi.vcf"
		done
	else
		for sam in $samples
		do			
			echo -e "filtered variant file for $sam : $sam.variants.raw.multi.vcf"
			echo -e "unfiltered multi allelic variant file for $sam :$sam.variants.filter.vcf"
		done
	fi
	if [ $multi == "NO" ]
	then
		echo -e "\nDirectory with per sample Gene Summary reports (<sample>.Gene.Summary.txt): Reports_per_Sample/ANNOT/"
	else
		echo -e "\nDirectory with per sample Gene Summary reports (<group>.<TUMOR sample>.Gene.Summary.txt): Reports_per_Sample/ANNOT/"	
	fi
	
	if [ $tool == "whole_genome" ]
	then
		if [ $multi == "NO" ]
		then
			echo -e "\nDirectory with per sample SV and CNV reports\n-CNVs (<sample>.CNV.annotated.txt)\n-SVs (<sample>.SV.annotated.txt): Reports_per_Sample/ANNOT/"	
		else
			echo -e "\nDirectory with per sample SV and CNV reports\n-CNVs (<group>.<TUMOR sample>.CNV.annotated.txt)\n-SVs (<group>.<TUMOR sample>.SV.annotated.txt): Reports_per_Sample/ANNOT/"	
		fi
		if [ $multi == "NO" ]
		then
			echo -e "\nDirectory with Structural variation(SV) and Copy Number Variations(CNV) files in VCF format : Reports_per_Sample/SV/"
			echo -e "-CNVs\nRaw CNV calls in vcf format : <sample>.cnv.vcf"
			echo -e "Filtered CNV calls in vcf format : <sample>.cnv.filter.vcf"  
			echo -e "-SVs\nFiltered SV calls from BREAKDANCER in vcf format : <sample>.break.vcf"
			echo -e "Filtered SV calls from CREST in vcf format : <sample>.filter.crest.vcf"
		else
			echo -e "\nDirectory with SOMATIC Structural variation(SV) and SOMATIC Copy Number Variations(CNV) files in VCF format : Reports_per_Sample/SV/"
			echo -e "-CNVs\nRaw SOMATIC CNV calls in vcf format : <group>.<TUMOR sample>.cnv.vcf"
			echo -e "Filtered SOMATIC CNV calls in vcf format : <group>.<TUMOR sample>.cnv.filter.vcf"  
			echo -e "-SVs\nFiltered SOMATIC SV calls from BREAKDANCER in vcf format : <group>.<TUMOR sample>.somatic.break.vcf"
			echo -e "Filtered SOMATIC SV calls from CREST in vcf format : <group>.<TUMOR sample>.somatic.filter.crest.vcf"
		fi		
	fi		
	echo `date`
fi	

	
