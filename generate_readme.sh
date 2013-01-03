#!/bin/bash

if [ $# != 2 ]
then
	echo "Usage: </path/to/output folder><path/to/run info file>"
else
	echo `date`
	output=$1
	run_info=$2
	file=$output/README
	##Reports Reports_per_Sample 
	multi=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
	somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
	tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[A-Z]" "[a-z]")
	analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" " ")
	groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2 | tr ":" " ")
	
	##
	echo -e "*********************************************\nREADME starts\n*********************************************" > $file
	echo -e "`date`" >> $file
	echo -e "\n###File Descriptions\n" >> $file
	echo -e "* Session file to visualize alignment BAMs, coverage and UCSC tracks using IGV: igv_session.xml" >> $file
	echo -e "* Help guide to set up and visualize in IGV: IGV_Setup.doc" >> $file
	echo -e "* Summary report on analysis performed and statistics obtained: Main_Document.html" >> $file
	echo -e "* Description of statistics table in Main_Document.html: StatisticsDescription.html" >> $file
	if [[ $multi == "YES"  && $somatic_calling == "YES" ]]
	then
		echo -e "* Plain text version of statistics table in Main_Document.html for single samples: SampleStatistics.pair.tsv" >> $file
	fi
	echo -e "* Plain text version of statistics table in Main_Document.html for single samples: SampleStatistics.tsv" >> $file
	echo -e "* Description of columns in SNV and INDEL reports: ColumnDescription_Reports.xls" >> $file
	if [  $tool == "whole_genome" ]
	then
		echo -e "* Whole genome workflow diagram: whole_genome_workflow.png" >> $file
	else
		echo -e "* Exome workflow diagram: exome_workflow.png"	 >> $file
	fi	
	if [[ $analysis != "annotation" && $analysis != "ontarget" && $analysis != "alignment" ]]
	then
		echo -e "* Per sample coverage plot: Coverage.JPG" >> $file
	fi	
	echo -e "\n* Compressed file with the statistical numbers in GENOME_GPS for all the samples : numbers.tar.gz" >> $file
  	
	echo -e "\n* Configuration files for GENOME_GPS (under config folder):\n" >> $file
  	echo -e "-run_info.txt\n-sample_info.txt\n-tool_info.txt\n-memory_info.txt" >> $file
 	echo -e "\n#### Annotation Reports\n" >> $file
	echo -e "Directory with merged SNV and INDEL reports: : Reports/" >> $file
	echo -e "\n**ANNOTATION files (filtered - single annotation line per variant call(most impacting transcript))" >> $file
	echo -e "\n-SNVs\nSNV.xls\nSNV.filtered.xls\n\n-INDELs\nINDEL.xls\nINDEL.filtered.xls" >> $file 
	if [[ $multi == "YES" && $somatic_calling == "YES" ]]
	then
		echo -e "\n-SNVs (somatic merged report)\nTUMOR.SNV.xls\nTUMOR.SNV.filtered.xls\n\n-INDELs (somatic merged report)\nTUMOR.INDEL.xls\nTUMOR.INDEL.filtered.xls" >> $file
	fi	
	
	echo -e "\nDirectory with per sample or per group SNV, INDEL reports: Reports_per_Sample/" >> $file
	echo -e "\n-SNVs" >> $file
	if [ $multi == "NO" ]
	then
		for sam in $samples
		do
			echo -e "\nsample : $sam" >> $file
			echo -e "$sam.SNV.xls - There could be multiple annotation lines for each variant for $sam\n$sam.SNV.filtered.xls - just one annotation for each variant for $sam (most impacting one)" >> $file
		done
	else
		if [ $soamtic_calling == "NO" ]
		then
			for pair in $groups
			do
				echo -e "\ngroup : $pair" >> $file
				echo -e "$pair.SNV.xls - There could be multiple annotation lines for each variant for $pair\n$pair.SNV.filtered.xls - just one annotation for each variant for $pair (most impacting one)" >> $file
			done
		else
			for pair in $groups
			do
				echo -e "\ngroup : $pair" >> $file
				echo -e "$pair.SNV.xls - There could be multiple annotation lines for each variant for $pair\n$pair.SNV.filtered.xls - just one annotation for each variant for $pair (most impacting one)" >> $file
				echo -e "TUMOR.$pair.SNV.xls - There could be multiple annotation lines for each somatic variant for $pair\nTUMOR.$pair.SNV.filtered.xls - just one annotation for each somatic variant for $pair (most impacting one)" >> $file
			done
		fi	
	fi	
	echo -e "\n-INDELs" >> $file
	if [ $multi == "NO" ]
	then
		for sam in $samples
		do
			echo -e "\nsample : $sam" >> $file
			echo -e "$sam.INDEL.xls - There could be multiple annotation lines for each variant for $sam\n$sam.INDEL.filtered.xls - just one annotation for each variant for $sam (most impacting one)" >> $file
		done
	else
		for pair in $groups
		do
			echo -e "\ngroup : $pair" >> $file
			echo -e "$pair.INDEL.xls - There could be multiple annotation lines for each variant for $pair\n$pair.INDEL.filtered.xls - just one annotation for each variant for $pair (most impacting one)" >> $file 
			if [ $somatic_calling == "YES" ]
			then
				echo -e " TUMOR.$pair.INDEL.xls - There could be multiple annotation lines for each somatic variant for $pair\nTUMOR.$pair.INDEL.filtered.xls - just one annotation for each somatic variant for $pair (most impacting one)" >> $file
			fi	
		done
	fi	
	
	echo -e "\n**VCF files (Variant Calling format files. Could be opened in IGV)" >> $file
	if [ $multi == "YES" ]
	then
		for pair in $groups
		do
			echo -e "\ngroup : $pair" >> $file
			if [ $somatic_calling == "YES" ]
			then
				echo -e "filtered somatic variant file for $pair : $pair.somatic.variants.filter.vcf" >> $file
			fi
			echo -e "filtered multi sample variant file for $pair : $pair.variants.filter.vcf" >> $file
		done
	else
		for sam in $samples
		do			
			echo -e "\nsample : $sam" >> $file
			echo -e "filtered variant file for $sam :$sam.variants.filter.vcf" >> $file
		done
	fi
	if [ $multi == "NO" ]
	then
		echo -e "\nDirectory with per sample Gene Summary reports (<sample>.Gene.Summary.txt): Reports_per_Sample/ANNOT/" >> $file
	else
		if [ $somatic_calling == "YES" ]
		then
			echo -e "\nDirectory with per sample Gene Summary reports (<group>.<TUMOR sample>.Gene.Summary.txt): Reports_per_Sample/ANNOT/"	 >> $file
		else
			echo -e "\nDirectory with per sample Gene Summary reports (<sample>.Gene.Summary.txt): Reports_per_Sample/ANNOT/" >> $file
		fi	
	fi
	
	if [ $tool == "whole_genome" ]
	then
		if [ $multi == "NO" ]
		then
			echo -e "\nDirectory with per sample SV and CNV reports\n-CNVs (<sample>.CNV.annotated.txt)\n-SVs (<sample>.SV.annotated.txt): Reports_per_Sample/ANNOT/"	 >> $file
		else
			if [ $somatic_calling == "YES" ]
			then
				echo -e "\nDirectory with per sample SV and CNV reports\n-CNVs (<group>.<TUMOR sample>.CNV.annotated.txt)\n-SVs (<group>.<TUMOR sample>.SV.annotated.txt): Reports_per_Sample/ANNOT/"	 >> $file
			else
				echo -e "\nDirectory with per sample SV and CNV reports\n-CNVs (<sample>.CNV.annotated.txt)\n-SVs (<sample>.SV.annotated.txt): Reports_per_Sample/ANNOT/"	 >> $file
			fi
		fi
		if [ $multi == "NO" ]
		then
			echo -e "\nDirectory with Structural variation(SV) and Copy Number Variations(CNV) files in VCF format : Reports_per_Sample/SV/" >> $file
			echo -e "-CNVs\nRaw CNV calls in vcf format : <sample>.cnv.vcf" >> $file
			echo -e "Filtered CNV calls in vcf format : <sample>.cnv.filter.vcf"   >> $file
			echo -e "-SVs\nFiltered SV calls from BREAKDANCER in vcf format : <sample>.break.vcf" >> $file
			echo -e "Filtered SV calls from CREST in vcf format : <sample>.filter.crest.vcf" >> $file
		else
			if [ $somatic_calling == "YES" ]
			then
				echo -e "\nDirectory with SOMATIC Structural variation(SV) and SOMATIC Copy Number Variations(CNV) files in VCF format : Reports_per_Sample/SV/" >> $file
				echo -e "-CNVs\nRaw SOMATIC CNV calls in vcf format : <group>.<TUMOR sample>.cnv.vcf" >> $file
				echo -e "Filtered SOMATIC CNV calls in vcf format : <group>.<TUMOR sample>.cnv.filter.vcf"   >> $file
				echo -e "-SVs\nFiltered SOMATIC SV calls from BREAKDANCER in vcf format : <group>.<TUMOR sample>.somatic.break.vcf" >> $file
				echo -e "Filtered SOMATIC SV calls from CREST in vcf format : <group>.<TUMOR sample>.somatic.filter.crest.vcf">> $file
			else
				echo -e "\nDirectory with Structural variation(SV) and Copy Number Variations(CNV) files in VCF format : Reports_per_Sample/SV/" >> $file
				echo -e "-CNVs\nRaw CNV calls in vcf format : <sample>.cnv.vcf" >> $file
				echo -e "Filtered CNV calls in vcf format : <sample>.cnv.filter.vcf"   >> $file
				echo -e "-SVs\nFiltered SV calls from BREAKDANCER in vcf format : <sample>.break.vcf" >> $file
				echo -e "Filtered SV calls from CREST in vcf format : <sample>.filter.crest.vcf" >> $file
			fi
		fi		
	fi		
	TO=`id |awk -F '(' '{print $2}' | cut -f1 -d ')'`
	email=`finger $TO | head -1 | cut -d ';' -f2`
	name=`finger $TO | head -1 | cut -d ';' -f1 | awk -F':' '{print $NF}'|sed 's/^[ \t]*//;s/[ \t]*$//'`
	echo -e "\nFurther Questions:\nContact:\n$name\n$email\n" >> $file
	echo -e "`date`" >> $file
	echo -e "*********************************************\nREADME Ends\n*********************************************" >> $file
	echo `date`
fi	

	
