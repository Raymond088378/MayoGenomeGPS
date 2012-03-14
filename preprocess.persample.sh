#!/bin/sh
##	INFO
##	script is a wrapper script to annotate variants per sample 

#########################################
#	$1		=		Sample Name
#	$2		=		SIFT directory	
#	$3		=		SSEQ directory	
#	$4		=		TempReports directory
#	$5		=		run information file to get filenames
#	$6		=		input directory run dir
#	$7		=		chromosome
###############################################	

if [ $# != 7 ];
then
    echo "Usage:<sample> <sift dir> <sseq dir> <configuration file> <tempReport dir><email><analysis type><run info><reports per sample> <input variant folder><chromosome>";
else			
    set -x
    echo `date`
    sample=$1
    sift=$2 
    sseq=$3 		
    TempReports=$4 	#Tempreport output folder
    run_info=$5
    input_dir=$6
    which_chr=$7
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    SNV_caller=$( cat $run_info | grep -w '^SNV_CALLER' | cut -d '=' -f2 )
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    GeneIdMap=$( cat $tool_info | grep -w '^GeneIdMap' | cut -d '=' -f2)
    analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    
   
	if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
	then
		## SNVs
		echo -e "\t\t\t\t\t${sample}" >> $TempReports/${sample}.chr${which_chr}.temp.file
		var=${sample}.chr${which_chr}.raw.snvs.bed.i.ToMerge
		if [ $SNV_caller == "SNVmix" ]
		then
			echo -e "Chr\tPosition\tInCaptureKit\tRef\tAlt\tGenotypeClass\tAlt-SupportedReads\tRef-SupportedReads\tReadDepth\tProbability\tCloseToIndel" >> $TempReports/${sample}.chr${which_chr}.temp.file
		elif [ $SNV_caller == "GATK" ]
		then
			echo -e "Chr\tPosition\tInCaptureKit\tRef\tAlt\tGenotypeClass\tAlt-SupportedReads\tRef-SupportedReads\tReadDepth\tQuality\tCloseToIndel" >> $TempReports/${sample}.chr${which_chr}.temp.file
		else
			echo -e "Chr\tPosition\tInCaptureKit\tRef\tAlt" >> $TempReports/${sample}.chr${which_chr}.temp.file
		fi
		
		cat $input_dir/$var | awk 'BEGIN { OFS="\t"} {print $1,$2,$10,$3,$4,$5,$6,$7,$8,$9,$11}' > $input_dir/$var.tmp
		cat $input_dir/$var.tmp >> $TempReports/${sample}.chr${which_chr}.temp.file
		rm $input_dir/$var.tmp
		mv $TempReports/${sample}.chr${which_chr}.temp.file $TempReports/${sample}.chr${which_chr}.snv
	fi
	if [ $variant_type == "BOTH" -o $variant_type == "INDEL" ]
	then
		## INDELs
		indel=${sample}.chr${which_chr}.raw.indels.bed.i.ToMerge
		echo -e "\t\t\t\t\t\t\t${sample}" >> $TempReports/${sample}.chr${which_chr}.temp.file
		echo -e "Chr\tStart\tStop\tInCaptureKit\tRef\tAlt\tBase-Length\tIndel-supportedRead\tReadDepth" >> $TempReports/${sample}.chr${which_chr}.temp.file
		cat $input_dir/$indel | awk 'BEGIN { OFS="\t"} {print $1,$2,$3,$9,$4,$5,$6,$7,$8}' > $input_dir/$indel.tmp
		cat $input_dir/$indel.tmp >> $TempReports/${sample}.chr${which_chr}.temp.file
		mv $TempReports/${sample}.chr${which_chr}.temp.file $TempReports/$sample.chr${which_chr}.indel
		rm $input_dir/$indel.tmp
	fi
    echo `date`
fi	