#!/bin/sh
##	INFO
##	to merge the scripts to call in one script to control sapnning jobs
if [ $# -le 7 ]
then
    echo "Usage: <run info> <sample> <Temp reports> <Output OnTarget> <sift> <snpeff> <polyphen> <output dir><group optional>";
else
    set -x
    echo `date`			
    run_info=$1 
    sample=$2
    TempReports=$3 
    output_OnTarget=$4 
    sift=$5 
    snpeff=$6
    poly=$7
    output_dir=$8
	if [ $9 ]
	then
		group=$9
    fi
	#SGE_TASK_ID=1
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    which_chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    
	
    ## prepocessing the input file from variant module or user added 
    $script_path/preprocess.persample.sh $sample $TempReports $run_info $output_OnTarget $which_chr $group
    
    if [ $variant_type == "BOTH" -o $variant_type == "SNV" ]
    then
        if [ $9 ]
		then
			snv_var=$group.${sample}.chr${which_chr}.snv
        else
			snv_var=${sample}.chr${which_chr}.snv
		fi
		## add rsids
        $script_path/add.rsids_snvs.sh $TempReports $snv_var $which_chr $run_info
        ## add allele frequency
        $script_path/add.frequencies.sh $TempReports $snv_var $which_chr $run_info
        ## merge sift sseq codon UCSC tracks
        if [ $9 ]
		then
			$script_path/merge.snv.sh $TempReports $group.$sample $which_chr $sift $snpeff $poly $snv_var $run_info
		else
			$script_path/merge.snv.sh $TempReports $sample $which_chr $sift $snpeff $poly $snv_var $run_info
		fi
	fi
	
    if [ $variant_type == "BOTH" -o $variant_type = "INDEL" ]
    then	
		if [ $9 ]
		then
			indel_var=$group.${sample}.chr${which_chr}.indel
		else
			indel_var=${sample}.chr${which_chr}.indel
        fi
		## merge sseq to indel report
		$script_path/add.rsids_indels.sh $TempReports $indel_var $which_chr $run_info
		if [ $9 ]
		then
			$script_path/merge.indel.sh $TempReports $group.$sample $which_chr $snpeff $indel_var $run_info
		else
			$script_path/merge.indel.sh $TempReports $sample $which_chr $snpeff $indel_var $run_info	
		fi
	fi
    echo `date`
fi	
    



            
