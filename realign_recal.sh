#!/bin/sh
#	INFO
#	wrapper for 2A and 2B

if [ $# -le 5 ]
then
    echo -e "Usage: wrapper script to do realignment and variant calling \nMulti-Samples\n<input ':' sep> <bam ':' sep[normal:tumor1:tumor2:tumorN]> <samples ':' sep[normal:tumor1:tumor2:tumorN]> <outputdir bams> <run_info><1 for realign-recal or 0 for recal-realign>\nelse\n<input> <bam > <samples> <outputdir bams><run_info><1 for realign-recal or 0 for recal-realign>\n"
else
    set -x
    echo `date`
    input=$1    
    bam=$2
    samples=$3
    output_bam=$4
    run_info=$5
    flag=$6
	if [ $7 ]
	then
		SGE_TASK_ID=$7
	fi	
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	
	export JAVA_HOME=$javahome
	export PATH=$javahome/bin:$PATH
	
    ### update dash board    
    if [ $SGE_TASK_ID == 1 ]
    then
        $script_path/dashboard.sh $samples $run_info Realignment started
    fi    
    
    if [ $flag == 1 ]
    then
        $script_path/realign_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        if [ ! -s $output_bam/chr${chr}.realigned.bam ]
        then
            echo "ERROR: [`date`] realign_recal.sh File $output_bam/chr${chr}.realigned.bam not created" 
            exit 1
        fi

        $script_path/recal_per_chr.sh $output_bam chr${chr}.realigned.bam $output_bam $run_info 1 0 multi $chr
        if [ ! -s $output_bam/chr${chr}.cleaned.bam ]
        then
            echo "ERROR: [`date`] realign_recal.sh File $output_bam/chr${chr}.recalibrated.bam not created"
            exit 1
        fi
    else
        $script_path/recal_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        if [ ! -s $output_bam/chr${chr}.recalibrated.bam ]
        then
            echo "ERROR: [`date`] realign_recal.sh File $output_bam/chr${chr}.recalibrated.bam not created"
            exit 1
        fi

        $script_path/realign_per_chr.sh  $output_bam chr${chr}.recalibrated.bam $output_bam $run_info 1 0 multi $chr
        if [ ! -s $output_bam/chr${chr}.cleaned.bam ]
        then
            echo "ERROR: [`date`] realign_recal.sh File $output_bam/chr${chr}.realigned.bam not created"
            exit 1
        fi
    fi
	
    ## update the dash board
    if [ $SGE_TASK_ID == 1 ]
    then
		$script_path/dashboard.sh $samples $run_info Realignment complete
    fi    
    echo `date`
fi	