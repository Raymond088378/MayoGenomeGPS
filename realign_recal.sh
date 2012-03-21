#!/bin/sh
#	INFO
#	wrapper for 2A and 2B

if [ $# != 7 ]
then
    echo -e "Usage: wrapper script to do realiagnemnt and varaint calling \nMulti-Samples\n<input ':' sep> <bam ':' sep[normal:tumor1:tumor2:tumorN]> <samples ':' sep[normal:tumor1:tumor2:tumorN]> <outputdir bams> <outputdir variants> <run_info><1 for realign-recal or 0 for recal-realign>\nelse\n<input> <bam > <samples> <outputdir bams> <outputdir variants> <run_info><1 for realign-recal or 0 for recal-realign>\n"
else
    set -x
    echo `date`
    input=$1    
    bam=$2
    samples=$3
    output_bam=$4
    output_var=$5
    run_info=$6
    flag=$7
   
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
	java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    bed=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    out_dir=$out/$PI/$tool/$run_num
    PATH=$bed/:$PATH
	
        
    if [ $SGE_TASK_ID == 1 ]
    then
        if [ $analysis == "mayo" ]
        then
            s=`echo $samples | tr ":" " "`
            for sam in $s
            do
                pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sam | cut -d ":" -f1)
                lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
				i=1
				for lane in $lanes
				do
					index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
					if [ $index == "-" ]
					then
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -r $run_num -s Realignment -a WholeGenome -v $version
					else
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -f $flowcell -i $index -r $run_num -s Realignment -a WholeGenome -v $version
					fi
					let i=i+1	
				done		
            done
        fi
    fi    
    

    if [ $flag == 1 ]
    then
        $script_path/realign_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        if [ ! -s $output_bam/chr${chr}.realigned.bam ]
        then
            echo "ERROR: realign_wrapper.sh File $output_bam/chr${chr}.realigned.bam not created" 
            exit 1
        fi

        $script_path/recal_per_chr.sh $output_bam chr${chr}.realigned.bam $output_bam $run_info 1 0 multi $chr
        if [ ! -s $output_bam/chr${chr}.cleaned.bam ]
        then
            echo "ERROR: realign_wrapper File $output_bam/chr${chr}.recalibrated.bam not created"
            exit 1
        fi
    else
        $script_path/recal_per_chr.sh $input $bam $output_bam $run_info 0 1 $samples $chr
        if [ ! -s $output_bam/chr${chr}.recalibrated.bam ]
        then
            echo "ERROR: realign_wrapper File $output_bam/chr${chr}.recalibrated.bam not created"
            exit 1
        fi

        $script_path/realign_per_chr.sh  $output_bam chr${chr}.recalibrated.bam $output_bam $run_info 1 0 multi $chr
        if [ ! -s $output_bam/chr${chr}.cleaned.bam ]
        then
            echo "ERROR: realign_wrapper File $output_bam/chr${chr}.realigned.bam not created"
            exit 1
        fi
    fi
	
	## update the dash board
	if [ $SGE_TASK_ID == 1 ]
    then
        if [ $analysis == "mayo" ]
        then
            s=`echo $samples | tr ":" " "`
            for sam in $s
            do
                pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sam | cut -d ":" -f1)
                lanes=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," " ")
                i=1
				for lane in $lanes
				do
					index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tail -n 1 | tr "," "\n" | head -n $i | tail -n 1)
					if [ $index == "-" ]
					then
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -r $run_num -s Realignment -a WholeGenome -v $version
					else
						$java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -i $index -r $run_num -s Realignment -a WholeGenome -v $version
					fi
					let i=i+1
				done		
            done
        fi
    fi    
	echo `date`
fi	