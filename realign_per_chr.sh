#!/bin/bash
## this scripts work per chr and accepts an array job paramter to extract the chr information
## this scripts checks for sorted and rad group information for a abam and do as per found
## GATK version using GenomeAnalysisTK-1.2-4-gd9ea764
## here we consider if chopped is 1 means all the sample BAM are chopped and same with 0 
 
if [ $# != 8 ]
then
    echo -e "Usage:\nIf user wants to do realignment fist \n<input dir ':' sep><input bam ':' sep><outputdir><run_info><1 or 0 if bam is per chr><1 for realign first><sample ':' sep>\nelse\n<input dir><input bam><output dir><run_info> <1 or 0 if bam is per chr> < 0 for realign second><sample(add a dummy sample name as we dont care about the sample name (example:multi))>  ";
else	
    set -x
    echo `date`
    input=$1    
    bam=$2
    output=$3
    run_info=$4
    chopped=$5
    realign=$6
    samples=$7
    chr=$8
    # get job array ID
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    maxreads=$( cat $tool_info | grep -w '^MAX_READS_REALIGN' | cut -d '=' -f2 ) 
    maxreadsmem=$( cat $tool_info | grep -w '^MAX_READS_MEM_REALIGN' | cut -d '=' -f2 ) 
    
    if [ ${#dbSNP} -ne 0 ]
    then
        param="-known $dbSNP" 
    fi
    
    if [ ${#Kgenome} -ne 0 ]
    then
        param=$param" -known $Kgenome" 
    fi
    
    

    if [ $realign == 1 ]
    then
        inputDirs=$( echo $input | tr ":" "\n" )
        bamNames=$( echo $bam | tr ":" "\n" )
        sampleNames=$( echo $samples | tr ":" "\n" )
        i=1
        for inp in $inputDirs
        do
            inputArray[$i]=$inp
            let i=i+1
        done
        i=1
        for ba in $bamNames
        do
            bamArray[$i]=$ba
            let i=i+1
        done
        i=1
        for sa in $sampleNames
        do
            sampleArray[$i]=$sa
            let i=i+1
        done
        if [ ${#inputArray[@]} != ${#bamArray[@]} -o ${#inputArray[@]} != ${#sampleArray[@]} ]
        then
            echo "ERROR : realign_per_chr ':' sep parameters are not matching" 
            exit 1;
        else    
            for i in $(seq 1 ${#sampleArray[@]})
            do
                sample=${sampleArray[$i]}
                input=${inputArray[$i]}
                bam=${bamArray[$i]}
                ##extracting and checking the BAM for specific chromosome
    
                if [ ! -s $input/$bam ]
                then
					$script_path/errorlog.sh $input/$bam realign_per_chr.sh ERROR "does not exist"
                    exit 1;
                fi
                $script_path/samplecheckBAM.sh $input $bam $output $run_info $sample $chopped $chr
            done
        fi
        input_bam=""
        for i in $(seq 1 ${#sampleArray[@]})
        do
            input_bam="${input_bam} -I $output/${sampleArray[$i]}.chr${chr}-sorted.bam"
        done
    else
        if [ ! -s $input/$bam ]
        then
            $script_path/errorlog.sh $input/$bam realign_per_chr.sh ERROR "does not exist"
            exit 1;
        fi
        $script_path/samplecheckBAM.sh $input $bam $output $run_info $samples $chopped $chr
        input_bam="-I $output/$samples.chr${chr}-sorted.bam"
    fi	
    
	if [ ! -d $output/temp/ ]
	then
		mkdir -p $output/temp/
	fi
	
    ## GATK Target Creator
    $java/java -Xmx5g -Xms512m -Djava.io.tmpdir=$output/temp/ -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -K $gatk/Hossain.Asif_mayo.edu.key \
    $param -L chr${chr} \
    -T RealignerTargetCreator \
    $input_bam \
    -o $output/chr${chr}.forRealigner.intervals

    if [ ! -s $output/chr${chr}.forRealigner.intervals ]
    then
        echo "WARNING : realign_per_chr. File $output/chr${chr}.forRealigner.intervals not created"
        bams=`echo $input_bam | sed -e '/-I/s///g'`
        num_bams=`echo $bams | tr " " "\n" | wc -l`
        if [ $num_bams -eq 1 ]
        then
            cp $bams $output/chr${chr}.realigned.bam
            cp $bams.bai $output/chr${chr}.realigned.bam.bai
        else
            INPUTARGS=`echo $bams | tr " " "\n" | awk '{print "I="$1}'` 
            $script_path/MergeBam.sh $INPUTARGS $output/chr${chr}.realigned.bam $output true $run_info 
		fi
    else
        ## Realignment
        $java/java -XX:+UseConcMarkSweepGC -Xmx5g -Xms512m -Djava.io.tmpdir=$output/temp/ \
        -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -T IndelRealigner \
        -K $gatk/Hossain.Asif_mayo.edu.key \
        $param -L chr${chr} \
        $input_bam \
        --maxReadsForRealignment $maxreads \
        --maxReadsInMemory $maxreadsmem \
        --out $output/chr${chr}.realigned.bam  \
        -targetIntervals $output/chr${chr}.forRealigner.intervals
        mv $output/chr${chr}.realigned.bai $output/chr${chr}.realigned.bam.bai
    fi

    if [ -s $output/chr${chr}.realigned.bam ]
    then
        if [ $realign == 0 ]
		then
            cp $output/chr${chr}.realigned.bam	$output/chr${chr}.cleaned.bam
            cp $output/chr${chr}.realigned.bam.bai $output/chr${chr}.cleaned.bam.bai
            $samtools/samtools flagstat $output/chr${chr}.cleaned.bam > $output/chr${chr}.flagstat
		fi	
    else
        $script_path/errorlog.sh $output/chr${chr}.realigned.bam realign_per_chr.sh ERROR "does not exist"
        exit 1;
    fi

    ## deleting the internediate files
    if [ $realign == 1 ]
    then
        for i in $(seq 1 ${#sampleArray[@]})
        do
            rm $output/${bamArray[$i]}.$chr.bam
            rm $output/${bamArray[$i]}.$chr.bam.bai
            rm $output/${sampleArray[$i]}.chr${chr}.bam
            rm $output/${sampleArray[$i]}.chr${chr}.bam.bai
            rm $output/${sampleArray[$i]}.chr${chr}-sorted.bam
            rm $output/${sampleArray[$i]}.chr${chr}-sorted.bam.bai
        done
    else
        rm $output/$bam.$chr.bam
        rm $output/$bam.$chr.bam.bai
        rm $output/$samples.chr${chr}.bam
        rm $output/$samples.chr${chr}.bam.bai
        rm $output/$samples.chr${chr}-sorted.bam
        rm $output/$samples.chr${chr}-sorted.bam.bai
    fi		
    rm $output/chr${chr}.forRealigner.intervals
    echo  `date`	
fi
