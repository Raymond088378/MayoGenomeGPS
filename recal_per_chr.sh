#!/bin/bash
## this scripts work per chr and accepts an array job paramter to extract the chr information
## creat a folder name temp in the output folder befor euisng this script
## GATK version using GenomeAnalysisTK-1.2-4-gd9ea764
if [ $# != 8 ]
then
    echo -e "Usage:\nIf user wants to do recalibration fist \n<input dir ':' sep><input bam ':' sep><outputdir><run_info><1 or 0 if bam is per chr><1 for recalibrate first ><sample ':' sep>\nelse\n<input dir><input bam><output dir><run_info> <1 or 0 if bam is per chr> < 0 for recal second><sample (a dummy sampel name i would say just type multi as sample>  ";
else	
    set -x
    echo `date`
    input=$1    
    bam=$2
    output=$3
    run_info=$4
    chopped=$5
    recal=$6
    samples=$7
    chr=$8

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2|tr "[a-z]" "[A-Z]")
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    
    if [ ${#dbSNP} -ne 0 ]
    then
        param="--knownSites $dbSNP" 
    fi
    
    if [ ${#Kgenome} -ne 0 ]
    then
        param=$param" --knownSites $Kgenome" 
    fi
    
    
    if [ $recal == 1 ]
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
        
        i=1
        for sa in $sampleNames
        do
            sampleArray[$i]=$sa
            let i=i+1
        done
        if [ ${#inputArray[@]} != ${#bamArray[@]} -o ${#inputArray[@]} != ${#sampleArray[@]} ]
        then
            echo "ERROR : ':' sep parameters are not matching check the $run_info file";
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
                    $script_path/errorlog.sh $input/$bam recal_per_chr.sh ERROR "does not exist"
                    exit 1
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
	
	if [ $tool == "whole_genome" ]
    then
    	region="-L chr${chr}"
    else
    	cat $TargetKit | grep -w chr$chr > $output/chr$chr.bed
    	region="-L $output/chr$chr.bed"
    fi	
    
    ## Recal metrics file creation
    $java/java -Xmx5g -Xms512m -Djava.io.tmpdir=$output/temp/ -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -K $gatk/Hossain.Asif_mayo.edu.key \
    $param $input_bam \
    $region \
    -T CountCovariates \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate \
    -recalFile $output/chr${chr}.recal_data.csv 

    if [ ! -s $output/chr${chr}.recal_data.csv ]
    then
        echo "WARNING : recal_per_chr. File $output/chr${chr}.recal_data.csv not created"
        bams=`echo $input_bam | sed -e '/-I/s///g'`
        num_bams=`echo $bams | tr " " "\n" | wc -l`
        if [ $num_bams -eq 1 ]
        then
            cp $bams $output/chr${chr}.recalibrated.bam
            cp $bams.bai $output/chr${chr}.recalibrated.bam.bai
        else
            INPUTARGS=`echo $bams | tr " " "\n" | awk '{print "I="$1}'` 
            $script_path/MergeBam.sh $INPUTARGS $output/chr${chr}.recalibrated.bam $output true $run_info 
        fi
    else	
        rm $output/chr$chr.bed
    	## recailbartion
        $java/java -Xmx5g -Xms512m -Djava.io.tmpdir=$output/temp/ -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -K $gatk/Hossain.Asif_mayo.edu.key \
        -L chr${chr} \
        $input_bam \
        -T TableRecalibration \
        --out $output/chr${chr}.recalibrated.bam \
        -recalFile $output/chr${chr}.recal_data.csv 
        mv $output/chr${chr}.recalibrated.bai $output/chr${chr}.recalibrated.bam.bai
    fi
    
    if [ -s $output/chr${chr}.recalibrated.bam ]
    then
        rm $input/$bam $input/$bam.bai
        if [ $recal == 0 ]
        then
            mv $output/chr${chr}.recalibrated.bam $output/chr${chr}.cleaned.bam
            mv $output/chr${chr}.recalibrated.bam.bai $output/chr${chr}.cleaned.bam.bai
            $samtools/samtools flagstat $output/chr${chr}.cleaned.bam > $output/chr$chr.flagstat
        fi		
    else
        $script_path/errorlog.sh output/chr${chr}.recalibrated.bam recal_per_chr.sh ERROR "does not exist"
        exit 1;
    fi
    
    ## deleting internediate files
    if [ $recal == 1 ]
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
    rm $output/chr${chr}.recal_data.csv 		
    echo `date`	
fi
