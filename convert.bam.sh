#!/bin/sh

########################################################
###### 	SORT BAM, adds RG & REMOVE DUPLICATES SCRIPT FOR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			convert.bam.sh
######		Date:				07/27/2011
######		Summary:			Using PICARD to sort and mark duplicates in bam 
######		Input files:		$1	=	/path/to/input directory
######							$2	=	name of BAM to sort
######							$3	=	sample name
######							$4	=	/path/to/run_info.txt
######		Output files:		Sorted and clean BAM
######		TWIKI:				http://bioinformatics.mayo.edu/BMI/bin/view/Main/BioinformaticsCore/Analytics/WholeGenomeWo
########################################################

if [ $# != 4 ];
then
    echo -e "Usage: wrapper to add read group and sort the bam </path/to/input directory> <name of BAM to sort> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    input_bam=$2
    sample=$3
    run_info=$4
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2| tr "[A-Z]" "[a-z]")
	samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    dup=$( cat $run_info | grep -w '^MARKDUP' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    reorder=$( cat $run_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    version=$( cat $run_info | grep -w '^VERSION' | cut -d '=' -f2)
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
    output=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    flowcell=`echo $run_num | awk -F'_' '{print $NF}' | sed 's/.\(.*\)/\1/'`
    
########################################################	
######		PICARD to sort raw BAM file

    if [ ! -s $input/$input_bam ]
    then
        echo "ERROR: $0 File $input/$input_bam does not exist"
        exit 1
    fi

    ## check if BAM is sorted
    SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$input_bam -s $samtools`
    if [ $SORT_FLAG == 1 ]
    then
        ln -s $input/$input_bam $input/$sample.sorted.bam
    else
        $java/java -Xmx6g -Xms512m \
        -jar $picard/SortSam.jar \
        INPUT=$input/$input_bam \
        OUTPUT=$input/$sample.sorted.bam \
        MAX_RECORDS_IN_RAM=2000000 \
        SO=coordinate \
        TMP_DIR=$input/ \
        VALIDATION_STRINGENCY=SILENT
        
        if [ -s $input/$sample.sorted.bam ]
        then
            rm $input/$sample.bam
        else
            echo "ERROR: convert.bam.sh File $input/$sample.sorted.bam not generated"
            exit 1
        fi
    fi

#############################################################	
######		PICARD to check availability of ReadGroup and platform info


    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`
    sam=`echo $sample | awk -F'.' '{print $1}'`
    if [ "$RG_ID" != "$sam" ]
    then
        $java/java -Xmx6g -Xms512m \
        -jar $picard/AddOrReplaceReadGroups.jar \
        INPUT=$input/$sample.sorted.bam OUTPUT=$input/$sample.sorted.rg.bam \
        PL=$platform SM=$sample CN=$center ID=$sample PU=$sample LB=$GenomeBuild \
        TMP_DIR=$input \
        VALIDATION_STRINGENCY=SILENT

        if [ -s $input/$sample.sorted.rg.bam ]
        then
            mv $input/$sample.sorted.rg.bam $input/$sample.sorted.bam
        else
            echo "ERROR: convert.bam.sh File $input/$sample.sorted.rg.bam not generated"
            exit 1
        fi
    fi


########################################################	
######		PICARD to mark duplicates
    if [ $dup == "YES" ]
    then
        $java/java -Xmx6g -Xms512m \
        -jar $picard/MarkDuplicates.jar \
        INPUT=$input/$sample.sorted.bam \
        OUTPUT=$input/$sample.sorted.rmdup.bam \
        METRICS_FILE=$input/$sample.dup.metrics \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$input/

        if [ -s $input/$sample.sorted.rmdup.bam ]
        then
            mv $input/$sample.sorted.rmdup.bam $input/$sample.sorted.bam
        else
            echo "ERROR: convert.bam.sh $input/$sample.sorted.rmdup.bam not generated"
            exit 1
        fi
    fi	

    $samtools/samtools index $input/$sample.sorted.bam $input/$sample.sorted.bam.bai

    if [ ! -s $input/$sample.sorted.bam.bai ]
    then
        echo "ERROR: convert.bam.sh File $input/$sample.sorted.bam.bai not generated"
        exit 1
    fi

########################################################	
######		PICARD to reorder BAM/SAM
    if [ $reorder == "YES" ]
    then
        $java/java -Xmx6g -Xms512m \
        -jar $picard/ReorderSam.jar \
        I=$input/$sample.sorted.bam \
        O=$input/$sample.sorted.tmp.bam \
        R=$ref_path \
        VALIDATION_STRINGENCY=SILENT \
        CREATE_INDEX=true

        if [ -s $input/$sample.sorted.tmp.bam ]
        then
            echo "ERROR: convert.bam.sh File $input/$sample.sorted.tmp.bam not generated"
            exit 1
        else
            mv $input/$sample.sorted.tmp.bam $input/$sample.sorted.bam
        fi

        if [ -s $input/$sample.sorted.tmp.bam.bai ]
        then
            echo "ERROR: convert.bam.sh File $input/$sample.sorted.tmp.bam.bai not generated"
            exit 1
        else
            mv $input/$sample.sorted.tmp.bai $input/$sample.sorted.bam.bai
        fi
    fi
    
    ### flagstat on each bam file
    $samtools/samtools flagstat $input/$sample.sorted.bam > $input/$sample.flagstat
    if [ ! -s $input/$sample.flagstat ]
    then
        echo "ERROR: convert.bam.sh flagstat failed for $input/$sample.flagstat"
    fi
    ## update secondary dahboard
    if [ $analysis == "mayo" -o $analysis == "realign-mayo"  ]
	then
		id=`echo $sample | awk -F'.' '{print $NF}'`
		sam=`echo $sample | awk -F'.' '{print $1}'`
		pos=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | grep -n $sam | cut -d ":" -f1)
		lane=$( cat $run_info | grep -w '^LANEINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $SGE_TASK_ID | tail -n 1)
		index=$( cat $run_info | grep -w '^LABINDEXES' | cut -d '=' -f2 | tr ":" "\n" | head -n $pos | tr "," "\n" | head -n $SGE_TASK_ID | tail -n 1)
        if [ $index == "-" ]
        then
            $java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -r $run_num -s Alignment -a WholeGenome -v $version
        else
            $java/java -jar $script_path/AddSecondaryAnalysis.jar -p $script_path/AddSecondaryAnalysis.properties -l $lane -c -f $flowcell -i $index -r $run_num -s Alignment -a WholeGenome -v $version
        fi    
    fi
    echo `date`
fi
