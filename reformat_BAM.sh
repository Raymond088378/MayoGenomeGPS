#!/bin/sh
if [ $# != 3 ];
then
    echo -e "Usage: wrapper to merge bam files and validate the bam for downstream analysis \n merge_align.bam.sh </path/to/input directory> <name of BAM to sort> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    input=$1
    sample=$2
    run_info=$3
	
########################################################	
######		Reading run_info.txt and assigning to variables
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 )
    reorder=$( cat $run_info | grep -w '^REORDERSAM' | cut -d '=' -f2| tr "[a-z]" "[A-Z]")
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2 )
    ref_path=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
    platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
    GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
	
	########################################################	
######		PICARD to merge BAM file
    INPUTARGS="";
    files=""
    cd $input
    for file in $input/*sorted.bam
    do
      INPUTARGS="INPUT="$file" "$INPUTARGS;
      files=$file" $files";
    done
    
	num_times=`echo $INPUTARGS | tr " " "\n" | grep -c -w 'INPUT'`
	if [ $num_times == 1 ]
	then
		bam=`echo $INPUTARGS | cut -d '=' -f2`
		mv $bam $input/$sample.bam
        SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $input/$sample.bam -s $samtools`
        if [ $SORT_FLAG == 1 ]
        then
            mv $input/$sample.bam $input/$sample.sorted.bam
        else
            $java/java -Xmx6g -Xms512m \
            -jar $picard/SortSam.jar \
            INPUT=$input/$sample.bam \
            OUTPUT=$input/$sample.sorted.bam \
            MAX_RECORDS_IN_RAM=2000000 \
            SO=coordinate \
            TMP_DIR=$input/ \
            VALIDATION_STRINGENCY=SILENT
            if [ -s $input/$sample.sorted.bam ]
            then
                rm $input/$sample.bam
            else
                echo "sorting failed for $sample"
            fi
        fi
		$samtools/samtools index $input/$sample.sorted.bam
	else	
		$java/java -Xmx6g -Xms512m \
		-jar $picard/MergeSamFiles.jar \
		$INPUTARGS\
		OUTPUT=$input/$sample.sorted.bam \
		TMP_DIR=$input \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY="SILENT"
	fi

    if [ ! -s $input/$sample.sorted.bam ]
    then
        echo "ERROR: merge_align.bam.sh File $input/$sample.sorted.bam was not created"
        exit 1
    else
        if [ $num_times -gt 1 ]
		then
			rm $files
        fi
        if [ $num_times -lt 1 ]
		then
			mv $input/$sample.sorted.bai $input/$sample.sorted.bam.bai
		fi
	fi
    
    ### add read grouup information
    RG_ID=`$samtools/samtools view -H $input/$sample.sorted.bam | grep "^@RG" | tr '\t' '\n' | grep "^ID"| cut -f 2 -d ":"`

    if [ "$RG_ID" == "$sample" ]
    then
		echo "no need to convert same read group"
    else	
        $java/java -Xmx6g -Xms512m \
        -jar $picard/AddOrReplaceReadGroups.jar \
        INPUT=$input/$sample.sorted.bam \
        OUTPUT=$input/$sample.sorted.rg.bam \
        PL=$platform SM=$sample CN=$center ID=$sample PU=$sample LB=$GenomeBuild \
        TMP_DIR=$input \
        VALIDATION_STRINGENCY=SILENT
        if [ -s $input/$sample.sorted.rg.bam ]
        then
            mv $input/$sample.sorted.rg.bam $input/$sample.sorted.bam
        else
            echo "ERROR: merge_align.bam.sh File $input/$sample.sorted.rg.bam not generated"
            exit 1
        fi
    fi
	if [ $reorder == "YES" ]
	then
		$java/java -Xmx6g -Xms512m \
		-jar $picard/ReorderSam.jar \
		I=$input/$sample.sorted.bam \
		O=$input/$sample.sorted.tmp.bam \
		R=$ref_path \
		TMP_DIR=$input \
		VALIDATION_STRINGENCY=SILENT
		
		if [ -s $input/$sample.sorted.tmp.bam ]
		then
			mv $input/$sample.sorted.tmp.bam $input/$sample.sorted.bam
		else
			echo "ERROR: reordering $sample failed $input/$sample.sorted.tmp.bam "
			exit 1
		fi    
	fi
	
    echo `date`
fi	
	