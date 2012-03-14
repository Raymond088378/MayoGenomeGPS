if [ $# != 7 ]
then
    echo -e "Usage: wrapper to validate bam file\n samplecheckBAM.sh <input dir><input bam><outputdir><run_info><sample><1 or 0 if bam is per chr><which _chr>";
else	
    set -x
    echo `date`
    input=$1
    bam=$2
    output=$3
    run_info=$4
    sample=$5
    chopped=$6
    chr=$7
    
    if [ -d $output/temp ]
    then
        echo "already there"
    else
        mkdir $output/temp
    fi
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )

    
    check=`[ -f $input/$bam.bai ] && echo "1" || echo "0"`
    if [ $check -eq 0 ]
    then
        ln -s $input/$bam $output/$bam.$chr.bam
        $samtools/samtools index $output/$bam.$chr.bam
    else
        ln -s $input/$bam $output/$bam.$chr.bam
        ln -s $input/$bam.bai $output/$bam.$chr.bam.bai	
    fi	
            
    ## extracting the BAM for specific chromosome
    if [ $chopped == 0 ]
    then
        $samtools/samtools view -b  $output/$bam.$chr.bam chr${chr} > $output/$sample.chr${chr}.bam
        $samtools/samtools index $output/$sample.chr${chr}.bam
    else
        ln -s  $output/$bam.$chr.bam $output/$sample.chr${chr}.bam
        $samtools/samtools index $output/$sample.chr${chr}.bam
    fi
    
    ## check if BAM is sorted
    SORT_FLAG=`perl $script_path/checkBAMsorted.pl -i $output/$sample.chr${chr}.bam -s $samtools`
    if [ $SORT_FLAG == 1 ]
    then
        ln -s $output/$sample.chr${chr}.bam $output/$sample.chr${chr}-sorted.bam
    else
        $java/java -Xmx6g -Xms512m \
        -jar $picard/SortSam.jar \
        INPUT=$output/$sample.chr${chr}.bam \
        OUTPUT=$output/$sample.chr${chr}-sorted.bam \
        SO=coordinate \
        MAX_RECORDS_IN_RAM=5000000 \
        TMP_DIR=$output/temp \
        VALIDATION_STRINGENCY=SILENT
    fi
    
    ## check if read group and platform is availbale in BAM
    RG_FLAG=`perl $script_path/checkBAMreadGroup.pl -i $output/$sample.chr${chr}-sorted.bam -s $samtools`
    if [ $RG_FLAG == 0 ]
    then
        center=$( cat $run_info | grep -w '^CENTER' | cut -d '=' -f2 )
        platform=$( cat $run_info | grep -w '^PLATFORM' | cut -d '=' -f2 )
        GenomeBuild=$( cat $run_info | grep -w '^GENOMEBUILD' | cut -d '=' -f2 )
        Aligner=$( cat $run_info | grep -w '^ALIGNER' | cut -d '=' -f2 )
        Aligner=`echo "$Aligner" | tr "[a-z]" "[A-Z]"`
        $java/java -Xmx6g -Xms512m \
        -jar $picard/AddOrReplaceReadGroups.jar \
        INPUT=$output/$sample.chr${chr}-sorted.bam \
        OUTPUT=$output/$sample.chr${chr}-sorted.bam.rg.bam \
        PL=$platform SM=$sample PU=$sample CN=$center ID=$sample LB=$GenomeBuild \
        TMP_DIR=$output/temp \
        VALIDATION_STRINGENCY=SILENT
        if [ -s $output/$sample.chr${chr}-sorted.bam.rg.bam ]
        then
            mv $output/$sample.chr${chr}-sorted.bam.rg.bam $output/$sample.chr${chr}-sorted.bam
        else
            echo "ERROR: $output/$sample.chr${chr}-sorted.bam.rg.bam failed to generate"
            exit 1 
        fi    
    fi	
    
    ### after this point BAM is good to go with the GATk tool (any module)
    $samtools/samtools index $output/$sample.chr${chr}-sorted.bam
    echo `date`
fi	
	