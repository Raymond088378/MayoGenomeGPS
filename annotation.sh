#!/bin/sh
	
##	INFO
## 	This module is used for annotating variants

########################### 
#       $1      =       OUtput directroy
#		$2		=		run information
###########################

if [ $# != 2 ];
then
    echo -e "Usage: SCRIPT to annotate variants \n <output dir> <run info>"
else	
    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)	
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    analysis=$( cat $run_info | grep -w '^ANALYSIS' | cut -d '=' -f2)
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    dbsnp_rsids=$( cat $tool_info | grep -w '^dbSNP_SNV_rsIDs' | cut -d '=' -f2 )
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2 )
    variant_type=$( cat $run_info | grep -w '^VARIANT_TYPE' | cut -d '=' -f2)
    chrs=$( cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2)
    analysis=`echo "$analysis" | tr "[A-Z]" "[a-z]"`
    variant_type=`echo "$variant_type" | tr "[a-z]" "[A-Z]"`
    
    ## creating the folder structure
    mkdir $output_dir/annotation
    output_annot=$output_dir/annotation
    mkdir $output_annot/SIFT
    sift=$output_annot/SIFT
    mkdir $output_annot/SSEQ
    sseq=$output_annot/SSEQ
    output_OnTarget=$output_dir/OnTarget
    #extracting samples and chr
    sampleNames=$( echo $samples | tr ":" "\n" )
    chrIndexes=$( echo $chrs | tr ":" "\n" )
    i=1
    for sample in $sampleNames
    do
            sampleArray[$i]=$sample
            let i=i+1
    done
    i=1
    for chr in $chrIndexes
    do
            chrArray[$i]=$chr
            let i=i+1
    done
    array_jobs=${#chrArray[@]}
    array_sample_jobs=${#sampleArray[@]}
    #path=$script_path/annotation
    path=$script_path
    arg="-V -wd $output_dir/logs -q $queue -m a -M $email -l h_stack=10M"
    
    if [ $analysis != "annotation" ]
    then
        job_ids=$(cat $output_dir/job_ids/VARIANTS | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ",")
        #Call sift and sseq per sample for SNVs
        SIFT=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sift.sh $sift $output_OnTarget $run_info` 
        SSEQ=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sseq.sh $sseq $output_OnTarget $email $run_info`	
    #this module take care that the user wants to annotate the SNV or INDEL or BOTH
    elif [ $analysis == "annotation" ]
    then
        input=$( cat $run_info | grep -w INPUT_DIR | cut -d '=' -f2)
        mkdir $output_dir/OnTarget
        output_OnTarget=$output_dir/OnTarget
        touch $output_dir/job_ids/all_annot_jobs
        if [ $variant_type == "BOTH" ]
        then
                REFORMAT=`qsub $arg -t 1-$array_sample_jobs -l h_vmem=4G $path/reformat.VARIANTs.sh $output_OnTarget $input $sample_info $run_info 2`
                job_ids=`echo $REFORMAT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
                SIFT=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sift.sh $sift $output_OnTarget $run_info`
                SSEQ=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sseq.sh $sseq $output_OnTarget $email $run_info`		
        elif [ $variant_type == "SNV" ]
        then
            REFORMAT=`qsub $arg -t 1-$array_sample_jobs -l h_vmem=4G $path/reformat.VARIANTs.sh $output_OnTarget $input $sample_info $run_info 1`
            job_ids=`echo $REFORMAT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
            SIFT=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sift.sh $sift $output_OnTarget $run_info`
            SSEQ=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sseq_SNV.sh $sseq $output_OnTarget $email $run_info`	
        elif [ $variant_type == "INDEL" ]
        then
            REFORMAT=`qsub $arg -t 1-$array_sample_jobs -l h_vmem=4G $path/reformat.VARIANTs.sh $output_OnTarget $input $sample_info $run_info 1`
            job_ids=`echo $REFORMAT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
            SSEQ=`qsub $arg -hold_jid $job_ids -t 1-$array_sample_jobs -l h_vmem=4G $path/sseq_INDEL.sh $sseq $output_OnTarget $email $run_info`
        fi
    fi	
    ## merge the sift ids
    job_ids_sift=`echo $SIFT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
    job_ids_sseq=`echo $SSEQ | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
    #SIFTID=`qsub $arg -hold_jid $job_ids_sseq,$job_ids_sift -l h_vmem=4G $path/merge.siftid.sh $sift` 
    #job_ids_SIFTID=`echo $SIFTID | cut -d ' ' -f3`
    
    ## creating folder structure	
    mkdir $output_dir/TempReports
    TempReports=$output_dir/TempReports
    mkdir $output_dir/Reports_per_Sample
    mkdir $output_dir/VariantDatabase
    job_dir=$output_dir/job_ids
    ## making reports per sample
    for i in $(seq 1 ${#sampleArray[@]})
    do 
        sample=${sampleArray[$i]}	
        qsub $arg -hold_jid $job_ids_sseq,$job_ids_sift -t 1-$array_jobs -l h_vmem=8G $path/per.sample.reports.per.chr.sh $run_info $sample $TempReports $output_OnTarget $sift $sseq $output_dir >> $job_dir/ANNOTATION.forsamples
    done	
    job_ids_per_sample=$( cat $job_dir/ANNOTATION.forsamples | cut -d ' ' -f3 | cut -d '.' -f1 | tr "\n" ",")
    rm $job_dir/ANNOTATION.forsamples
    ## to merge reports
    MERGE=`qsub $arg -hold_jid $job_ids_per_sample -t 1-$array_sample_jobs -l h_vmem=4G $path/merge.per.sample.report.sh $output_dir $TempReports $run_info`
    job_ids_merge=`echo $MERGE | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
    
    ## annoatate all the per sample reports
    qsub $arg -hold_jid $job_ids_merge -l h_vmem=8G $path/annotate.per.sample.sh $output_dir $run_info >> $output_dir/job_ids/ANNOTATION 
    
    ## file for variant database
    job_ids_ann=$( cat $output_dir/job_ids/ANNOTATION | cut -d ' ' -f3  | tr "\n" "," )
    qsub $arg -hold_jid $job_ids_merge -t 1-$array_sample_jobs -l h_vmem=4G $path/variantDatabase.sh $TempReports $output_dir/VariantDatabase $run_info
    
    ## to gerneate merge reports
    if [ $analysis != "annotation" ]
    then
        mkdir $output_dir/Reports
        REPORT=`qsub $arg -hold_jid $job_ids_sseq,$job_ids_sift -t 1-$array_jobs -l h_vmem=8G $path/reports.per.chr.sh $sift $sseq $TempReports $run_info $output_dir/OnTarget` 
        job_ids_report=`echo $REPORT | cut -d ' ' -f3 | tr "\n" "," | sed -e "s/\..*,//g"`
        ## annotate merged report
        qsub $arg -hold_jid $job_ids_report -l h_vmem=8G $path/merge.merged.report.sh $output_dir $TempReports $run_info >> $output_dir/job_ids/ANNOTATION
        ### variant distance
        qsub $arg -hold_jid $job_ids_report -l h_vmem=4G $path/variant.distance.sh $TempReports $output_dir $run_info
    fi		
    echo `date`
fi	
	
## end of annoatation module script	
	
	
		
		
			
		
	
