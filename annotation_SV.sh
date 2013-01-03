#!/bin/bash
    
########################################################
###### 	SV ANNOTATION FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			annotation.SV.sh
######		Date:				11/09/2011
######		Summary:			Annotates CREST and BREAKDANCER outputs
######		Input 
######		$1	=	structural directory
######		$2	=	/path/to/run_info.txt
########################################################

if [ $# != 4 ]
then
    echo -e "script to get annotated report for SNV reported for a sample\nUsage: </path/to/output dir> </path/to/run_info.txt>c < output folder> <smaple/group name>";
else
    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    report_dir=$3
	sample=$4
	group=$sample
    
########################################################	
######		Reading run_info.txt and assigning to variables

    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	Rsoft=$( cat $tool_info | grep -w '^R_SOFT' | cut -d '=' -f2 )
	export PATH=$Rsoft:$PATH
	somatic_calling=$( cat $tool_info | grep -w '^SOMATIC_CALLING' | cut -d '=' -f2 )
##############################################################		
        
    SV_dir=$output_dir/Reports_per_Sample/
    
    if [[ $multi_sample != "YES"  || $somatic_calling == "NO" ]]
    then
        echo "Single sample"
        ### preparing breakdancer files
        break=$SV_dir/SV
        crest=$SV_dir/SV
        
        if [ ! -s $break/$sample.break.vcf ]
        then
			$script_path/email.sh $break/$sample.break.vcf "not exist" summaryze_struct_single.sh $run_info
			touch $break/$sample.break.vcf.fix.log
			$script_path/wait.sh $break/$sample.break.vcf.fix.log
        fi
        
        cat $crest/$sample.filter.crest.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $crest/$sample.crest.tmp
        cat $crest/$sample.crest.tmp | awk 'gsub($5, "crest_"$5,$5)1' | tr " " "\t" > $crest/$sample.crest.txt
        cat $break/$sample.break.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $break/$sample.break.tmp
        cat $break/$sample.break.tmp | awk 'gsub($5, "breakdancer_"$5,$5)1' | tr " " "\t" > $break/$sample.breakdancer.txt
        
        if [ ! -s $break/$sample.breakdancer.txt ]
        then
            $script_path/errorlog.sh $break/$sample.breakdancer.txt annotation_SV.sh WARNING "no SV "
            touch $report_dir/$sample.SV.annotated.txt
			rm $crest/$sample.crest.tmp $break/$sample.break.tmp $break/$sample.breakdancer.txt $crest/$sample.crest.txt 
        else
            ### concatenating & formatting crest and breakdancer files
            touch $SV_dir/$sample.SV.tmp
            # cat $break/$sample.breakdancer.txt $crest/$sample.crest.bed > $SV_dir/$sample.tmp
            cat $break/$sample.breakdancer.txt $crest/$sample.crest.txt > $SV_dir/$sample.tmp
            cat -n $SV_dir/$sample.tmp > $SV_dir/$sample.tmp1
            cat $SV_dir/$sample.tmp1 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' >> $SV_dir/$sample.SV.tmp
            cat $SV_dir/$sample.SV.tmp |  awk '{if ($3>10000) print $0}' > $SV_dir/$sample.SV.tmp1
            cat $SV_dir/$sample.SV.tmp1 | awk '{print $2"\t"$3"\t"($3+10000)"\t"$1}' > $SV_dir/$sample.SV.leftend.bed
            cat $SV_dir/$sample.SV.tmp1 | awk '{print $5"\t"($6-10000)"\t"$6"\t"$4}' > $SV_dir/$sample.SV.rightend.bed                            
            if [ `cat $SV_dir/$sample.SV.leftend.bed | wc -l` -gt 0 ]
			then
				$bedtools/intersectBed -b $SV_dir/$sample.SV.leftend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.leftend.intersect.bed
            else
				touch $SV_dir/$sample.SV.leftend.intersect.bed
			fi
			if [ `cat $SV_dir/$sample.SV.rightend.bed | wc -l` -gt 0 ]
			then
				$bedtools/intersectBed -b $SV_dir/$sample.SV.rightend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.rightend.intersect.bed
            else
				touch $SV_dir/$sample.SV.rightend.intersect.bed
			fi
			
			$script_path/GeneAnnotation.SV.pl $SV_dir/$sample.SV.tmp $SV_dir/$sample.SV.leftend.intersect.bed $SV_dir/$sample.SV.rightend.intersect.bed $SV_dir/$sample.SV.leftend.annotation.tmp $SV_dir/$sample.SV.rightend.annotation.tmp
            cat $SV_dir/$sample.SV.leftend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.leftend.annotation.txt
            cat $SV_dir/$sample.SV.rightend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.rightend.annotation.txt
            cat $SV_dir/$sample.SV.rightend.annotation.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > $SV_dir/$sample.SV.rightend.annotation.final.txt
            paste $SV_dir/$sample.SV.leftend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.final.txt > $SV_dir/$sample.SV.annotated.tmp
            $Rsoft/Rscript $script_path/format_SV_annotation.r $SV_dir/$sample.SV.annotated.tmp $SV_dir/$sample.SV.annot.txt
            cat $SV_dir/$sample.SV.annot.txt | sed '/NOGENE/s//NA/g' > $SV_dir/$sample.SV.annot.tmp1
            cat $SV_dir/$sample.SV.annot.tmp1 | cut -f2,3,4,6,7,8,9,10 > $SV_dir/$sample.SV.txt
            cat $SV_dir/$sample.SV.txt | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8"\t"$3"_"$7}' > $SV_dir/$sample.SV.temp
            cat $SV_dir/$sample.SV.temp | awk '$0 !~ "NA_NA"' > $SV_dir/$sample.SV.temp1
            #  generating per sample SV annotation files
            touch $report_dir/$sample.SV.annotated.txt
            if [ ! -f $SV_dir/$sample.SV.temp1 ]
            then
                $script_path/errorlog.sh $SV_dir/$sample.SV.temp1 annotation_SV.sh ERROR "failed to create"
                exit 1
            fi    
            echo -e "ChrA\tPosA\tChrB\tPosB\tSV_Type\tSpanning_Reads\tGeneA_GeneB" >> $report_dir/$sample.SV.annotated.txt
            cat $SV_dir/$sample.SV.temp1 >> $report_dir/$sample.SV.annotated.txt
            
            if [ -s $report_dir/$sample.SV.annotated.txt ]
            then
                # removing intermediate files
                rm $SV_dir/$sample.tmp $SV_dir/$sample.tmp1
                rm $break/$sample.breakdancer.txt $crest/$sample.crest.txt $break/$sample.break.tmp $break/$sample.crest.tmp
				rm $SV_dir/$sample.SV.leftend.bed $SV_dir/$sample.SV.rightend.bed $SV_dir/$sample.SV.tmp $SV_dir/$sample.SV.temp $SV_dir/$sample.SV.tmp1 $SV_dir/$sample.SV.temp1
                rm $SV_dir/$sample.SV.leftend.intersect.bed $SV_dir/$sample.SV.rightend.intersect.bed
                rm $SV_dir/$sample.SV.leftend.annotation.tmp $SV_dir/$sample.SV.rightend.annotation.tmp
                rm $SV_dir/$sample.SV.leftend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.final.txt
                rm $SV_dir/$sample.SV.annotated.tmp $SV_dir/$sample.SV.annot.txt $SV_dir/$sample.SV.annot.tmp1 $SV_dir/$sample.SV.txt
            else
                $script_path/errorlog.sh $report_dir/$sample.SV.annotated.txt annotation_SV.sh ERROR "failed to create"
				exit 1;
            fi    
        fi
    else
        echo "Multi sample"
		samples=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2)
		let num_tumor=`echo $samples|tr " " "\n"|wc -l`-1
		tumor_list=`echo $samples | tr " " "\n" | tail -$num_tumor`
        break=$SV_dir/SV
        crest=$SV_dir/SV
		for tumor in $tumor_list
		do
			if [ ! -s $break/$group.$tumor.somatic.break.vcf ]
			then
				$script_path/email.sh $break/$group.$tumor.somatic.break.vcf "not exist" summaryze_struct_group.sh $run_info
				touch $break/$group.$tumor.somatic.break.vcf.fix.log
				$script_path/wait.sh $break/$group.$tumor.somatic.break.vcf.fix.log
			fi
			cat $crest/$group.$tumor.somatic.filter.crest.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $crest/$group.$tumor.crest.tmp
			cat $crest/$group.$tumor.crest.tmp | awk 'gsub($5, "crest_"$5,$5)1' | tr " " "\t" > $crest/$group.$tumor.crest.txt
			cat $break/$group.$tumor.somatic.break.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $break/$group.$tumor.break.tmp
			cat $break/$group.$tumor.break.tmp | awk 'gsub($5, "breakdancer_"$5,$5)1' | tr " " "\t" > $break/$group.$tumor.breakdancer.txt
							
			if [ ! -s $break/$group.$tumor.breakdancer.txt ]
            then
                $script_path/errorlog.sh $break/$group.$tumor.breakdancer.txt annotation_SV.sh WARNING "not SV"
                touch $report_dir/$group.$tumor.SV.annotated.txt
				rm $crest/$group.$tumor.crest.tmp $crest/$group.$tumor.crest.txt $break/$group.$tumor.break.tmp $break/$group.$tumor.breakdancer.txt
            else
				### concatenating & formatting crest and breakdancer files
                touch $SV_dir/$group.$tumor.SV.tmp
                cat $break/$group.$tumor.breakdancer.txt $crest/$group.$tumor.crest.txt > $SV_dir/$group.$tumor.tmp
                cat -n $SV_dir/$group.$tumor.tmp > $SV_dir/$group.$tumor.tmp1
                cat $SV_dir/$group.$tumor.tmp1 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' >> $SV_dir/$group.$tumor.SV.tmp
                cat $SV_dir/$group.$tumor.SV.tmp |  awk '{if ($3>10000) print $0}' > $SV_dir/$group.$tumor.SV.tmp1
                cat $SV_dir/$group.$tumor.SV.tmp1 | awk '{print $2"\t"$3"\t"($3+10000)"\t"$1}' > $SV_dir/$group.$tumor.SV.leftend.bed
                cat $SV_dir/$group.$tumor.SV.tmp1 | awk '{print $5"\t"($6-10000)"\t"$6"\t"$4}' > $SV_dir/$group.$tumor.SV.rightend.bed                             
                if [ `cat $SV_dir/$group.$tumor.SV.leftend.bed | wc -l` -gt 0 ]
				then
					$bedtools/intersectBed -b $SV_dir/$group.$tumor.SV.leftend.bed -a $master_gene_file -wb > $SV_dir/$group.$tumor.SV.leftend.intersect.bed
                else
					touch $SV_dir/$group.$tumor.SV.leftend.intersect.bed
				fi
				if [ `cat $SV_dir/$group.$tumor.SV.rightend.bed | wc -l` -gt 0 ]
				then
					$bedtools/intersectBed -b $SV_dir/$group.$tumor.SV.rightend.bed -a $master_gene_file -wb > $SV_dir/$group.$tumor.SV.rightend.intersect.bed
                else
					touch $SV_dir/$group.$tumor.SV.rightend.intersect.bed	
				fi
				$script_path/GeneAnnotation.SV.pl $SV_dir/$group.$tumor.SV.tmp $SV_dir/$group.$tumor.SV.leftend.intersect.bed $SV_dir/$group.$tumor.SV.rightend.intersect.bed $SV_dir/$group.$tumor.SV.leftend.annotation.tmp $SV_dir/$group.$tumor.SV.rightend.annotation.tmp
                cat $SV_dir/$group.$tumor.SV.leftend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$group.$tumor.SV.leftend.annotation.txt
                cat $SV_dir/$group.$tumor.SV.rightend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$group.$tumor.SV.rightend.annotation.txt
                cat $SV_dir/$group.$tumor.SV.rightend.annotation.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > $SV_dir/$group.$tumor.SV.rightend.annotation.final.txt
                paste $SV_dir/$group.$tumor.SV.leftend.annotation.txt $SV_dir/$group.$tumor.SV.rightend.annotation.final.txt > $SV_dir/$group.$tumor.SV.annotated.tmp
                $Rsoft/Rscript $script_path/format_SV_annotation.r $SV_dir/$group.$tumor.SV.annotated.tmp $SV_dir/$group.$tumor.SV.annot.txt
                cat $SV_dir/$group.$tumor.SV.annot.txt | sed '/NOGENE/s//NA/g' > $SV_dir/$group.$tumor.SV.annot.tmp1
                cat $SV_dir/$group.$tumor.SV.annot.tmp1 | cut -f2,3,4,6,7,8,9,10 > $SV_dir/$group.$tumor.SV.txt
                cat $SV_dir/$group.$tumor.SV.txt | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8"\t"$3"_"$7}' > $SV_dir/$group.$tumor.SV.temp
                cat $SV_dir/$group.$tumor.SV.temp | awk '$0 !~ "NA_NA"' > $SV_dir/$group.$tumor.SV.temp1
                #  generating per sample SV annotation files
                touch $report_dir/$group.$tumor.SV.annotated.txt
                if [ ! -f $SV_dir/$group.$tumor.SV.temp1 ]
                then
                    $script_path/errorlog.sh $SV_dir/$group.$tumor.SV.temp1 annotation_SV.sh WARNING "failed to create"
                    exit 1
                fi    
                echo -e "ChrA\tPosA\tChrB\tPosB\tSV_Type\tSpanning_Reads\tGeneA_GeneB" >> $report_dir/$group.$tumor.SV.annotated.txt
                cat $SV_dir/$group.$tumor.SV.temp1 >> $report_dir/$group.$tumor.SV.annotated.txt
            
                if [ -s $report_dir/$group.$tumor.SV.annotated.txt ]
                then
                    # removing intermediate files
                    rm $SV_dir/$group.$tumor.tmp $SV_dir/$group.$tumor.tmp1
                    rm $break/$group.$tumor.breakdancer.txt $crest/$group.$tumor.crest.txt $break/$group.$tumor.break.tmp $break/$group.$tumor.crest.tmp
                    rm $SV_dir/$group.$tumor.SV.leftend.bed $SV_dir/$group.$tumor.SV.rightend.bed $SV_dir/$group.$tumor.SV.tmp $SV_dir/$group.$tumor.SV.temp $SV_dir/$group.$tumor.SV.tmp1 $SV_dir/$group.$tumor.SV.temp1
                    rm $SV_dir/$group.$tumor.SV.leftend.intersect.bed $SV_dir/$group.$tumor.SV.rightend.intersect.bed
                    rm $SV_dir/$group.$tumor.SV.leftend.annotation.tmp $SV_dir/$group.$tumor.SV.rightend.annotation.tmp
                    rm $SV_dir/$group.$tumor.SV.leftend.annotation.txt $SV_dir/$group.$tumor.SV.rightend.annotation.txt $SV_dir/$group.$tumor.SV.rightend.annotation.final.txt
                    rm $SV_dir/$group.$tumor.SV.annotated.tmp $SV_dir/$group.$tumor.SV.annot.txt $SV_dir/$group.$tumor.SV.annot.tmp1 $SV_dir/$group.$tumor.SV.txt
                else
                    $script_path/errorlog.sh $report_dir/$group.$tumor.SV.annotated.txt annotation_SV.sh WARNING "failed to create"
					exit 1;
                fi  
            fi
        done
	fi
	echo `date`
fi
