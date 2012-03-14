#!/bin/sh
    
########################################################
###### 	SV ANNOTATION FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			annotation.SV.sh
######		Date:				11/09/2011
######		Summary:			Annotates CREST and BREAKDANCER outputs
######		Input 
######		$1	=	structural directory
######		$2	=	/path/to/run_info.txt
########################################################

if [ $# != 3 ]
then
    echo "\nUsage: </path/to/output dir> </path/to/run_info.txt>c < output folder>";
else
    set -x
    echo `date`
    output_dir=$1
    run_info=$2
    report_dir=$3
    
########################################################	
######		Reading run_info.txt and assigning to variables

    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    sample=$(cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)  
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
	PATH=$bedtools/:$PATH
##############################################################		
        
    SV_dir=$output_dir/Reports_per_Sample/
    
    if [ $multi_sample != "YES" ]
    then
        echo "Single sample"
        ### preparing breakdancer files
        break=$SV_dir/SV
        crest=$SV_dir/SV
        
        if [ ! -s $break/$sample.break.vcf ]
        then
            echo "ERROR :anotation.SV.sh $break/$sample.break.vcf is empty"
            exit 1
        fi
        
        cat $crest/$sample.filter.crest.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $crest/$sample.crest.tmp
        cat $crest/$sample.crest.tmp | awk 'gsub($5, "crest_"$5,$5)1' | tr " " "\t" > $crest/$sample.crest.txt
        cat $break/$sample.break.vcf | awk '$0 !~ /#/' | tr ";" "\t" | tr "=" "\t" | awk '{ print $1"\t"$2"\t"$1"\t"$12"\t"$10"\t"$16}' > $break/$sample.break.tmp
        cat $break/$sample.break.tmp | awk 'gsub($5, "breakdancer_"$5,$5)1' | tr " " "\t" > $break/$sample.breakdancer.txt
                        
        ### concatenating & formatting crest and breakdancer files
        touch $SV_dir/$sample.SV.tmp
        # cat $break/$sample.breakdancer.txt $crest/$sample.crest.bed > $SV_dir/$sample.tmp
        cat $break/$sample.breakdancer.txt $crest/$sample.crest.txt > $SV_dir/$sample.tmp

        cat -n $SV_dir/$sample.tmp > $SV_dir/$sample.tmp1
        cat $SV_dir/$sample.tmp1 | awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' >> $SV_dir/$sample.SV.tmp
        cat $SV_dir/$sample.SV.tmp |  awk '{if ($3>10000) print $0}' > $SV_dir/$sample.SV.tmp1
        cat $SV_dir/$sample.SV.tmp1 | awk '{print $2"\t"$3"\t"($3+10000)"\t"$1}' > $SV_dir/$sample.SV.leftend.bed
        cat $SV_dir/$sample.SV.tmp1 | awk '{print $5"\t"($6-10000)"\t"$6"\t"$4}' > $SV_dir/$sample.SV.rightend.bed
                                        
        $bedtools/intersectBed -b $SV_dir/$sample.SV.leftend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.leftend.intersect.bed
        $bedtools/intersectBed -b $SV_dir/$sample.SV.rightend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.rightend.intersect.bed
        perl $script_path/GeneAnnotation.SV.pl $SV_dir/$sample.SV.tmp $SV_dir/$sample.SV.leftend.intersect.bed $SV_dir/$sample.SV.rightend.intersect.bed $SV_dir/$sample.SV.leftend.annotation.tmp $SV_dir/$sample.SV.rightend.annotation.tmp
        cat $SV_dir/$sample.SV.leftend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.leftend.annotation.txt
        cat $SV_dir/$sample.SV.rightend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.rightend.annotation.txt
        cat $SV_dir/$sample.SV.rightend.annotation.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > $SV_dir/$sample.SV.rightend.annotation.final.txt
        paste $SV_dir/$sample.SV.leftend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.final.txt > $SV_dir/$sample.SV.annotated.tmp
        Rscript $script_path/format_SV_annotation.r $SV_dir/$sample.SV.annotated.tmp $SV_dir/$sample.SV.annot.txt
        cat $SV_dir/$sample.SV.annot.txt | sed '/NOGENE/s//NA/g' > $SV_dir/$sample.SV.annot.tmp1
        cat $SV_dir/$sample.SV.annot.tmp1 | cut -f2,3,4,6,7,8,9,10 > $SV_dir/$sample.SV.txt
        cat $SV_dir/$sample.SV.txt | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8"\t"$3"_"$7}' > $SV_dir/$sample.SV.temp
        cat $SV_dir/$sample.SV.temp | awk '$0 !~ "NA_NA"' > $SV_dir/$sample.SV.temp1
        #  generating per sample SV annotation files
        touch $report_dir/$sample.SV.annotated.txt
        if [ ! -s $SV_dir/$sample.SV.temp1 ]
        then
            echo "ERROR : annotation.SV.sh file failed for $sample"
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
            echo "ERROR : annotation.SV.sh file $report_dir/$sample.SV.annotated.txt is empty"
        fi    
    else
        
                echo "Multi sample"
                
                $script_path/breakdancer_tumor-normal.sh $output_dir $run_info
                
                for group in `echo $groups | tr ":" " "`
                do
                        sampleNames=$( cat $sample_info | grep -w "^$group" | cut -d '=' -f2)

                        i=1
                        for sample in $sampleNames
                        do
                                sampleArray[$i]=$sample
                                let i=i+1
                        done
                        for i in $(seq 1 ${#sampleArray[@]})
                        do
                                if [ $(($i % 2)) -eq 0 ]
                                then
                                        sample=${sampleArray[$i]}
                                        crest=$SV_dir/crest/$group
                                        touch $crest/$sample.predSV.txt
                                        for j in `seq 1 22` X Y M
                                        do
                                                cat $crest/$sample.$j.predSV.txt >> $crest/$sample.predSV.txt
                                        done
                                        
                                        cat $crest/$sample.predSV.txt | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"($4+$8)"\t""crest_"$9}' > $crest/$sample.crest.bed

                                        ### preparing breakdancer files
                                        break=$SV_dir/break/$group/$sample
                                        
                                        for j in `seq 1 22` X Y M
                                        do
                                                cat $break/$sample.$j.break >> $break/$sample.break
                                        done
                                        
                                        cat $break/$sample.break | cut -f1,2,3,4,5,6,7,10 > $break/$sample.tmp
                                        cat $break/$sample.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' > $break/$sample.tmp1
                                        sed -e '/#/d' $break/$sample.tmp1 | awk 'gsub($8, "breakdancer_"$8,$8)1' > $break/$sample.tmp2
                                        cat $break/$sample.tmp2 | tr " " "\t" > $break/$sample.breakdancer.txt
                        
                                        ### concatenating & formatting crest and breakdancer files
                                        touch $SV_dir/$sample.SV.tmp
                                        cat $break/$sample.breakdancer.txt $crest/$sample.crest.bed > $SV_dir/$sample.tmp
                                        cat -n $SV_dir/$sample.tmp > $SV_dir/$sample.tmp1
                                        cat $SV_dir/$sample.tmp1 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >> $SV_dir/$sample.SV.tmp
                                        cat $SV_dir/$sample.SV.tmp |  awk '{if ($3>10000) print $0}' > $SV_dir/$sample.SV.tmp1
                                        cat $SV_dir/$sample.SV.tmp1 | awk '{print $2"\t"$3"\t"($3+10000)"\t"$1}' > $SV_dir/$sample.SV.leftend.bed
                                        cat $SV_dir/$sample.SV.tmp1 | awk '{print $6"\t"($7-10000)"\t"$7"\t"$5}' > $SV_dir/$sample.SV.rightend.bed
                                        $bedtools/intersectBed -b $SV_dir/$sample.SV.leftend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.leftend.intersect.bed
                                        $bedtools/intersectBed -b $SV_dir/$sample.SV.rightend.bed -a $master_gene_file -wb > $SV_dir/$sample.SV.rightend.intersect.bed
                                        perl $script_path/GeneAnnotation.SV.pl $SV_dir/$sample.SV.tmp $SV_dir/$sample.SV.leftend.intersect.bed $SV_dir/$sample.SV.rightend.intersect.bed $SV_dir/$sample.SV.leftend.annotation.tmp $SV_dir/$sample.SV.rightend.annotation.tmp
                                        cat $SV_dir/$sample.SV.leftend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.leftend.annotation.txt
                                        cat $SV_dir/$sample.SV.rightend.annotation.tmp | tr "^" "\t" | sort -n > $SV_dir/$sample.SV.rightend.annotation.txt
                                        cat $SV_dir/$sample.SV.rightend.annotation.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$5"\t"$6}' > $SV_dir/$sample.SV.rightend.annotation.final.txt
                                        paste $SV_dir/$sample.SV.leftend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.final.txt > $SV_dir/$sample.SV.annotated.tmp
                                        Rscript $script_path/format_SV_annotation.r $SV_dir/$sample.SV.annotated.tmp $SV_dir/$sample.SV.annot.txt
                                        cat $SV_dir/$sample.SV.annot.txt | sed '/NOGENE/s//NA/g' > $SV_dir/$sample.SV.annot.tmp1
                                        cat $SV_dir/$sample.SV.annot.tmp1 | cut -f2,3,4,5,7,8,9,10,11,12 > $SV_dir/$sample.SV.txt
                                        cat $SV_dir/$sample.SV.txt | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10"\t"$4"_"$8}' > $SV_dir/$sample.SV.tmp

                                        ###  generating per sample SV annotation files
                                        touch $report_dir/$sample.SV.annotated.txt
                                        echo -e "ChrA\tPosA\tOrientationA/StrandA\tChrB\tPosB\tOrientationB/StrandB\tSpanning_Reads\tSV_Type\tGeneA_GeneB" >> $report_dir/$sample.SV.annotated.txt
                                        cat $SV_dir/$sample.SV.tmp >> $report_dir/$sample.SV.annotated.txt
                                        
                                        ### removing intermediate files
                                        rm $break/$sample.tmp $break/$sample.tmp1 $break/$sample.tmp2
                                        rm $SV_dir/$sample.tmp $SV_dir/$sample.tmp1 $SV_dir/$sample.SV.leftend.bed $SV_dir/$sample.SV.rightend.bed
                                        rm $SV_dir/$sample.SV.tmp $SV_dir/$sample.SV.leftend.intersect.bed $SV_dir/$sample.SV.rightend.intersect.bed
                                        rm $SV_dir/$sample.SV.leftend.annotation.tmp $SV_dir/$sample.SV.rightend.annotation.tmp $SV_dir/$sample.SV.tmp1
                                        rm $SV_dir/$sample.SV.leftend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.txt $SV_dir/$sample.SV.rightend.annotation.final.txt
                                        rm $SV_dir/$sample.SV.annotated.tmp $SV_dir/$sample.SV.annot.txt $SV_dir/$sample.SV.annot.tmp1 $SV_dir/$sample.SV.txt $SV_dir/$sample.SV.tmp
                                else
                                        echo "Normal sample....skipping"
                                fi
                        done
                done
        fi
fi
