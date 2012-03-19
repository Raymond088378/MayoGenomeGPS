#!/bin/sh
	
########################################################
###### 	SNV ANNOTATION FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE

######		Program:			annotation.SNV.sh
######		Date:				11/09/2011
######		Summary:			Annotates GATK SNV and INDEL outputs
######		Input 
######		$1	=	structural directory
######		$2	=	/path/to/run_info.txt
########################################################

if [ $# != 4 ]
then
    echo -e "\nUsage: </path/to/output directory for variants> </path/to/output directory for OnTarget> <sample name> </path/to/run_info.txt>";
else
    set -x
    echo `date`
    variants=$1
    OnTarget=$2
    sample=$3
    run_info=$4
    
#SGE_TASK_ID=1
########################################################	
######		Reading run_info.txt and assigning to variables

    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
    master_gene_file=$( cat $tool_info | grep -w '^MASTER_GENE_FILE' | cut -d '=' -f2 )
    email=$( cat $run_info | grep -w '^EMAIL' | cut -d '=' -f2)
    queue=$( cat $run_info | grep -w '^QUEUE' | cut -d '=' -f2)
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    chrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n")
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    script_path=$script_path
    tool=`echo "$tool" | tr "[A-Z]" "[a-z]"`
	PATH=$bedtools/:$PATH
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
##############################################################		
    
    
    if [ $multi_sample != "YES" ]
    then
        echo "Single sample"
        input=$variants/$sample
        ### parse the file to generate seperate SNV and INDEL vcf
        if [ ! -s $input/$sample.variants.chr$chr.filter.vcf ]
        then
            echo "ERROR: $input/$sample.variants.chr$chr.filter.vcf file is empty "
            exit 1
        fi    
        cat $input/$sample.variants.chr$chr.filter.vcf | awk '(length($4) == 1 && length($5) == 1 && $7 ~ /PASS/) || $0 ~ /#/' > $input/$sample.variants.chr$chr.SNV.filter.vcf
        cat $input/$sample.variants.chr$chr.filter.vcf | awk '(length($4) > 1 || length($5) > 1 && $7 ~ /PASS/) || $0 ~ /#/' > $input/$sample.variants.chr$chr.INDEL.filter.vcf
        
        intersect_file=$TargetKit
        $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$sample.chr$chr.SNV.bed.i.vcf 
        $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf 
        

        # convert to bed format and getting the coding variants or variants in target region
# cat $input/$sample.variants.chr$chr.SNV.filter.vcf | sed -e "/^#/d" | awk '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > $OnTarget/$sample.chr$chr.SNV.bed
#       cat $input/$sample.variants.chr$chr.INDEL.filter.vcf | sed -e "/^#/d" | awk '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > $OnTarget/$sample.chr$chr.INDEL.bed
        
        
		#if [ $tool == "exome" ]
        #then
         #   intersect_file=$TargetKit
        #elif [ $tool == "whole_genome" ]
        #then
        #    intersect_file=$master_gene_file
        #fi    
        
        ##intersecting the variant file using the interesect file
# $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.SNV.bed -b $intersect_file |sort |uniq | cut -f 1,3,4,5,6,7,8,9,10,11 > $OnTarget/$sample.chr$chr.SNV.bed.i
#       $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.INDEL.bed -b $intersect_file | sort | uniq | cut -f 1,3,4,5,6,7,8,9,10,11 > $OnTarget/$sample.chr$chr.INDEL.bed.i
#       rm $OnTarget/$sample.chr$chr.SNV.bed $OnTarget/$sample.chr$chr.INDEL.bed
        
         ## make bed files as VCF format
#       echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample}" > $OnTarget/$sample.chr$chr.header
#       cat $OnTarget/$sample.chr$chr.header $OnTarget/$sample.chr$chr.SNV.bed.i > $OnTarget/$sample.chr$chr.SNV.bed.i.vcf
#       cat $OnTarget/$sample.chr$chr.header $OnTarget/$sample.chr$chr.INDEL.bed.i > $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf
#       rm $OnTarget/$sample.chr$chr.SNV.bed.i $OnTarget/$sample.chr$chr.INDEL.bed.i
 
        perl $script_path/parse.vcf.INDEL.pl -i $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf -o $OnTarget/$sample.chr$chr.indels -s $sample
        perl $script_path/parse.vcf.SNV.pl -i $OnTarget/$sample.chr$chr.SNV.bed.i.vcf -o $OnTarget/$sample.chr$chr.snvs -s $sample
        #rm $OnTarget/$sample.chr$chr.header
        rm $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf
        rm $OnTarget/$sample.chr$chr.SNV.bed.i.vcf
        rm $input/$sample.variants.chr$chr.SNV.filter.vcf
        rm $input/$sample.variants.chr$chr.INDEL.filter.vcf
        ### interesect with capture kit to see if the vaariant is found in capture kit and annotate teh variant with 1 or 0 and for whole genome just 1
        if [ $tool == "exome" ]
        then
            ### SNVs
            awk '{print $1"\t"$2-1"\t"$2"\t"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9}' $OnTarget/$sample.chr$chr.snvs > $OnTarget/$sample.chr$chr.snvs.bed
            rm $OnTarget/$sample.chr$chr.snvs
            $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.snvs.bed -b $CaptureKit -c > $OnTarget/$sample.chr$chr.snvs.bed.i
            rm $OnTarget/$sample.chr$chr.snvs.bed
            perl -an -e '@a=split(/_/,$F[-2]); print join("\t",@a); print "\t$F[$#F]\n";' $OnTarget/$sample.chr$chr.snvs.bed.i > $OnTarget/$sample.chr$chr.raw.snvs.bed.i.ToMerge
            rm $OnTarget/$sample.chr$chr.snvs.bed.i 
            
            ### INDELs
            #perl -an -e 'if($F[3] =~ /^\+/){$F[2]=$F[2]+1;print join ("\t",@F),"\n"}else { print join ("\t",@F),"\n";}' $OnTarget/$sample.chr$chr.indels > $OnTarget/$sample.chr$chr.indels.bed
            #rm $OnTarget/$sample.chr$chr.indels 
            awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' $OnTarget/$sample.chr$chr.indels > $OnTarget/$sample.chr$chr.indels.bed
            rm $OnTarget/$sample.chr$chr.indels 
            $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.indels.bed -b $CaptureKit -c > $OnTarget/$sample.chr$chr.indels.bed.i
            rm $OnTarget/$sample.chr$chr.indels.bed
            #perl -an -e 'if($F[-2] =~ /^\+/){$F[-3]=$F[-4];print "$F[-5]\t$F[-4]\t$F[-3]\t$F[-2]\t$F[$#F]\n";}else { print "$F[-5]\t$F[-4]\t$F[-3]\t$F[-2]\t$F[$#F]\n";}' $OnTarget/$sample.chr$chr.indels.bed.i > $OnTarget/$sample.chr$chr.raw.indels.bed.i.ToMerge 
            perl -an -e '@a=split(/_/,$F[-2]); print join("\t",@a); print "\t$F[$#F]\n";' $OnTarget/$sample.chr$chr.indels.bed.i > $OnTarget/$sample.chr$chr.raw.indels.bed.i.ToMerge  
            rm $OnTarget/$sample.chr$chr.indels.bed.i
        elif [ $tool == "whole_genome" ]
        then
            cat $OnTarget/$sample.chr$chr.snvs | awk '{print $0"\t1"}' > $OnTarget/$sample.chr$chr.raw.snvs.bed.i.ToMerge
            rm $OnTarget/$sample.chr$chr.snvs
            cat $OnTarget/$sample.chr$chr.indels | awk '{print $0"\t1"}' > $OnTarget/$sample.chr$chr.raw.indels.bed.i.ToMerge
            rm $OnTarget/$sample.chr$chr.indels
        fi
        
        if [ ! -s $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge ]
        then
            #touch $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge
            echo "WARNING: $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge is empty "
        fi    
        ### intersect with indel call to flag the snv calls 0/1
        perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge -i $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge -n 10 -p 2 -o $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos
        mv $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge
        
        if [ ! -s $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge ]
        then
            echo "WARNING: $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge is empty"
        fi    
        
        
############# files are in the correct format required for TREAT annotation module to work
    else
        echo "Multi sample"
        #### we need to generate three files for a group for SNVs and INDELs
        ### somatic tumor normal
        group=$sample
        samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)
        ##  variants.chr.filter.vcf    
        input=$variants/$group
        #intersect_file=$master_gene_file
        intersect_file=$TargetKit
		for sample in $samples
        do  
            col=`cat $input/$group.variants.chr$chr.filter.vcf | awk '$0 ~ /#CHROM/' | awk -v num=$sample -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == num) {print i} } }'`
            
            cat $input/$group.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '(length($4) == 1 && length($5) == 1 && $7 ~ /PASS/) || $0 ~ /#/' | grep -v -w "\./\." > $input/$sample.variants.chr$chr.SNV.filter.vcf
            cat $input/$group.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '((length($4) > 1 || length($5) > 1 ) && $7 ~ /PASS/) || $0 ~ /#/' | grep -v -w "\./\."  > $input/$sample.variants.chr$chr.INDEL.filter.vcf
            
            # convert to bed format and getting the coding variants or variants in target region
            cat $input/$sample.variants.chr$chr.SNV.filter.vcf | sed -e "/^#/d" | awk 'BEGIN {OFS="\t"} { print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > $OnTarget/$sample.chr$chr.SNV.bed
            cat $input/$sample.variants.chr$chr.INDEL.filter.vcf | sed -e "/^#/d" | awk 'BEGIN {OFS="\t"} { print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > $OnTarget/$sample.chr$chr.INDEL.bed
            $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.SNV.bed -b $intersect_file |sort |uniq | awk 'BEGIN {OFS="\t"} { print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $OnTarget/$sample.chr$chr.SNV.bed.i
            $bedtools/intersectBed -a $OnTarget/$sample.chr$chr.INDEL.bed -b $intersect_file | sort | uniq | awk 'BEGIN {OFS="\t"} { print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $OnTarget/$sample.chr$chr.INDEL.bed.i
            rm $OnTarget/$sample.chr$chr.SNV.bed $OnTarget/$sample.chr$chr.INDEL.bed
            
            ## make bed files as VCF format
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample}" > $OnTarget/$sample.chr$chr.header
            cat $OnTarget/$sample.chr$chr.header $OnTarget/$sample.chr$chr.SNV.bed.i > $OnTarget/$sample.chr$chr.SNV.bed.i.vcf
            cat $OnTarget/$sample.chr$chr.header $OnTarget/$sample.chr$chr.INDEL.bed.i > $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf
            rm $OnTarget/$sample.chr$chr.SNV.bed.i $OnTarget/$sample.chr$chr.INDEL.bed.i
            perl $script_path/parse.vcf.INDEL.pl -i $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf -o $OnTarget/$sample.chr$chr.indels -s $sample
            perl $script_path/parse.vcf.SNV.pl -i $OnTarget/$sample.chr$chr.SNV.bed.i.vcf -o $OnTarget/$sample.chr$chr.snvs -s $sample
            rm $OnTarget/$sample.chr$chr.header
            rm $OnTarget/$sample.chr$chr.INDEL.bed.i.vcf
            rm $OnTarget/$sample.chr$chr.SNV.bed.i.vcf
            
            cat $OnTarget/$sample.chr$chr.snvs | awk '{print $0"\t1"}' > $OnTarget/$sample.chr$chr.raw.snvs.bed.i.ToMerge
            rm $OnTarget/$sample.chr$chr.snvs
            cat $OnTarget/$sample.chr$chr.indels | awk '{print $0"\t1"}' > $OnTarget/$sample.chr$chr.raw.indels.bed.i.ToMerge
            rm $OnTarget/$sample.chr$chr.indels
            
            if [ ! -s $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge ]
            then
                echo "WARNING: $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge is empty "
            fi    
            ### intersect with indel call to flag the snv calls 0/1
            perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge -i $OnTarget/$sample.chr${chr}.raw.indels.bed.i.ToMerge -n 10 -p 2 -o $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos
            mv $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge.pos $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge

            if [ ! -s $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge ]
            then
                echo "WARNING: $OnTarget/$sample.chr${chr}.raw.snvs.bed.i.ToMerge is empty"
            fi 
        done
        ### for somatic calls
        samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2 )
        i=1
        for sample in $samples
        do
            sampleArray[$i]=$sample
            let i=i+1
        done
        for i in $(seq 2 ${#sampleArray[@]})
        do  
            tumor=${sampleArray[$i]}
            col=`cat $input/$group.somatic.variants.chr$chr.filter.vcf | awk '$0 ~ /#CHROM/' | awk -v num=$tumor -F '\t' '{ for(i=1;i<=NF;i++){ if ($i == num) {print i} } }'`
            #cat $input/$group.somatic.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '(length($4) == 1 && length($5) == 1 && $7 ~ /PASS/) || $0 ~ /#/' | grep -v -w "0/0" > $input/$group.$tumor.somatic.variants.chr$chr.SNV.filter.vcf
            cat $input/$group.somatic.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '(length($4) == 1 && length($5)) || $0 ~ /#/' | grep -v -w "0/0" > $input/$group.$tumor.somatic.variants.chr$chr.SNV.filter.vcf
            #cat $input/$group.somatic.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '((length($4) > 1 || length($5) > 1) && $7 ~ /PASS/) || $0 ~ /#/' | grep -v -w "0/0"  > $input/$group.$tumor.somatic.variants.chr$chr.INDEL.filter.vcf
            cat $input/$group.somatic.variants.chr$chr.filter.vcf | awk '$0 !~ /##/' | awk -v num=$col 'BEGIN {OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$num } ' | awk '((length($4) > 1 || length($5) > 1)) || $0 ~ /#/' | grep -v -w "0/0"  > $input/$group.$tumor.somatic.variants.chr$chr.INDEL.filter.vcf
            # convert to bed format and getting the coding variants or variants in target region
            cat $input/$group.$tumor.somatic.variants.chr$chr.SNV.filter.vcf | sed -e "/^#/d" | awk 'BEGIN {OFS="\t"} { print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > $OnTarget/$group.$tumor.chr$chr.SNV.bed
            cat $input/$group.$tumor.somatic.variants.chr$chr.INDEL.filter.vcf| sed -e "/^#/d" | awk 'BEGIN {OFS="\t"} { print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > $OnTarget/$group.$tumor.chr$chr.INDEL.bed
            $bedtools/intersectBed -a $OnTarget/$group.$tumor.chr$chr.SNV.bed -b $intersect_file |sort |uniq | awk 'BEGIN {OFS="\t"} { print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $OnTarget/$group.$tumor.chr$chr.SNV.bed.i
            $bedtools/intersectBed -a $OnTarget/$group.$tumor.chr$chr.INDEL.bed -b $intersect_file | sort | uniq | awk 'BEGIN {OFS="\t"} { print $1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i
            rm $OnTarget/$group.$tumor.chr$chr.SNV.bed $OnTarget/$group.$tumor.chr$chr.INDEL.bed
            
            ## make bed files as VCF format
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${tumor}" > $OnTarget/$group.$tumor.chr$chr.header
            cat $OnTarget/$group.$tumor.chr$chr.header $OnTarget/$group.$tumor.chr$chr.SNV.bed.i > $OnTarget/$group.$tumor.chr$chr.SNV.bed.i.vcf
            cat $OnTarget/$group.$tumor.chr$chr.header $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i > $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i.vcf
            rm $OnTarget/$group.$tumor.chr$chr.SNV.bed.i $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i
            perl $script_path/parse.vcf.INDEL.pl -i $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i.vcf -o $OnTarget/$group.$tumor.chr$chr.indels -s $sample
            perl $script_path/parse.vcf.SNV.pl -i $OnTarget/$group.$tumor.chr$chr.SNV.bed.i.vcf -o $OnTarget/$group.$tumor.chr$chr.snvs -s $sample
            rm $OnTarget/$group.$tumor.chr$chr.header
            rm $OnTarget/$group.$tumor.chr$chr.INDEL.bed.i.vcf
            rm $OnTarget/$group.$tumor.chr$chr.SNV.bed.i.vcf

            cat $OnTarget/$group.$tumor.chr$chr.snvs | awk '{print $0"\t1"}' > $OnTarget/$group.$tumor.chr$chr.raw.snvs.bed.i.ToMerge
            rm $OnTarget/$group.$tumor.chr$chr.snvs
            cat $OnTarget/$group.$tumor.chr$chr.indels | awk '{print $0"\t1"}' > $OnTarget/$group.$tumor.chr$chr.raw.indels.bed.i.ToMerge
            rm $OnTarget/$group.$tumor.chr$chr.indels
            if [ ! -s $OnTarget/$group.$tumor.chr${chr}.raw.indels.bed.i.ToMerge ]
            then
                echo "WARNING: $OnTarget/$group.$tumor.chr${chr}.raw.indels.bed.i.ToMerge is empty "
            fi    
            ### intersect with indel call to flag the snv calls 0/1
            perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge -i $OnTarget/$group.$tumor.chr${chr}.raw.indels.bed.i.ToMerge -n 10 -p 2 -o $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge.pos
            mv $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge.pos $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge

            if [ ! -s $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge ]
            then
                echo "WARNING: $OnTarget/$group.$tumor.chr${chr}.raw.snvs.bed.i.ToMerge is empty"
            fi 
        done
        echo `date`
    fi
fi
