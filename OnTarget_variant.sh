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
    tool=`echo "$tool" | tr "[A-Z]" "[a-z]"`
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    distance=$( cat $tool_info | grep -w '^SNP_DISTANCE_INDEL' | cut -d '=' -f2 )
##############################################################		
    
    intersect_file=$TargetKit
    if [ $multi_sample != "YES" ]
    then
        echo "Single sample"
        input=$variants/$sample
        if [ ! -s $input/$sample.variants.chr$chr.filter.vcf ]
        then
            $script_path/errorlog.sh $input/$sample.variants.chr$chr.filter.vcf OnTarget_variants.sh ERROR "not exist"
            exit 1;
        fi    
        perl $script_path/vcf_to_variant_vcf.pl -i $input/$sample.variants.chr$chr.filter.vcf -v $input/$sample.variants.chr$chr.SNV.filter.vcf -l $input/$sample.variants.chr$chr.INDEL.filter.vcf -f 
        $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
        rm $input/$sample.variants.chr$chr.SNV.filter.vcf 
        $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
        rm $input/$sample.variants.chr$chr.INDEL.filter.vcf
        ### interesect with capture kit to see if the vaariant is found in capture kit and annotate teh variant with 1 or 0 and for whole genome just 1
        if [ $tool == "exome" ]
        then
            $bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
            rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
            cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf |  awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp
            mv $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
            $bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            rm $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
            cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.temp
            mv $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.temp  $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
        elif [ $tool == "whole_genome" ]
        then
            cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
            rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
            cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            rm $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
        fi
        if [ `cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
        then
            $script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			cp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL=0",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
        else
			perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
        fi
		rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
        if [ `cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
        then
            $script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
        fi      
        ### add format field tags to vcf
		perl $script_path/add_format_field_vcf.pl $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf SNV > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp
		mv $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
		perl $script_path/add_format_field_vcf.pl $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf INDEL > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
        mv $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
    else
        echo "Multi sample"
        #### we need to generate three files for a group for SNVs and INDELs
        ### somatic tumor normal
        group=$sample
        samples=$( cat $sample_info| grep -w "^$group" | cut -d '=' -f2)    
        input=$variants/$group
        intersect_file=$TargetKit
		for sample in $samples
        do  
            perl $script_path/vcf_to_variant_vcf.pl -i $input/$group.variants.chr$chr.filter.vcf -v $input/$sample.variants.chr$chr.SNV.filter.vcf -l $input/$sample.variants.chr$chr.INDEL.filter.vcf -s $sample -f
            cat $input/$sample.variants.chr$chr.SNV.filter.vcf | awk '$0 ~ /^#/ || ($10 !~ /^\./ && $10 !~ /^0\/0/)' > $input/$sample.variants.chr$chr.SNV.filter.vcf.temp
            mv $input/$sample.variants.chr$chr.SNV.filter.vcf.temp $input/$sample.variants.chr$chr.SNV.filter.vcf
            cat $input/$sample.variants.chr$chr.INDEL.filter.vcf | awk '$0 ~ /^#/ || ($10 !~ /^\./ && $10 !~ /^0\/0/)' > $input/$sample.variants.chr$chr.INDEL.filter.vcf.temp
            mv $input/$sample.variants.chr$chr.INDEL.filter.vcf.temp $input/$sample.variants.chr$chr.INDEL.filter.vcf
            $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
            rm $input/$sample.variants.chr$chr.SNV.filter.vcf 
            $bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
            rm $input/$sample.variants.chr$chr.INDEL.filter.vcf
            if [ $tool == "exome" ]
            then
                ### SNV
                $bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp
                mv $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.temp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
                $bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
                cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.temp
                mv $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.temp $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            else    
            cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
            cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
            fi
            rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
            rm $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
            if [ `cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
            then
                $script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
				cp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf 
				cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL="$NF,$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
			else
				perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
				cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL=0",$9,$10;}' > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
			fi
			rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
            if [ `cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
            then
				$script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
            fi
            perl $script_path/add_format_field_vcf.pl $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf SNV > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp
			mv $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf.tmp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
			perl $script_path/add_format_field_vcf.pl $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf INDEL > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
			mv $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf.tmp $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf	
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
            perl $script_path/vcf_to_variant_vcf.pl -i $input/$group.somatic.variants.chr$chr.filter.vcf -v $input/$group.$tumor.variants.chr$chr.SNV.filter.vcf -l $input/$group.$tumor.variants.chr$chr.INDEL.filter.vcf -s $tumor -f 
            $bedtools/intersectBed -header -a $input/$group.$tumor.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.vcf
            rm $input/$group.$tumor.variants.chr$chr.SNV.filter.vcf
            $bedtools/intersectBed -header -a $input/$group.$tumor.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.vcf
            rm $input/$group.$tumor.variants.chr$chr.INDEL.filter.vcf
            if [ $tool == "exome" ]
            then
                ##SNV
                $bedtools/intersectBed -header -a $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
                cat $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf.temp
                mv $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf.temp $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
                ##INDEL
                $bedtools/intersectBed -header -a $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c > $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf
                cat $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE="$NF,$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf.temp
                mv $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf.temp $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf 
            else    
                cat $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
                cat $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.vcf| awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CAPTURE=1",$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf
            fi
            rm $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.vcf
            rm $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.vcf
            if [ `cat $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
            then
                $script_path/errorlog.sh $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
				cp $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.pos.vcf
				cat $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL=0",$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
			else   
				perl $script_path/markSnv_IndelnPos.pl -s $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.pos.vcf
                cat $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.pos.vcf | awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8";CLOSE2INDEL="$NF,$9,$10;}' > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
            fi
			rm $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.pos.vcf
            if [ `cat $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
            then
                $script_path/errorlog.sh $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
            fi   
			perl $script_path/add_format_field_vcf.pl $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf SNV > $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf.tmp
			mv $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf.tmp $OnTarget/$group.$tumor.variants.chr$chr.SNV.filter.i.c.vcf
			perl $script_path/add_format_field_vcf.pl $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf INDEL > $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf.tmp
			mv $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf.tmp $OnTarget/$group.$tumor.variants.chr$chr.INDEL.filter.i.c.vcf	
        done
    fi
    echo `date`
fi
