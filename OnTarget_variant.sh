#!/bin/bash
	
########################################################
###### 	SNV ANNOTATION FOR TUMOR/NORMAL PAIR WHOLE GENOME ANALYSIS PIPELINE
######		Program:			annotation.SNV.sh
######		Date:				11/09/2011
######		Summary:			Annotates GATK SNV and INDEL outputs
######		Input 
######		$1	=	structural directory
######		$2	=	/path/to/run_info.txt
########################################################

if [ $# -le 3 ]
then
    echo -e "script to get ontarget varaints uisng the capture kit specified\nUsage: </path/to/output directory for variants> </path/to/output directory for OnTarget> <sample name> </path/to/run_info.txt><SGE_TASK_ID(optional)>";
else
    set -x
    echo `date`
    variants=$1
    OnTarget=$2
    sample=$3
    run_info=$4
    if [ $5 ]
    then
        SGE_TASK_ID=$5
    fi
    ########################################################	
    ######		Reading run_info.txt and assigning to variables

    input=$( cat $run_info | grep -w '^INPUT_DIR' | cut -d '=' -f2)
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    bedtools=$( cat $tool_info | grep -w '^BEDTOOLS' | cut -d '=' -f2 )
    samples=$( cat $run_info | grep -w '^SAMPLENAMES' | cut -d '=' -f2)
    groups=$( cat $run_info | grep -w '^GROUPNAMES' | cut -d '=' -f2)
    multi_sample=$( cat $run_info | grep -w '^MULTISAMPLE' | cut -d '=' -f2)
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    chrs=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n")
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]")
    sample_info=$( cat $run_info | grep -w '^SAMPLE_INFO' | cut -d '=' -f2)
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    CaptureKit=$( cat $tool_info | grep -w '^CAPTUREKIT' | cut -d '=' -f2 )
    chr=$(cat $run_info | grep -w '^CHRINDEX' | cut -d '=' -f2 | tr ":" "\n" | head -n $SGE_TASK_ID | tail -n 1)
    distance=$( cat $tool_info | grep -w '^SNP_DISTANCE_INDEL' | cut -d '=' -f2 )
	gene_body=$( cat $tool_info | grep -w '^MATER_GENE_BODY' | cut -d '=' -f2 )
	
##############################################################		
    
    if [ $tool == "exome" ]
    then
		intersect_file=$TargetKit
    else
		intersect_file=$gene_body
    fi	
    
    if [ $multi_sample != "YES" ]
    then
		echo "Single sample"
		input=$variants/$sample
		if [ ! -s $input/$sample.variants.chr$chr.filter.vcf ]
		then
			$script_path/email.sh $input/$sample.variants.chr$chr.filter.vcf "doesn't exist" $run_info
			touch $input/$sample.variants.chr$chr.filter.vcf.fix.log
			$script_path/wait.sh $input/$sample.variants.chr$chr.filter.vcf.fix.log
		fi    
		$script_path/vcf_to_variant_vcf.pl -i $input/$sample.variants.chr$chr.filter.vcf -v $input/$sample.variants.chr$chr.SNV.filter.vcf -l $input/$sample.variants.chr$chr.INDEL.filter.vcf -t both -f 
		$bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
		rm $input/$sample.variants.chr$chr.SNV.filter.vcf 
		$bedtools/intersectBed -header -a $input/$sample.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
		rm $input/$sample.variants.chr$chr.INDEL.filter.vcf
		### interesect with capture kit to see if the vaariant is found in capture kit and annotate teh variant with 1 or 0 and for whole genome just 1
		if [ $tool == "exome" ]
		then
			$bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c |  $script_path/add.info.capture.vcf.pl  > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
			rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
			$bedtools/intersectBed -header -a $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c |  $script_path/add.info.capture.vcf.pl  > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
			rm $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
		elif [ $tool == "whole_genome" ]
		then
			cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf |  $script_path/add.info.capture.vcf.pl > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
			cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf |  $script_path/add.info.capture.vcf.pl > $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf
			rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.vcf
			rm $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.vcf
		fi
		if [ `cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf | wc -l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
			cp $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
		else
			$script_path/markSnv_IndelnPos.pl -s $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/$sample.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf
		fi
		rm $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.pos.vcf
		if [ `cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf |wc-l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/$sample.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
		fi      
    else
		echo "Multi sample"
		### multi sample calling
		group=$sample 
		input=$variants/$group 
		if [ ! -s $input/$group.variants.chr$chr.filter.vcf ]
		then
			$script_path/email.sh $input/$group.variants.chr$chr.filter.vcf "doesn't exist" $run_info
			touch $input/$group.variants.chr$chr.filter.vcf.fix.log
			$script_path/wait.sh $input/$group.variants.chr$chr.filter.vcf.fix.log
		fi   
		$script_path/vcf_to_variant_vcf.pl -i $input/$group.variants.chr$chr.filter.vcf -v $input/$group.variants.chr$chr.SNV.filter.vcf -l $input/$group.variants.chr$chr.INDEL.filter.vcf -t both -f
		$bedtools/intersectBed -header -a $input/$group.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/$group.variants.chr$chr.SNV.filter.i.vcf
		rm $input/$group.variants.chr$chr.SNV.filter.vcf 
		$bedtools/intersectBed -header -a $input/$group.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/$group.variants.chr$chr.INDEL.filter.i.vcf
		rm $input/$group.variants.chr$chr.INDEL.filter.vcf
		if [ $tool == "exome" ]
		then
			$bedtools/intersectBed -header -a $OnTarget/$group.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c |  $script_path/add.info.capture.vcf.pl  > $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf
			$bedtools/intersectBed -header -a $OnTarget/$group.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c | $script_path/add.info.capture.vcf.pl  > $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf
		else    
			cat $OnTarget/$group.variants.chr$chr.SNV.filter.i.vcf |  $script_path/add.info.capture.vcf.pl > $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf
			cat $OnTarget/$group.variants.chr$chr.INDEL.filter.i.vcf |  $script_path/add.info.capture.vcf.pl> $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf
		fi
		rm $OnTarget/$group.variants.chr$chr.SNV.filter.i.vcf
		rm $OnTarget/$group.variants.chr$chr.INDEL.filter.i.vcf
		if [ `cat $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf | wc -l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
			cp $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.pos.vcf 
			cat $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf
		else
			$script_path/markSnv_IndelnPos.pl -s $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/$group.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf
		fi
		rm $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.pos.vcf
		if [ `cat $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf | wc -l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/$group.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
		fi	
		
		### for somatic calls
		if [ ! -s $input/$group.somatic.variants.chr$chr.filter.vcf ]
		then
			$script_path/email.sh $input/$group.somatic.variants.chr$chr.filter.vcf "doesn't exist" $run_info
			touch $input/$group.somatic.variants.chr$chr.filter.vcf.fix.log
			$script_path/wait.sh $input/$group.somatic.variants.chr$chr.filter.vcf.fix.log
		fi   
		$script_path/vcf_to_variant_vcf.pl -i $input/$group.somatic.variants.chr$chr.filter.vcf -v $input/TUMOR.$group.variants.chr$chr.SNV.filter.vcf -l $input/TUMOR.$group.variants.chr$chr.INDEL.filter.vcf -t both -f
		$bedtools/intersectBed -header -a $input/TUMOR.$group.variants.chr$chr.SNV.filter.vcf -b $intersect_file > $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.vcf
		rm $input/TUMOR.$group.variants.chr$chr.SNV.filter.vcf
		$bedtools/intersectBed -header -a $input/TUMOR.$group.variants.chr$chr.INDEL.filter.vcf -b $intersect_file > $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.vcf
		rm $input/TUMOR.$group.variants.chr$chr.INDEL.filter.vcf
		if [ $tool == "exome" ]
		then
			$bedtools/intersectBed -header -a $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.vcf -b $CaptureKit -c |  $script_path/add.info.capture.vcf.pl  > $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf
			$bedtools/intersectBed -header -a $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.vcf -b $CaptureKit -c | $script_path/add.info.capture.vcf.pl  > $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf 
		else    
			cat $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.vcf | $script_path/add.info.capture.vcf.pl > $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf
			cat $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.vcf | $script_path/add.info.capture.vcf.pl > $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf
		fi
		rm $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.vcf
		rm $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.vcf
		if [ `cat $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf | wc -l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
			cp $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf
		else   
			$script_path/markSnv_IndelnPos.pl -s $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf -i $OnTarget/TUMOR.$group.variants.chr$chr.INDEL.filter.i.c.vcf -n $distance -o $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.pos.vcf
			cat $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.pos.vcf | $script_path/add.info.close2indel.vcf.pl > $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf
		fi
		rm $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.pos.vcf
		if [ `cat $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf | awk '$0 !~ /^#/' | wc -l` -lt 1 ]
		then
			if [ `cat $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf | wc -l` -lt 1 ]
			then
				$script_path/errorlog.sh $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh ERROR "failed to generate"
				exit 1;
			else	
				$script_path/errorlog.sh $OnTarget/TUMOR.$group.variants.chr$chr.SNV.filter.i.c.vcf OnTarget_variants.sh WARNING "no variant calls"
			fi
		fi   
    fi
    echo `date`
fi
