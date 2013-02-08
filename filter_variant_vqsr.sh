#!/bin/bash

### called from merge_variant_group and merge_variant_single

if [ $# -le  3 ]
then
    echo -e "Script to filter the variants using VQSR\nUsage: ./filter_variant_vqsr.sh <raw vcf complete path><outputfile><type={BOTH,SNP,INDEL}><run info>"
else
    set -x
    echo `date`
    inputvcf=$1
    outfile=$2
    variant=$3
    run_info=$4
    
    if [ $5 ]
    then
        somatic=$5
    fi
    
    input_dir=`dirname $inputvcf`
    input_name=`basename $inputvcf`
    output=`dirname $outfile`
	
    raw_calls=`cat $inputvcf | awk '$0 !~ /#/' | wc -l`
	
    if [ -d $output/temp ]
    then
        echo "temp already there"
    else
        mkdir -p $output/temp
    fi
    
    if [ -d $output/plot ]
    then
        echo "plot already there"
    else
        mkdir -p $output/plot
    fi

    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    samtools=$( cat $tool_info | grep -w '^SAMTOOLS' | cut -d '=' -f2)	
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    Kgenome=$( cat $tool_info | grep -w '^KGENOME_REF' | cut -d '=' -f2)
    mills=$( cat $tool_info | grep -w '^MILLS_REF' | cut -d '=' -f2)
    hapmap=$( cat $tool_info | grep -w '^HAPMAP_VCF' | cut -d '=' -f2)
    omni=$( cat $tool_info | grep -w '^OMNI_VCF' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    picard=$( cat $tool_info | grep -w '^PICARD' | cut -d '=' -f2 ) 
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2 ) 
    out=$( cat $run_info | grep -w '^BASE_OUTPUT_DIR' | cut -d '=' -f2)
    PI=$( cat $run_info | grep -w '^PI' | cut -d '=' -f2)
    tool=$( cat $run_info | grep -w '^TYPE' | cut -d '=' -f2 | tr "[A-Z]" "[a-z]" )
    run_num=$( cat $run_info | grep -w '^OUTPUT_FOLDER' | cut -d '=' -f2)
    r_soft=$( cat $tool_info | grep -w '^R_SOFT' | cut -d '=' -f2)
	vqsr_snv_param=$( cat $tool_info | grep -w '^VQSR_params_SNV' | cut -d '=' -f2 )
	vqsr_indel_param=$( cat $tool_info | grep -w '^VQSR_params_INDEL' | cut -d '=' -f2 )
	memory_info=$( cat $run_info | grep -w '^MEMORY_INFO' | cut -d '=' -f2)
	export PATH=$r_soft:$java:$PATH
	### split the vcf file for indels and snvs to appy the vqsr seperately
    $script_path/vcf_to_variant_vcf.pl -i $inputvcf -v $inputvcf.SNV.vcf -l $inputvcf.INDEL.vcf -t both   

	raw_calls_snvs=`cat $inputvcf.SNV.vcf | awk '$0 !~ /#/' | wc -l`
	raw_calls_indels=`cat $inputvcf.INDEL.vcf | awk '$0 !~ /#/' | wc -l`
	
    if [ ! -s $hapmap ]
    then
        $script_path/errorlog.sh $hapmap filter_variant_vqsr.sh ERROR "not exist"
        exit 1;
    fi

    if [ ! -s $omni ]
    then
        $script_path/errorlog.sh $omni filter_variant_vqsr.sh ERROR "not exist"
        exit 1;
    fi

    if [ ! -s $dbSNP ]
    then
        $script_path/errorlog.sh $dbSNP filter_variant_vqsr.sh ERROR "not exist"
        exit 1;
    fi

    if [ ! -s $inputvcf ]
    then
        $script_path/errorlog.sh $inputvcf filter_variant_vqsr.sh ERROR "not exist"
        exit 1;
    fi
	
	if [ ! -s $mills ]
	then
		$script_path/errorlog.sh $mills filter_variant_vqsr.sh ERROR "not exist"
		exit 1;
	fi	
    
    if [ $raw_calls_snvs -gt 0 ]
    then
		## Variant Recalibrator SNP
		if [ ! $5 ]
		then
			let check=0
			let count=0
			while [[ $check -eq 0 && $count -le 3 ]]
			do
                mem=$( cat $memory_info | grep -w '^VariantRecalibrator_JVM' | cut -d '=' -f2)
				$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-K $gatk/Hossain.Asif_mayo.edu.key \
				-T VariantRecalibrator \
				-mode SNP  \
				-input $inputvcf.SNV.vcf \
				-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
				-resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
				-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbSNP \
				-recalFile $output/temp/$input_name.recal \
				-tranchesFile $output/temp/$input_name.tranches \
				-rscriptFile $output/plot/$input_name.plots.R $vqsr_snv_param
				sleep 5
				check=`[ -s $output/temp/$input_name.tranches.pdf ] && echo "1" || echo "0"`
				if [ $check -eq 0 ]
				then
					if [[  `find . -name '*.log'` ]]
					then
						if [ `grep -l $output/plot/$input_name.plots.R *.log` ]
						then
							rm `grep -l $output/plot/$input_name.plots.R *.log`
							rm core.*
						fi
					fi
				fi
				let count=count+1
			done

			## Apply Recalibrator
			if [[ `cat $output/temp/$input_name.recal | wc -l` -gt 2  && -s $output/temp/$input_name.tranches ]]
			then
				mem=$( cat $memory_info | grep -w '^ApplyRecalibration_JVM' | cut -d '=' -f2)
				$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-K $gatk/Hossain.Asif_mayo.edu.key \
				-mode SNP \
				-T ApplyRecalibration \
				-nt $threads \
				-input $inputvcf.SNV.vcf \
				-recalFile $output/temp/$input_name.recal \
				-tranchesFile $output/temp/$input_name.tranches \
				-o $outfile.SNV.vcf
				sleep 5
			fi
			if [ -s $outfile.SNV.vcf ]
			then
				filter_calls_snvs=`cat $outfile.SNV.vcf | awk '$0 !~ /#/' | wc -l`
			else
				let filter_calls_snvs=0
			fi	
			if [[ ! -s $outfile.SNV.vcf || ! -s ${outfile}.SNV.vcf.idx || $raw_calls_snv != $filter_calls_snv ]]
			then
				$script_path/errorlog.sh $outfile.SNV.vcf filter_variant_vqsr.sh WARNING "failed to appy VQSR"
				### V3 best practice from GATK
				$script_path/hardfilters_variants.sh $inputvcf.SNV.vcf $outfile.SNV.vcf $run_info SNP
			else
				rm $inputvcf.SNV.vcf $inputvcf.SNV.vcf.idx        
			fi
			
			## remove the file
			if [ -s $output/temp/$input_name.recal ]
			then
				rm $output/temp/$input_name.recal $output/temp/$input_name.recal.idx
			fi
			
			if [ -s $output/temp/$input_name.tranches ]
			then
				rm $output/temp/$input_name.tranches $output/temp/$input_name.tranches.pdf 	
			fi
		else
			$script_path/hardfilters_variants.sh $inputvcf.SNV.vcf $outfile.SNV.vcf $run_info SNP
		fi
	else
		cp $inputvcf.SNV.vcf $outfile.SNV.vcf
        rm $inputvcf.SNV.vcf
	fi
    
	### applt vqsr for indels
	
	if [ $raw_calls_indels -gt 0 ]
	then
		if [ ! $5 ]
		then
			let count=0
			let check=0
			while [[ $check -eq 0 && $count -le 3 ]]
			do
				mem=$( cat $memory_info | grep -w '^VariantRecalibrator_JVM' | cut -d '=' -f2)
				$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-K $gatk/Hossain.Asif_mayo.edu.key \
				-T VariantRecalibrator \
				-mode INDEL  \
				-input $inputvcf.INDEL.vcf \
				-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $mills \
				-recalFile $output/temp/$input_name.recal \
				-tranchesFile $output/temp/$input_name.tranches \
				-rscriptFile $output/plot/$input_name.plots.R $vqsr_indel_param
				sleep 5 
				check=`[ -s $output/temp/$input_name.tranches.pdf ] && echo "1" || echo "0"`
				if [ $check -eq 0 ]
				then
					if [[  `find . -name '*.log'` ]]
					then
						if [ `grep -l $output/plot/$input_name.plots.R *.log` ]
						then
							rm `grep -l $output/plot/$input_name.plots.R *.log`
							rm core.*
						fi
					fi
				fi
				let count=count+1
			done	
			## Apply Recalibrator
			if [[ `cat $output/temp/$input_name.recal | wc -l` -gt 2  && -s $output/temp/$input_name.tranches ]]
			then
				mem=$( cat $memory_info | grep -w '^ApplyRecalibration_JVM' | cut -d '=' -f2)
				$java/java $mem -jar $gatk/GenomeAnalysisTK.jar \
				-R $ref \
				-et NO_ET \
				-K $gatk/Hossain.Asif_mayo.edu.key \
				-mode INDEL \
				-T ApplyRecalibration \
				-nt $threads \
				-input $inputvcf.INDEL.vcf \
				-recalFile $output/temp/$input_name.recal \
				-tranchesFile $output/temp/$input_name.tranches \
				-o $outfile.INDEL.vcf
				sleep 5 
			fi
			if [ -s $outfile.INDEL.vcf ]
			then
				filter_calls_indels=`cat $outfile.INDEL.vcf | awk '$0 !~ /#/' | wc -l`
			else
				let filter_calls_indels=0
			fi
			if [[ ! -s $outfile.INDEL.vcf || ! -s ${outfile}.INDEL.vcf.idx || $raw_calls_indels != $filter_calls_indels ]]
			then
				$script_path/errorlog.sh $outfile.INDEL.vcf filter_variant_vqsr.sh WARNING "failed to create"
				### V3 best practice from GATK
				$script_path/hardfilters_variants.sh $inputvcf.INDEL.vcf $outfile.INDEL.vcf $run_info INDEL
			else
				rm $inputvcf.INDEL.vcf $inputvcf.INDEL.vcf.idx        
			fi	
			## remove the file
			if [ -s $output/temp/$input_name.recal ]
			then
				rm $output/temp/$input_name.recal $output/temp/$input_name.recal.idx
			fi
			
			if [ -s $output/temp/$input_name.tranches ]
			then
				rm $output/temp/$input_name.tranches $output/temp/$input_name.tranches.pdf	
			fi
		else
			$script_path/hardfilters_variants.sh $inputvcf.INDEL.vcf $outfile.INDEL.vcf $run_info INDEL
		fi    
    else
        cp $inputvcf.INDEL.vcf $outfile.INDEL.vcf
        rm $inputvcf.INDEL.vcf
    fi    
	### combine both the indels and SNVs
	#in="-V $outfile.INDEL.vcf -V $outfile.SNV.vcf"
	in="$outfile.INDEL.vcf $outfile.SNV.vcf"
	$script_path/concatvcf.sh "$in" $outfile $run_info YES
	if [ ! -s $outfile ]
	then
		$script_path/errorlog.sh $outfile filter_variant_vqsr.sh ERROR "failed to create"
		exit 1;
	fi	
	echo `date`
fi	