#!/bin/sh

if [ $# != 4 ]
then
	echo -e "Script to appy hard filters on vcf variant files\nUsage: script to appy hard filters \n <input vcf > <output vcf ><run info file><type SNP/INDEL>"
else
    set -x
    echo `date`
    inputvcf=$1
    outputvcf=$2
    run_info=$3
    type=$4        
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2 )
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)	
    	
	export PATH=$java:$PATH
    
    if [ $type == "SNP" ]
	then
		num_snvs=`cat $inputvcf | awk '$0 !~ /^#/' | wc -l`
	
		if [ $num_snvs -ge 1 ]
		then
			### apply filters for SNV
			$java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
			-R $ref \
			-et NO_ET \
			-K $gatk/Hossain.Asif_mayo.edu.key \
			-l INFO \
			-T VariantFiltration \
			-V $inputvcf \
			-o $outputvcf \
			--filterExpression "MQ < 40.0" --filterName MQFilter \
			--filterExpression "FS > 60.0" --filterName FSFilter \
			--filterExpression "HaplotypeScore > 13.0" --filterName Haplotypefilter \
			--filterExpression "MQRankSum < -12.5" --filterName MappingQaulityFilter \
			--filterExpression "QD < 2.0" --filterName QDFilter \
			--filterExpression "ReadPosRankSum < -8.0" --filterName ReadPosFilter
				
			filter_snvs=`cat $outputvcf | awk '$0 !~ /^#/' | wc -l`
			if [[ -s $outputvcf.idx && $num_snvs == $filter_snvs ]]
			then
				rm $inputvcf $inputvcf.idx
			else
				$script_path/errorlog.sh $outputvcf.idx hardfilters_variants.sh ERROR "failed to create"
				exit 1;
			fi
		fi	
	fi
	
	if [ $type == "INDEL" ]
	then
		num_indels=`cat $inputvcf | awk '$0 !~ /^#/' | wc -l`
		if [ $num_indels -ge 1 ]
		then
			#### apply filters for INDEL
			$java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
			-R $ref \
			-et NO_ET \
			-K $gatk/Hossain.Asif_mayo.edu.key \
			-l INFO \
			-T VariantFiltration -V $inputvcf -o $outputvcf \
			--filterExpression "ReadPosRankSum < -20.0" --filterName ReadPosFilter \
			--filterExpression "QD < 2.0" --filterName QDFilter \
			--filterExpression "FS > 200.0" --filterName FSFilter
        
			filter_indels=`cat $outputvcf | awk '$0 !~ /^#/' | wc -l`
			if [[ -s $outputvcf.idx && $num_indels == $filter_indels ]]
			then
				rm $inputvcf $inputvcf.idx 
			else
				$script_path/errorlog.sh $outputvcf.idx hardfilters_variants.sh ERROR "failed to create"
				exit 1;
			fi
		fi
    fi
	echo `date`
fi	
	
	