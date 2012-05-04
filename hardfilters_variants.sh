#!/bin/sh

if [ $# != 3 ]
then
	echo "Usage: script to appy hard filters \n <input vcf > <output vcf ><run info file>"
else
    set -x
    echo `date`
    inputvcf=$1
    outputvcf=$2
    run_info=$3
            
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)	
    
    cat $inputvcf | awk '(length($4) == 1 && length($5) == 1) || $0 ~ /#/' > $inputvcf.SNV.vcf
    cat $inputvcf | awk '(length($4) > 1 || length($5) > 1) || $0 ~ /#/' > $inputvcf.INDEL.vcf
    
    num_snvs=`cat $inputvcf.SNV.vcf | awk '$0 !~ /^#/' | wc -l`
	
    if [ $num_snvs -ge 1 ]
    then
        ### apply filters for SNV
        $java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -K $gatk/Hossain.Asif_mayo.edu.key \
        -l INFO \
        -T VariantFiltration \
        -V $inputvcf.SNV.vcf \
        -o $inputvcf.SNV.vcf.filter.vcf \
        --filterExpression "QD < 2.0" \
        --filterName QDFilter \
        --filterExpression "MQ < 40.0" \
        --filterExpression "FS > 60.0" \
        --filterName FSFilter \
        --filterExpression "HaplotypeScore > 13.0" \
        --filterName Haplotypefilter \
        --filterExpression "MQRankSum < -12.5" \
        --filterName MappingQaulityFilter \
        --filterExpression "ReadPosRankSum < -8.0" \
        --filterName ReadPosFilter 
            
        filter_snvs=`cat $inputvcf.SNV.vcf.filter.vcf | awk '$0 !~ /^#/' | wc -l`
        if [[ -s $inputvcf.SNV.vcf.filter.vcf && $num_snvs == $filter_snvs ]]
        then
            rm $inputvcf.SNV.vcf
        else
            echo "ERROR: failed to appply hard filters $inputvcf (SNV)"
            exit 1;
        fi
    fi	
	
    num_indels=`cat $inputvcf.INDEL.vcf | awk '$0 !~ /^#/' | wc -l`
    if [ $num_indels -ge 1 ]
    then
        #### apply filters for INDEL
        $java/java -Xmx1g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
        -R $ref \
        -et NO_ET \
        -K $gatk/Hossain.Asif_mayo.edu.key \
        -l INFO \
        -T VariantFiltration \
        -V $inputvcf.INDEL.vcf \	
        -o $inputvcf.INDEL.vcf.filter.vcf \
        --filterExpression "QD < 2.0" \
        --filterName QDFilter \
        --filterExpression "ReadPosRankSum < -20.0" \
        --filterName ReadPosFilter \
        --filterExpression "FS > 200.0" \
        --filterName FSFilter 
		
        filter_indels=`cat $inputvcf.INDEL.vcf.filter.vcf | awk '$0 !~ /^#/' | wc -l`
        if [[ -s $inputvcf.INDEL.vcf.filter.vcf && $num_indels == $filter_indels ]]
        then
            rm $inputvcf.INDEL.vcf
        else
            echo "ERROR: failed to appply hard filters $inputvcf (INDEL)"
            exit 1;
        fi
    fi
	
    if [[ $num_snvs -gt 1 && $num_indels -gt 1 ]]	
    then
        input="-V $inputvcf.SNV.vcf.filter.vcf -V $inputvcf.INDEL.vcf.filter.vcf"
        $script_path/combinevcf.sh $input $outputvcf $run_info
    elif [ $num_snvs -gt 1 ]
    then
        mv $inputvcf.SNV.vcf.filter.vcf $outputvcf
    elif [ $num_indels -gt 1 ]
    then
        mv $inputvcf.INDEL.vcf.filter.vcf $outputvcf
    else
        "ERROR: no snvs or indels in input vcf file $inputvcf"
    fi
    
    if [ ! -s $outputvcf ]
    then
        echo "Merging the varaint file failed for $outputvcf"
        exit 1;
    fi		
    echo` date`
fi	
	
	