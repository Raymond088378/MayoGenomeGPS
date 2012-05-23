#!/bin/sh

if [ $# != 4 ]
then
    echo -e "Usage: to annotate vcf files \n <input vcf file><chromosome><run info file></path/to/tumor bam> </path/to/normal bam> "
else
    set -x
    echo `date`
    
    vcf=$1
    chr=$2
    run_info=$3
    input=$4
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2)
    export JAVA_HOME=$java

    # ## annotate SNVs
    $java/java -Xmx6g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -K $gatk/Hossain.Asif_mayo.edu.key \
    -T VariantAnnotator \
    $input \
    -V $vcf \
    --dbsnp $dbSNP \
    -L $vcf \
    -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A HaplotypeScore -A DepthOfCoverage -A MappingQualityZero -A RMSMappingQuality -A FisherStrand \
    --out $vcf.temp	
    
    if [ -s $vcf.temp ]
    then
        mv $vcf.temp $vcf
    else		
        $script_path/errorlog.sh $vcf.temp annotate_vcf.sh ERROR "failed to create"
        exit 1
    fi

    if [ -s $vcf.temp.idx ]
    then
        mv $vcf.temp.idx $vcf.idx
    else	
        $script_path/errorlog.sh $vcf.temp.idx annotate_vcf.sh ERROR "failed to create"
        exit 1
    fi
    echo `date`
fi	