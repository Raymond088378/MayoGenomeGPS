#!/bin/sh

if [ $# != 3 ]
then
    echo -e "Usage: to annotate vcf files \n <input vcf file><chromosome><run info file></path/to/tumor bam> </path/to/normal bam> "
else
    set -x
    echo `date`
    vcf=$1
    run_info=$2
    input=$3
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    script_path=$( cat $tool_info | grep -w '^WORKFLOW_PATH' | cut -d '=' -f2)
    params=$( cat $tool_info | grep -w '^VCF_annotation_params' | cut -d '=' -f2)
    out=`dirname $vcf`
    if [ ! -d $out/temp ]
    then
        mkdir -p $out/temp
    fi
    
    export PATH=$java:$PATH
    # ## annotate SNVs or INDELs
    $java/java -Xmx4g -Xms512m -Djava.io.tmpdir=$out/temp/ -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -K $gatk/Hossain.Asif_mayo.edu.key \
    -T VariantAnnotator \
    $input \
    -V $vcf \
    --dbsnp $dbSNP \
    -L $vcf --out $vcf.temp $params
    	
    if [ -s $vcf.temp.idx ]
    then    
        mv $vcf.temp.idx $vcf.idx
        mv $vcf.temp $vcf
    else	
        $script_path/errorlog.sh $vcf.temp.idx annotate_vcf.sh WARNING "failed to annotate"
        rm $vcf.temp $vcf.temp.idx
    fi
    echo `date`
fi	
