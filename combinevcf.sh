#!/bin/sh

if [ $# != 4 ]
then
    echo -e "Usage: combine vcfs <input files><output vcf ><run info><to delete input files(yes/no)"
else
    set -x
    echo `date`
    
    input=$1
    output=$2
    run_info=$3
    flag=`echo $4 | tr "[a-z]" "[A-Z]"`
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    
    
    $java/java -Xmx2g -Xms512m -jar $gatk/GenomeAnalysisTK.jar \
    -R $ref \
    -et NO_ET \
    -K $gatk/Hossain.Asif_mayo.edu.key \
    -T CombineVariants \
    $input \
    -o $output
    
    if [ ! -s $output ]
    then
        echo "ERROR: combinevcf.sh failed to generate $output"
    else
        if [ $flag == "YES" ]
        then
            rm `echo $input | sed -e '/-V/s///g'`
            sample=`echo $input | sed -e '/-V/s///g' | tr " " "\n" | awk '{print $0".idx"}'`
			for i in $sample
			do
				if [ $i != ".idx" ]
				then
					rm $i
				fi
			done	
        fi
    fi
    echo `date`
fi	
	