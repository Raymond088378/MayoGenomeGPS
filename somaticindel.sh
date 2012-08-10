#!/bin/sh

if [ $# != 7 ]
then
    echo -e "Usage: script to run somtic indel caller \n <normal bam ><tumor bam><chromosome><tumor sample><output dir><output file vcf><run info>"
else
    set -x
    echo `date`
    tumor_bam=$1
    normal_bam=$2
    chr=$3
    tumor_sample=$4
    output=$5
    output_file=$6
    run_info=$7
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    gatk=$( cat $tool_info | grep -w '^GATK' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    window=$( cat $tool_info | grep -w '^INDEL_WINDOW_SIZE' | cut -d '=' -f2)
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	command_line_params=$( cat $tool_info | grep -w '^SOMATIC_INDEL_params' | cut -d '=' -f2 )
	export JAVA_HOME=$javahome
	export PATH=$javahome/bin:$PATH
	

    indel_v=$tumor_sample.chr$chr.indel.txt
	
    check=0
	count=0
	if [ ! -d $output/temp ]
	then
		mkdir $output/temp
	fi
	
	while [[ $check -eq 0 && $count -le 10 ]]
    do
		$java/java -Xmx3g -Xms512m -Djava.io.tmpdir=$output/temp/ -jar $gatk/GenomeAnalysisTK.jar \
		-R $ref \
		-et NO_ET \
		-K $gatk/Hossain.Asif_mayo.edu.key \
		-T SomaticIndelDetector \
		-L chr$chr \
		--window_size $window \
		-o $output/$output_file \
		-verbose $output/$indel_v \
		-I:normal $normal_bam \
		-I:tumor $tumor_bam $command_line_params
		sleep 1m
		check=`[ -s $output/$output_file.idx ] && echo "1" || echo "0"`
        if [ $check -eq 0 ]
        then
            rm `grep -l $output/$output_file *.log`
			rm core.*
        fi
		let count=count+1	
    done 
	
    if [ ! -s $output/$output_file.idx ]
    then
        $script_path/errorlog.sh $output/$output_file somaticindel.sh ERROR "failed to create"
		exit 1;
    else
        perl $script_path/convertvcf.pl $output/$output_file > $output/$output_file.tmp
		mv $output/$output_file.tmp $output/$output_file
		perl $script_path/fixindelAD.pl $output/$output_file $output/$output_file.temp
        mv $output/$output_file.temp $output/$output_file
        rm $output/$indel_v
        n=`cat $output/$output_file |  awk '$0 ~ /^##FORMAT=<ID=GT/' | wc -l`
		if [ $n == 0 ]
        then
            perl $script_path/add.gt.to.vcf.pl $output/$output_file | awk '$0 ~ /^#/ || $8 ~ /SOMATIC/' > $output/$output_file.tmp
        else
            cat $output/$output_file | awk '$0 ~ /^#/ || $8 ~ /SOMATIC/' > $output/$output_file.tmp
		fi
        cat $output/$output_file.tmp | awk '$0 ~ /^#/ || $5 ~ /,/' > $output/$output_file.multi.vcf
		cat $output/$output_file.tmp | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/$output_file
		rm $output/$output_file.tmp
        rm $output/$output_file.idx
    fi
    echo `date`
fi
