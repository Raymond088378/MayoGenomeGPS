#!/bin/sh
if [ $# != 8 ]
then
    echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name> <output vcf file name> <run info>"
else
    set -x
    echo `date`
    normal_bam=$1
    tumor_bam=$2
    output=$3
    chr=$4
    tumor_sample=$5
    normal_sample=$6
    output_file=$7
    run_info=$8
    
    tool_info=$( cat $run_info | grep -w '^TOOL_INFO' | cut -d '=' -f2)
    mutect=$( cat $tool_info | grep -w '^MUTECT' | cut -d '=' -f2)
    java=$( cat $tool_info | grep -w '^JAVA' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2)
    dbSNP=$( cat $tool_info | grep -w '^dbSNP_REF' | cut -d '=' -f2)
    threads=$( cat $tool_info | grep -w '^THREADS' | cut -d '=' -f2)
    TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
    script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
    only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	javahome=$( cat $tool_info | grep -w '^JAVA_HOME' | cut -d '=' -f2 )
	
	export JAVA_HOME=$javahome
	export PATH=$javahome/bin:$PATH
    
    if [ $only_ontarget == "YES" ]
    then
	cat $TargetKit | grep -w chr$chr > $output/$tumor_sample.chr$chr.target.bed
        len=`cat $output/$tumor_sample.chr$chr.target.bed |wc -l`
        if [ $len -gt 0 ]
        then
            param="-L $output/$tumor_sample.chr$chr.target.bed"
        else
            param="-L chr$chr"
        fi    
    else
        param="-L chr$chr"
    fi
    
    $java/java -Xmx3g -Xms512m -jar $mutect/muTect-1.0.27783.jar \
    -T MuTect \
    --reference_sequence $ref \
    $param \
    --input_file:tumor $tumor_bam \
    --input_file:normal $normal_bam \
    -B:dbsnp,VCF $dbSNP \
    -et NO_ET \
    -nt $threads \
    --out $output/$output_file \
    --coverage_file $output/$tumor_sample.chr$chr.coverage.wig.txt
    
	perl $script_path/mutect2vcf.pl -i $output/$output_file -o $output/$output_file.temp -ns $normal_sample -ts $tumor_sample
	mv $output/$output_file.temp  $output/$output_file
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 ~ /,/' > $output/$output_file.multi.vcf
	cat $output/$output_file | awk '$0 ~ /^#/ || $5 !~ /,/' > $output/$output_file.tmp
	mv $output/$output_file.tmp  $output/$output_file
	
	if [ ! -s  $output/$output_file ]
	then
		$script_path/errorlog.sh $output/$output_file mutect.sh ERROR "failed to create"
		exit 1;
	fi
	if [ $only_ontarget == "YES" ]
	then
		rm $output/$tumor_sample.chr$chr.target.bed
	fi
	echo `date`
fi    