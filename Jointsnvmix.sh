#!/bin/sh

if [ $# != 8 ]
then
    echo "Usage: <normal bam> <tumor bam > <output dir> <chromosome> <tumor sample name> <normal sample name ><output file> <run info>"
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
    jointsnvmix=$( cat $tool_info | grep -w '^JOINTSNVMIX' | cut -d '=' -f2)    
    python=$( cat $tool_info | grep -w '^PYTHON' | cut -d '=' -f2) 
    pythonpath=$( cat $tool_info | grep -w '^PYTHONLIB' | cut -d '=' -f2)
    ref=$( cat $tool_info | grep -w '^REF_GENOME' | cut -d '=' -f2) 
	script_path=$( cat $tool_info | grep -w '^WHOLEGENOME_PATH' | cut -d '=' -f2 )
	TargetKit=$( cat $tool_info | grep -w '^ONTARGET' | cut -d '=' -f2 )
	only_ontarget=$( cat $tool_info | grep -w '^TARGETTED' | cut -d '=' -f2 | tr "[a-z]" "[A-Z]" )
	mqual=$( cat $tool_info | grep -w 'MAPPING_QUALITY' | cut -d '=' -f2)
	bqual=$( cat $tool_info | grep -w 'BASE_QUALITY' | cut -d '=' -f2)
	
	if [ $only_ontarget == "YES" ]
	then
		cat $TargetKit | grep -w chr$chr > $output/$sample.$chr.target.bed
    fi
	
    export PYTHONPATH=$PYTHONPATH:$pythonpath
    export PATH=$PATH:$PYTHONPATH
    
    if [ ! -s $normal_bam ]
    then
        echo "$normal_bam normal bam doesn't exist"
        exit 1;
    fi
    
    if [ ! -s $tumor_bam ]
    then
        echo "$tumor_bam tumor bam doesn't exist"
        exit 1;
    fi
    
    ### make sure both the bams are sorted
	
   # $python $jointsnvmix/build/scripts-2.7/jsm.py train --model snvmix2  --skip_size $X --min_base_qual 20 --min_map_qual 20 --chromosome chr$chr \
	#	$ref $normal_bam $tumor_bam $output/$tumor_sample.$normal_sample.chr$chr.train.txt
		
		
    ### run joint snvmix classify to call teh somatic mutation
    
	$python $jointsnvmix/build/scripts-2.7/jsm.py classify --model snvmix2 --post_process --min_base_qual $bqual --min_map_qual $mqual --chromosome chr$chr --out_file $output/$output_file.txt --parameters_file $jointsnvmix/config/params.cfg $ref $normal_bam $tumor_bam
	
		
	### script to convert text output to vcf output
	perl $script_path/jsm2vcf.pl -i $output/$output_file.txt -o $output/$output_file -ns $normal_sample -ts $tumor_sample
	rm $output/$output_file.txt
    echo `date`
fi